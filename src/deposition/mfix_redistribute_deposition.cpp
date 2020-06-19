#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLTensorOp.H>
#include <AMReX_MLEBTensorOp.H>

//
// Redistribute solids volume fraction
//
void
mfix::mfix_redistribute_deposition (int lev,
                                    amrex::MultiFab & mf_eps,
                                    amrex::MultiFab & mf_to_redistribute,
                                    const amrex::MultiFab * volfrac,
                                    const amrex::FabArray<EBCellFlagFab>* flags_fab,
                                    amrex::Real max_eps)
{
   BL_PROFILE("mfix::mfix_redistribute_solids_volume");

   MultiFab mf_to_redist_copy(mf_to_redistribute.boxArray(),
                              mf_to_redistribute.DistributionMap(),
                              mf_to_redistribute.nComp(),
                              mf_to_redistribute.nGrow(),
                              MFInfo(),
                              mf_to_redistribute.Factory());
   MultiFab::Copy(mf_to_redist_copy, mf_to_redistribute, 0, 0,
                  mf_to_redistribute.nComp(), mf_to_redistribute.nGrow());

   MultiFab scale_fab(mf_eps.boxArray(), mf_eps.DistributionMap(), mf_eps.nComp(),
                      mf_eps.nGrow(), MFInfo(), mf_eps.Factory());
   scale_fab.setVal(0.);

   for (MFIter mfi(mf_eps,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

     // We don't want to do this for ghost cells
     const Box& bx = mfi.tilebox();

     const int ncomp = mf_to_redistribute.nComp();

     // We are only interested in redistributing excessive particle volume
     // from small cells.

     if ( (*flags_fab)[mfi].getType(amrex::grow(bx,0)) != FabType::covered &&
          (*flags_fab)[mfi].getType(amrex::grow(bx,0)) != FabType::regular ) {

       // We don't want to loop over ghost cells. Only redistribute
       // solids volume that is locally owned.
       const Box& grow_bx1 = amrex::grow(bx,1);

       Box domain(geom[lev].Domain());
       const amrex::Dim3 dom_low  = amrex::lbound(domain);
       const amrex::Dim3 dom_high = amrex::ubound(domain);

       Array4<const EBCellFlag> const& flags = (*flags_fab)[mfi].array();
       Array4<const Real> const& vfrac = volfrac->array(mfi);

       Array4<Real> const&      ep_s = mf_eps.array(mfi);

       const int cyclic_x = geom[0].isPeriodic(0);
       const int cyclic_y = geom[0].isPeriodic(1);
       const int cyclic_z = geom[0].isPeriodic(2);

       IntVect mask_box_lo(grow_bx1.smallEnd());
       IntVect mask_box_hi(grow_bx1.bigEnd());

       if(not cyclic_x) {
         mask_box_lo[0] = amrex::max(mask_box_lo[0], dom_low.x);
         mask_box_hi[0] = amrex::min(mask_box_hi[0], dom_high.x);
       }

       if(not cyclic_y) {
         mask_box_lo[1] = amrex::max(mask_box_lo[1], dom_low.y);
         mask_box_hi[1] = amrex::min(mask_box_hi[1], dom_high.y);
       }

       if(not cyclic_z) {
         mask_box_lo[2] = amrex::max(mask_box_lo[2], dom_low.z);
         mask_box_hi[2] = amrex::min(mask_box_hi[2], dom_high.z);
       }

       // Box "mask_box" is used to restrict were we redistribute the overflow.
       // The following is what we want to do:
       // -- Mask ghost cells when the BCs are not periodic
       // -- Mask cells we are going to redistribute (ep_s > max_eps)
       Box mask_box(mask_box_lo, mask_box_hi);

       Array4<Real> const& mf_redist = mf_to_redistribute.array(mfi);

       Array4<Real> const& duplicate = mf_to_redist_copy.array(mfi);
       Array4<Real> const& scale_array = scale_fab.array(mfi);

       amrex::ParallelFor(bx,
         [flags,ep_s,mf_redist,vfrac,duplicate,scale_array,max_eps,ncomp,
          mask_box] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         if(flags(i,j,k).isSingleValued() and ep_s(i,j,k) > max_eps)
         {
           amrex::Real sum_vfrac_eps = 0.0;
           amrex::Real sum_vfrac     = 0.0;

           for(int ii(-1); ii <= 1; ii++)
           for(int jj(-1); jj <= 1; jj++)
           for(int kk(-1); kk <= 1; kk++) {
             if((ii != 0 or jj != 0 or kk != 0 ) and
                flags(i,j,k).isConnected({ii,jj,kk}) and
                mask_box.contains(IntVect(i+ii,j+jj,k+kk)) and
                ((not flags(i+ii,j+jj,k+kk).isSingleValued()) or
                 (ep_s(i+ii,j+jj,k+kk) <= max_eps)))
             {
               sum_vfrac     += vfrac(i+ii,j+jj,k+kk);
               sum_vfrac_eps += vfrac(i+ii,j+jj,k+kk)*ep_s(i+ii,j+jj,k+kk);
             }
           }

           // Average volume fraction in the neighborhood around of the packed
           // cell. This value is capped by the user-defined max pack value.
           amrex::Real avg_eps = amrex::min(max_eps, sum_vfrac_eps / sum_vfrac);

           // Fraction of material we want to redistribute
           amrex::Real scale = amrex::max(0.0, ep_s(i,j,k) - avg_eps) / ep_s(i,j,k);
           scale_array(i,j,k) = scale;

           for(int n(0); n < ncomp; n++) {

             // This is the amount of the multifab we are redistributing
             // that we are moving into the packed cell's neighborhood.
             amrex::Real overflow =
               duplicate(i,j,k,n) * scale * vfrac(i,j,k) / sum_vfrac;

             for(int kk(-1); kk <= 1; kk++)
             for(int jj(-1); jj <= 1; jj++)
             for(int ii(-1); ii <= 1; ii++) {
               if((ii != 0 or jj != 0 or kk != 0) and
                  flags(i,j,k).isConnected({ii,jj,kk}) and
                  mask_box.contains(IntVect(i+ii,j+jj,k+kk)) and
                  ((not flags(i+ii,j+jj,k+kk).isSingleValued()) or
                   (ep_s(i+ii,j+jj,k+kk) <= max_eps)))
               {
                 Gpu::Atomic::Add(&mf_redist(i+ii,j+jj,k+kk,n), overflow);
               }
             }
           }
        }
      });

      amrex::ParallelFor(bx,
        [flags,ep_s,mf_redist,scale_array,max_eps,ncomp]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if(flags(i,j,k).isSingleValued() and ep_s(i,j,k) > max_eps)
        {
          const Real scale = scale_array(i,j,k);

          for(int n(0); n < ncomp; n++) {
            Real redist_val = mf_redist(i,j,k,n);

            // Account for the change in material in the source cell
            amrex::Real delta = -scale*redist_val;
            redist_val += delta;
            mf_redist(i,j,k,n) = redist_val;
          }
        }
      });

    }
  }
}
