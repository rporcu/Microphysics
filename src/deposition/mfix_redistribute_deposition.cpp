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

       FArrayBox  mask_fbx(grow_bx1);
       Array4<Real> const& mask = mask_fbx.array();

       Box domain(geom[lev].Domain());
       const amrex::Dim3 dom_low  = amrex::lbound(domain);
       const amrex::Dim3 dom_high = amrex::ubound(domain);

       Array4<const EBCellFlag> const& flags = (*flags_fab)[mfi].array();
       Array4<const Real> const& vfrac = volfrac->array(mfi);

       Array4<Real> const&      ep_s = mf_eps.array(mfi);

       const int cyclic_x = geom[0].isPeriodic(0);
       const int cyclic_y = geom[0].isPeriodic(1);
       const int cyclic_z = geom[0].isPeriodic(2);

       // Array "mask" is used to restrict were we redistribute the overflow.
       // -- Mask ghost cells when the BCs are not periodic
       // -- Mask cells we are going to redistribute (ep_s > max_eps)
       amrex::ParallelFor(grow_bx1,
         [mask,flags,ep_s,cyclic_x,cyclic_y,cyclic_z,dom_low,dom_high,max_eps]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept 
         {
           if(((not cyclic_x) and (i < dom_low.x or i > dom_high.x)) or
              ((not cyclic_y) and (j < dom_low.y or j > dom_high.y)) or
              ((not cyclic_z) and (k < dom_low.z or k > dom_high.z)) or
              (flags(i,j,k).isSingleValued() and  (ep_s(i,j,k) > max_eps)))
             mask(i,j,k) = 0.0;
           else
             mask(i,j,k) = 1.0;
         });

       Array4<Real> const& mf_redist = mf_to_redistribute.array(mfi);

       amrex::ParallelFor(bx,
         [flags,ep_s,mf_redist,mask,vfrac,max_eps,ncomp]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {

           if(flags(i,j,k).isSingleValued()){

             if( ep_s(i,j,k) > max_eps){

               amrex::Real sum_vfrac_eps = 0.0;
               amrex::Real sum_vfrac     = 0.0;

               for(int ii(-1); ii <= 1; ii++)
                 for(int jj(-1); jj <= 1; jj++)
                   for(int kk(-1); kk <= 1; kk++)
                     if ( (ii != 0 or jj != 0 or kk != 0 ) and
                          (flags(i,j,k).isConnected({ii,jj,kk}) == 1)) {
                       sum_vfrac     += mask(i+ii,j+jj,k+kk)*vfrac(i+ii,j+jj,k+kk);
                       sum_vfrac_eps += mask(i+ii,j+jj,k+kk)*vfrac(i+ii,j+jj,k+kk)*ep_s(i+ii,j+jj,k+kk);
                     }

               // Average volume fraction in the neighborhood around of the packed
               // cell. This value is capped by the user-defined max pack value.
               amrex::Real avg_eps = amrex::min(max_eps, sum_vfrac_eps / sum_vfrac );

               // Fraction of material we want to redistribute
               amrex::Real scale = amrex::max(0.0, ep_s(i,j,k) - avg_eps) / ep_s(i,j,k);

               for(int n(0); n < ncomp; n++) {

                 // This is the amount of the multifab we are redistributing
                 // that we are moving into the packed cell's neighborhood.
                 amrex::Real overflow = mf_redist(i,j,k,n) * scale * vfrac(i,j,k) / sum_vfrac;

                 for(int kk(-1); kk <= 1; kk++) {
                   for(int jj(-1); jj <= 1; jj++) {
                     for(int ii(-1); ii <= 1; ii++) {
                       if((ii != 0 or jj != 0 or kk != 0) and
                          (flags(i,j,k).isConnected({ii,jj,kk}))) {

                         amrex::Gpu::Atomic::Add(&mf_redist(i+ii,j+jj,k+kk,n),
                                                 mask(i+ii,j+jj,k+kk)*overflow);
                       }
                     }
                   }
                 }

                 // Account for the change in material in the source cell
                 amrex::Real delta = -scale*mf_redist(i,j,k,n);
                 amrex::Gpu::Atomic::Add(&mf_redist(i,j,k,n), delta);
               }

             }
           }
         });

       amrex::Gpu::synchronize();
     }

   }
}
