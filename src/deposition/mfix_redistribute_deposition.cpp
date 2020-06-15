#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLTensorOp.H>
#include <AMReX_MLEBTensorOp.H>

namespace {

struct get_index
{
  const Long jstride;
  const Long kstride;
  const Long nstride;
  const Dim3 begin;
  const Dim3 end;
  const int ncomp;

  AMREX_GPU_HOST_DEVICE
  constexpr get_index(Dim3 const& a_begin, Dim3 const& a_end, int a_ncomp) noexcept
    : jstride(a_end.x-a_begin.x)
    , kstride(jstride*(a_end.y-a_begin.y))
    , nstride(kstride*(a_end.z-a_begin.z))
    , begin(a_begin)
    , end(a_end)
    , ncomp(a_ncomp)
  {}

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  int operator() (int i, int j, int k, int n=0) const noexcept
  {
    return (i-begin.x)+(j-begin.y)*jstride+(k-begin.z)*kstride+n*nstride;
  }
};

}

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

#ifdef AMREX_USE_CUDA
       dim3 Block(4,4,4);
       dim3 Grid(amrex::Math::ceil(bx.length()[0] / float(Block.x)),
                 amrex::Math::ceil(bx.length()[1] / float(Block.y)),
                 amrex::Math::ceil(bx.length()[2] / float(Block.z)));

       const unsigned sm_size =
         sizeof(Real)*(Block.x+2)*(Block.y+2)*(Block.z+2) + // eps_sm
         sizeof(Real)*(Block.x+2)*(Block.y+2)*(Block.z+2) + // vfrac_sm
         sizeof(EBCellFlag)*(Block.x+2)*(Block.y+2)*(Block.z+2); //flags_sm

       if(sm_size > Gpu::Device::sharedMemPerBlock()) {
         amrex::Abort("Exceeding GPU shared memory in"
             "mfix_redistribute_deposition_gpu.cpp");
       }
         
       amrex::IntVect bx_lo(bx.loVect());
       amrex::IntVect bx_hi(bx.hiVect());

       amrex::launch_global<<<Grid,Block,sm_size,Gpu::gpuStream()>>>(
         [flags,ep_s,mf_redist,vfrac,duplicate,scale_array,max_eps,
          ncomp,mask_box,bx_lo,bx_hi] AMREX_GPU_DEVICE () noexcept
       {
         int tid_x = threadIdx.x;
         int tid_y = threadIdx.y;
         int tid_z = threadIdx.z;

         const int stride_x = blockIdx.x*blockDim.x + bx_lo[0];
         const int stride_y = blockIdx.y*blockDim.y + bx_lo[1];
         const int stride_z = blockIdx.z*blockDim.z + bx_lo[2];

         int i = tid_x + stride_x;
         int j = tid_y + stride_y;
         int k = tid_z + stride_z;

         IntVect sm_begin = IntVect::TheZeroVector();
         IntVect sm_end(blockDim.x+2, blockDim.y+2, blockDim.z+2);

         get_index idx_eps(sm_begin.dim3(), sm_end.dim3(), 1);
         get_index idx_vfrac(sm_begin.dim3(), sm_end.dim3(), 1);
         get_index idx_flags(sm_begin.dim3(), sm_end.dim3(), 1);

         extern __shared__ Real sm[];

         const Box sm_box(sm_begin, sm_end - IntVect::TheUnitVector());
         const Long BlocksSize = sm_box.numPts();

         Real* eps_sm = sm;
         Real* vfrac_sm = (Real*)&eps_sm[BlocksSize];
         EBCellFlag* flags_sm = (EBCellFlag*)&vfrac_sm[BlocksSize];

         if((i <= bx_hi[0]) and (j <= bx_hi[1]) and (k <= bx_hi[2]))
         {
           tid_x += 1; tid_y += 1; tid_z += 1;

           int loc_x[3] = {0, tid_x, amrex::min((int)blockDim.x, bx_hi[0]-stride_x+1)+1};
           int loc_y[3] = {0, tid_y, amrex::min((int)blockDim.y, bx_hi[1]-stride_y+1)+1};
           int loc_z[3] = {0, tid_z, amrex::min((int)blockDim.z, bx_hi[2]-stride_z+1)+1};

           int glob_x[3] = {stride_x-1, i, amrex::min(stride_x+(int)blockDim.x, bx_hi[0]+1)};
           int glob_y[3] = {stride_y-1, j, amrex::min(stride_y+(int)blockDim.y, bx_hi[1]+1)};
           int glob_z[3] = {stride_z-1, k, amrex::min(stride_z+(int)blockDim.z, bx_hi[2]+1)};

           for (int ii(0); ii < 3; ii++)
           for (int jj(0); jj < 3; jj++)
           for (int kk(0); kk < 3; kk++) {
             int t_x(loc_x[ii]);
             int t_y(loc_y[jj]);
             int t_z(loc_z[kk]);

             int g_x(glob_x[ii]);
             int g_y(glob_y[jj]);
             int g_z(glob_z[kk]);

             eps_sm[idx_eps(t_x,t_y,t_z)] = ep_s(g_x,g_y,g_z);
             vfrac_sm[idx_vfrac(t_x,t_y,t_z)] = vfrac(g_x,g_y,g_z);
             flags_sm[idx_flags(t_x,t_y,t_z)] = flags(g_x,g_y,g_z);
           }

           __syncthreads();
#else
       amrex::ParallelFor(bx,
         [flags,ep_s,mf_redist,vfrac,duplicate,scale_array,max_eps,ncomp,
          mask_box] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         const Real* eps_sm = ep_s.dataPtr();
         const Real* vfrac_sm = vfrac.dataPtr();
         const EBCellFlag* flags_sm = flags.dataPtr();

         get_index idx_eps(ep_s.begin, ep_s.end, 1);
         get_index idx_vfrac(vfrac.begin, vfrac.end, 1);
         get_index idx_flags(flags.begin, flags.end, 1);

         const int tid_x = i;
         const int tid_y = j;
         const int tid_z = k;

         {
#endif
           if(flags_sm[idx_flags(tid_x,tid_y,tid_z)].isSingleValued() and
              eps_sm[idx_eps(tid_x,tid_y,tid_z)] > max_eps)
           {
             amrex::Real sum_vfrac_eps = 0.0;
             amrex::Real sum_vfrac     = 0.0;

             for(int ii(-1); ii <= 1; ii++)
             for(int jj(-1); jj <= 1; jj++)
             for(int kk(-1); kk <= 1; kk++) {
               if((ii != 0 or jj != 0 or kk != 0 ) and
                  flags_sm[idx_flags(tid_x,tid_y,tid_z)].isConnected({ii,jj,kk}) and
                  mask_box.contains(IntVect(i+ii,j+jj,k+kk)) and
                  ((not flags_sm[idx_flags(tid_x+ii,tid_y+jj,tid_z+kk)].isSingleValued()) or
                   (eps_sm[idx_eps(tid_x+ii,tid_y+jj,tid_z+kk)] <= max_eps)))
               {
                 const Real vfrac_val = vfrac_sm[idx_vfrac(tid_x+ii,tid_y+jj,tid_z+kk)];
                 const Real eps_val = eps_sm[idx_eps(tid_x+ii,tid_y+jj,tid_z+kk)];

                 sum_vfrac     += vfrac_val;
                 sum_vfrac_eps += vfrac_val*eps_val;
               }
             }

             // Average volume fraction in the neighborhood around of the packed
             // cell. This value is capped by the user-defined max pack value.
             amrex::Real avg_eps = amrex::min(max_eps, sum_vfrac_eps / sum_vfrac);

             // Fraction of material we want to redistribute
             const Real eps_val = eps_sm[idx_eps(tid_x,tid_y,tid_z)];
             amrex::Real scale =  amrex::max(0.0, eps_val-avg_eps)/eps_val;
             scale_array(i,j,k) = scale;

             for(int n(0); n < ncomp; n++) {

               // This is the amount of the multifab we are redistributing
               // that we are moving into the packed cell's neighborhood.
               Real overflow =
                 duplicate(i,j,k,n)*scale*vfrac_sm[idx_vfrac(tid_x,tid_y,tid_z)]/sum_vfrac;

               for(int kk(-1); kk <= 1; kk++)
               for(int jj(-1); jj <= 1; jj++)
               for(int ii(-1); ii <= 1; ii++) {
                 if((ii != 0 or jj != 0 or kk != 0) and
                    flags_sm[idx_flags(tid_x,tid_y,tid_z)].isConnected({ii,jj,kk}) and
                    mask_box.contains(IntVect(i+ii,j+jj,k+kk)) and
                    ((not flags_sm[idx_flags(tid_x+ii,tid_y+jj,tid_z+kk)].isSingleValued()) or
                     (eps_sm[idx_eps(tid_x+ii,tid_y+jj,tid_z+kk)] <= max_eps)))
                 {
                   Gpu::Atomic::Add(&mf_redist(i+ii,j+jj,k+kk,n),
                                    overflow);
                 }
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
