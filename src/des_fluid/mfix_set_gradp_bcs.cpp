#include <mfix_set_gradp_bcs.hpp>

using namespace amrex;

void
set_gradp_bcs (const Box& bx,
               FArrayBox& gp_fab,
               IArrayBox& bct_ilo_fab, 
               IArrayBox& bct_ihi_fab,
               IArrayBox& bct_jlo_fab,
               IArrayBox& bct_jhi_fab,
               IArrayBox& bct_klo_fab,
               IArrayBox& bct_khi_fab,
               Box& domain,
               BcList& bc_list,
               const int* ng)
{
  // Extract the lower and upper boundaries of Box bx and Domain
  const IntVect bx_lo(bx.loVect()), bx_hi(bx.hiVect());
  const IntVect dom_lo(domain.loVect()), dom_hi(domain.hiVect());

  Array4<Real> const& gp = gp_fab.array();

  Array4<int> const& bct_ilo = bct_ilo_fab.array();
  Array4<int> const& bct_ihi = bct_ihi_fab.array();
  Array4<int> const& bct_jlo = bct_jlo_fab.array();
  Array4<int> const& bct_jhi = bct_jhi_fab.array();
  Array4<int> const& bct_klo = bct_klo_fab.array();
  Array4<int> const& bct_khi = bct_khi_fab.array();

  // Create a 2D Box collapsing bx on x-direction
  IntVect bx_yz_hi(bx_hi);
  bx_yz_hi[0] = bx_lo[0];
  const Box bx_yz(bx_lo, bx_yz_hi);

  // Create a 2D Box collapsing bx on y-direction
  IntVect bx_xz_hi(bx_hi);
  bx_xz_hi[1] = bx_lo[1];
  const Box bx_xz(bx_lo, bx_xz_hi);

  // Create a 2D Box collapsing bx on z-direction
  IntVect bx_xy_hi(bx_hi);
  bx_xy_hi[2] = bx_lo[2];
  const Box bx_xy(bx_lo, bx_xy_hi);

  if(bx_lo[0] <= dom_lo[0])
  {
    AMREX_CUDA_HOST_DEVICE_FOR_3D(bx_yz, i, j, k,
    {
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        gp(dom_lo[0]-1,j,k,0) = gp(dom_lo[0],j,k,0);
        gp(dom_lo[0]-1,j,k,1) = gp(dom_lo[0],j,k,1);
        gp(dom_lo[0]-1,j,k,2) = gp(dom_lo[0],j,k,2);
      }
      else if(bct == bc_list.minf)
      {
        gp(dom_lo[0]-1,j,k,0) = gp(dom_lo[0],j,k,0);
        gp(dom_lo[0]-1,j,k,1) = 0;
        gp(dom_lo[0]-1,j,k,2) = 0;
      }
    });
  }

  if(bx_hi[0] >= dom_hi[0]+1)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_3D(bx_yz, i, j, k,
    {
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        gp(dom_hi[0]+1,j,k,0) = gp(dom_hi[0],j,k,0);
        gp(dom_hi[0]+1,j,k,1) = gp(dom_hi[0],j,k,1);
        gp(dom_hi[0]+1,j,k,2) = gp(dom_hi[0],j,k,2);
      }
      else if(bct == bc_list.minf)
      {
        gp(dom_hi[0]+1,j,k,0) = gp(dom_hi[0],j,k,0);
        gp(dom_hi[0]+1,j,k,1) = 0;
        gp(dom_hi[0]+1,j,k,2) = 0;
      }
    });
  }

  if(bx_lo[1] <= dom_lo[1])
  {
    AMREX_CUDA_HOST_DEVICE_FOR_3D(bx_xz, i, j, k,
    {
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        gp(i,dom_lo[1]-1,k,0) = gp(i,dom_lo[1],k,0);
        gp(i,dom_lo[1]-1,k,1) = gp(i,dom_lo[1],k,1);
        gp(i,dom_lo[1]-1,k,2) = gp(i,dom_lo[1],k,2);
      }
      else if(bct == bc_list.minf)
      {
        gp(i,dom_lo[1]-1,k,0) = 0;
        gp(i,dom_lo[1]-1,k,1) = gp(i,dom_lo[1],k,1);
        gp(i,dom_lo[1]-1,k,2) = 0;
      }
    });
  }

  if(bx_hi[1] >= dom_hi[1]+1)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_3D(bx_xz, i, j, k,
    {
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        gp(i,dom_hi[1]+1,k,0) = gp(i,dom_hi[1],k,0);
        gp(i,dom_hi[1]+1,k,1) = gp(i,dom_hi[1],k,1);
        gp(i,dom_hi[1]+1,k,2) = gp(i,dom_hi[1],k,2);
      }
      else if(bct == bc_list.minf)
      {
        gp(i,dom_hi[1]+1,k,0) = 0;
        gp(i,dom_hi[1]+1,k,1) = gp(i,dom_hi[1],k,1);
        gp(i,dom_hi[1]+1,k,2) = 0;
      }
    });
  }

  if(bx_lo[2] <= dom_lo[2])
  {
    AMREX_CUDA_HOST_DEVICE_FOR_3D(bx_xy, i, j, k,
    {
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        gp(i,j,dom_lo[2]-1,0) = gp(i,j,dom_lo[2],0);
        gp(i,j,dom_lo[2]-1,1) = gp(i,j,dom_lo[2],1);
        gp(i,j,dom_lo[2]-1,2) = gp(i,j,dom_lo[2],2);
      }
      else if(bct == bc_list.minf)
      {
        gp(i,j,dom_lo[2]-1,0) = 0;
        gp(i,j,dom_lo[2]-1,1) = 0;
        gp(i,j,dom_lo[2]-1,2) = gp(i,j,dom_lo[2],2);
      }
    });
  }

  if(bx_hi[2] >= dom_hi[2]+1)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_3D(bx_xy, i, j, k,
    {
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        gp(i,j,dom_hi[2]+1,0) = gp(i,j,dom_hi[2],0);
        gp(i,j,dom_hi[2]+1,1) = gp(i,j,dom_hi[2],1);
        gp(i,j,dom_hi[2]+1,2) = gp(i,j,dom_hi[2],2);
      }
      else if(bct == bc_list.minf)
      {
        gp(i,j,dom_hi[2]+1,0) = 0;
        gp(i,j,dom_hi[2]+1,1) = 0;
        gp(i,j,dom_hi[2]+1,2) = gp(i,j,dom_hi[2],2);
      }
    });
  }
}
