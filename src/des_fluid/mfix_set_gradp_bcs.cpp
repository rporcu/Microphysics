#include <mfix_set_gradp_bcs.hpp>

using namespace amrex;

void
set_gradp_bcs (const Box& bx,
               Array4<Real> const& gp,
               Array4<int> const& bct_ilo, 
               Array4<int> const& bct_ihi,
               Array4<int> const& bct_jlo,
               Array4<int> const& bct_jhi,
               Array4<int> const& bct_klo,
               Array4<int> const& bct_khi,
               Box& domain,
               BcList& bc_list,
               const int* ng)
{
  // Extract the lower and upper boundaries of Box bx and Domain
  const IntVect bx_lo(bx.loVect()), bx_hi(bx.hiVect());
  const amrex::Dim3 dom_lo = amrex::lbound(domain);
  const amrex::Dim3 dom_hi = amrex::ubound(domain);

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

  if(bx_lo[0] <= dom_lo.x)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_3D(bx_yz, i, j, k,
    {
      if((bct_ilo(j,k,0) == bc_list.pinf) or (bct_ilo(j,k,0) == bc_list.pout))
      {
        gp(dom_lo.x-1,j,k,0) = gp(dom_lo.x,j,k,0);
        gp(dom_lo.x-1,j,k,1) = gp(dom_lo.x,j,k,1);
        gp(dom_lo.x-1,j,k,2) = gp(dom_lo.x,j,k,2);
      }
      else if(bct_ilo(j,k,0) == bc_list.minf)
      {
        gp(dom_lo.x-1,j,k,0) = gp(dom_lo.x,j,k,0);
        gp(dom_lo.x-1,j,k,1) = 0;
        gp(dom_lo.x-1,j,k,2) = 0;
      }
    });
  }

  if(bx_hi[0] >= dom_hi.x+1)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_3D(bx_yz, i, j, k,
    {
      if((bct_ihi(j,k,0) == bc_list.pinf) or (bct_ihi(j,k,0) == bc_list.pout))
      {
        gp(dom_hi.x+1,j,k,0) = gp(dom_hi.x,j,k,0);
        gp(dom_hi.x+1,j,k,1) = gp(dom_hi.x,j,k,1);
        gp(dom_hi.x+1,j,k,2) = gp(dom_hi.x,j,k,2);
      }
      else if(bct_ihi(j,k,0) == bc_list.minf)
      {
        gp(dom_hi.x+1,j,k,0) = gp(dom_hi.x,j,k,0);
        gp(dom_hi.x+1,j,k,1) = 0;
        gp(dom_hi.x+1,j,k,2) = 0;
      }
    });
  }

  if(bx_lo[1] <= dom_lo.y)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_3D(bx_xz, i, j, k,
    {
      if((bct_jlo(i,k,0) == bc_list.pinf) or (bct_jlo(i,k,0) == bc_list.pout))
      {
        gp(i,dom_lo.y-1,k,0) = gp(i,dom_lo.y,k,0);
        gp(i,dom_lo.y-1,k,1) = gp(i,dom_lo.y,k,1);
        gp(i,dom_lo.y-1,k,2) = gp(i,dom_lo.y,k,2);
      }
      else if(bct_jlo(i,k,0) == bc_list.minf)
      {
        gp(i,dom_lo.y-1,k,0) = 0;
        gp(i,dom_lo.y-1,k,1) = gp(i,dom_lo.y,k,1);
        gp(i,dom_lo.y-1,k,2) = 0;
      }
    });
  }

  if(bx_hi[1] >= dom_hi.y+1)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_3D(bx_xz, i, j, k,
    {
      if((bct_jhi(i,k,0) == bc_list.pinf) or (bct_jhi(i,k,0) == bc_list.pout))
      {
        gp(i,dom_hi.y+1,k,0) = gp(i,dom_hi.y,k,0);
        gp(i,dom_hi.y+1,k,1) = gp(i,dom_hi.y,k,1);
        gp(i,dom_hi.y+1,k,2) = gp(i,dom_hi.y,k,2);
      }
      else if(bct_jhi(i,k,0) == bc_list.minf)
      {
        gp(i,dom_hi.y+1,k,0) = 0;
        gp(i,dom_hi.y+1,k,1) = gp(i,dom_hi.y,k,1);
        gp(i,dom_hi.y+1,k,2) = 0;
      }
    });
  }

  if(bx_lo[2] <= dom_lo.z)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_3D(bx_xy, i, j, k,
    {
      if((bct_klo(i,j,0) == bc_list.pinf) or (bct_klo(i,j,0) == bc_list.pout))
      {
        gp(i,j,dom_lo.z-1,0) = gp(i,j,dom_lo.z,0);
        gp(i,j,dom_lo.z-1,1) = gp(i,j,dom_lo.z,1);
        gp(i,j,dom_lo.z-1,2) = gp(i,j,dom_lo.z,2);
      }
      else if(bct_klo(i,j,0) == bc_list.minf)
      {
        gp(i,j,dom_lo.z-1,0) = 0;
        gp(i,j,dom_lo.z-1,1) = 0;
        gp(i,j,dom_lo.z-1,2) = gp(i,j,dom_lo.z,2);
      }
    });
  }

  if(bx_hi[2] >= dom_hi.z+1)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_3D(bx_xy, i, j, k,
    {
      if((bct_khi(i,j,0) == bc_list.pinf) or (bct_khi(i,j,0) == bc_list.pout))
      {
        gp(i,j,dom_hi.z+1,0) = gp(i,j,dom_hi.z,0);
        gp(i,j,dom_hi.z+1,1) = gp(i,j,dom_hi.z,1);
        gp(i,j,dom_hi.z+1,2) = gp(i,j,dom_hi.z,2);
      }
      else if(bct_khi(i,j,0) == bc_list.minf)
      {
        gp(i,j,dom_hi.z+1,0) = 0;
        gp(i,j,dom_hi.z+1,1) = 0;
        gp(i,j,dom_hi.z+1,2) = gp(i,j,dom_hi.z,2);
      }
    });
  }
}
