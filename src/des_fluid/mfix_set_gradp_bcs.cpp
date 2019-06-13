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
  // Extract the lower and upper boundaries of Domain
  const IntVect dom_lo(domain.loVect()), dom_hi(domain.hiVect());

  Array4<Real> const& gp = gp_fab.array();
  const IntVect gp_lo(gp_fab.loVect()), gp_hi(gp_fab.hiVect());

  Array4<int> const& bct_ilo = bct_ilo_fab.array();
  Array4<int> const& bct_ihi = bct_ihi_fab.array();
  Array4<int> const& bct_jlo = bct_jlo_fab.array();
  Array4<int> const& bct_jhi = bct_jhi_fab.array();
  Array4<int> const& bct_klo = bct_klo_fab.array();
  Array4<int> const& bct_khi = bct_khi_fab.array();

  // Create 2D low and hi Boxes
  IntVect bx_lo_yz_lo(gp_lo), bx_lo_yz_hi(gp_hi);
  IntVect bx_hi_yz_lo(gp_lo), bx_hi_yz_hi(gp_hi);
  bx_lo_yz_lo[0] = dom_lo[0]-1;
  bx_lo_yz_hi[0] = dom_lo[0]-1;
  bx_hi_yz_lo[0] = dom_hi[0]+1;
  bx_hi_yz_hi[0] = dom_hi[0]+1;
  const Box bx_yz_lo(bx_lo_yz_lo, bx_lo_yz_hi);
  const Box bx_yz_hi(bx_hi_yz_lo, bx_hi_yz_hi);

  // Create 2D low and hi Boxes
  IntVect bx_lo_xz_lo(gp_lo), bx_lo_xz_hi(gp_hi);
  IntVect bx_hi_xz_lo(gp_lo), bx_hi_xz_hi(gp_hi);
  bx_lo_xz_lo[1] = dom_lo[1]-1;
  bx_lo_xz_hi[1] = dom_lo[1]-1;
  bx_hi_xz_lo[1] = dom_hi[1]+1;
  bx_hi_xz_hi[1] = dom_hi[1]+1;
  const Box bx_xz_lo(bx_lo_xz_lo, bx_lo_xz_hi);
  const Box bx_xz_hi(bx_hi_xz_lo, bx_hi_xz_hi);

  // Create 2D low and hi Boxes
  IntVect bx_lo_xy_lo(gp_lo), bx_lo_xy_hi(gp_hi);
  IntVect bx_hi_xy_lo(gp_lo), bx_hi_xy_hi(gp_hi);
  bx_lo_xy_lo[2] = dom_lo[2]-1;
  bx_lo_xy_hi[2] = dom_lo[2]-1;
  bx_hi_xy_lo[2] = dom_hi[2]+1;
  bx_hi_xy_hi[2] = dom_hi[2]+1;
  const Box bx_xy_lo(bx_lo_xy_lo, bx_lo_xy_hi);
  const Box bx_xy_hi(bx_hi_xy_lo, bx_hi_xy_hi);

  //TODO possible data races!!!
  // one solution could be to make a copy of gp

  if(gp_lo[0] <= dom_lo[0])
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_yz_lo, i, j, k,
    {
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        gp(i,j,k,0) = gp(i+1,j,k,0);
        gp(i,j,k,1) = gp(i+1,j,k,1);
        gp(i,j,k,2) = gp(i+1,j,k,2);
      }
      else if(bct == bc_list.minf)
      {
        gp(i,j,k,0) = gp(i+1,j,k,0);
        gp(i,j,k,1) = 0;
        gp(i,j,k,2) = 0;
      }
    });
  }

  if(gp_hi[0] >= dom_hi[0]+1)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_yz_hi, i, j, k,
    {
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        gp(i,j,k,0) = gp(i-1,j,k,0);
        gp(i,j,k,1) = gp(i-1,j,k,1);
        gp(i,j,k,2) = gp(i-1,j,k,2);
      }
      else if(bct == bc_list.minf)
      {
        gp(i,j,k,0) = gp(i-1,j,k,0);
        gp(i,j,k,1) = 0;
        gp(i,j,k,2) = 0;
      }
    });
  }

  if(gp_lo[1] <= dom_lo[1])
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xz_lo, i, j, k,
    {
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        gp(i,j,k,0) = gp(i,j+1,k,0);
        gp(i,j,k,1) = gp(i,j+1,k,1);
        gp(i,j,k,2) = gp(i,j+1,k,2);
      }
      else if(bct == bc_list.minf)
      {
        gp(i,j,k,0) = 0;
        gp(i,j,k,1) = gp(i,j+1,k,1);
        gp(i,j,k,2) = 0;
      }
    });
  }

  if(gp_hi[1] >= dom_hi[1]+1)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xz_hi, i, j, k,
    {
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        gp(i,j,k,0) = gp(i,j-1,k,0);
        gp(i,j,k,1) = gp(i,j-1,k,1);
        gp(i,j,k,2) = gp(i,j-1,k,2);
      }
      else if(bct == bc_list.minf)
      {
        gp(i,j,k,0) = 0;
        gp(i,j,k,1) = gp(i,j-1,k,1);
        gp(i,j,k,2) = 0;
      }
    });
  }

  if(gp_lo[2] <= dom_lo[2])
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xy_lo, i, j, k,
    {
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        gp(i,j,k,0) = gp(i,j,k+1,0);
        gp(i,j,k,1) = gp(i,j,k+1,1);
        gp(i,j,k,2) = gp(i,j,k+1,2);
      }
      else if(bct == bc_list.minf)
      {
        gp(i,j,k,0) = 0;
        gp(i,j,k,1) = 0;
        gp(i,j,k,2) = gp(i,j,k+1,2);
      }
    });
  }

  if(gp_hi[2] >= dom_hi[2]+1)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xy_hi, i, j, k,
    {
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        gp(i,j,k,0) = gp(i,j,k-1,0);
        gp(i,j,k,1) = gp(i,j,k-1,1);
        gp(i,j,k,2) = gp(i,j,k-1,2);
      }
      else if(bct == bc_list.minf)
      {
        gp(i,j,k,0) = 0;
        gp(i,j,k,1) = 0;
        gp(i,j,k,2) = gp(i,j,k-1,2);
      }
    });
  }
}
