#include <mfix_set_velocity_bcs.hpp>
#include <bc_mod_F.H>

//              
//  This subroutine sets the BCs for the velocity components only.
//

//void
//set_velocity_bcs(Real* time,
//                 const BcList& bc_list,
//                 FArrayBox& vel_fab,
//                 const IArrayBox& bct_ilo_fab,
//                 const IArrayBox& bct_ihi_fab,
//                 const IArrayBox& bct_jlo_fab,
//                 const IArrayBox& bct_jhi_fab,
//                 const IArrayBox& bct_klo_fab,
//                 const IArrayBox& bct_khi_fab,
//                 const Box& domain,
//                 const int* ng,
//                 const int* extrap_dir_bcs)
//{
//  IntVect dom_lo(domain.loVect());
//  IntVect dom_hi(domain.hiVect());
//
//  Array4<Real> const& vel = vel_fab.array();
//
//  IntVect vel_lo(vel_fab.loVect());
//  IntVect vel_hi(vel_fab.hiVect());
//
//  Array4<const int> const& bct_ilo = bct_ilo_fab.array();
//  Array4<const int> const& bct_ihi = bct_ihi_fab.array();
//  Array4<const int> const& bct_jlo = bct_jlo_fab.array();
//  Array4<const int> const& bct_jhi = bct_jhi_fab.array();
//  Array4<const int> const& bct_klo = bct_klo_fab.array();
//  Array4<const int> const& bct_khi = bct_khi_fab.array();
//
//  const int nlft = std::max(0, dom_lo[0]-vel_lo[0]);
//  const int nbot = std::max(0, dom_lo[1]-vel_lo[1]);
//  const int ndwn = std::max(0, dom_lo[2]-vel_lo[2]);
//
//  const int nrgt = std::max(0, vel_hi[0]-dom_hi[0]);
//  const int ntop = std::max(0, vel_hi[1]-dom_hi[1]);
//  const int nup  = std::max(0, vel_hi[2]-dom_hi[2]);
//
//  // Create InVects for following 2D Boxes
//  IntVect bx_yz_lo_lo_2D(vel_lo), bx_yz_lo_hi_2D(vel_hi);
//  IntVect bx_yz_hi_lo_2D(vel_lo), bx_yz_hi_hi_2D(vel_hi);
//  IntVect bx_xz_lo_lo_2D(vel_lo), bx_xz_lo_hi_2D(vel_hi);
//  IntVect bx_xz_hi_lo_2D(vel_lo), bx_xz_hi_hi_2D(vel_hi);
//  IntVect bx_xy_lo_lo_2D(vel_lo), bx_xy_lo_hi_2D(vel_hi);
//  IntVect bx_xy_hi_lo_2D(vel_lo), bx_xy_hi_hi_2D(vel_hi);
//
//  // Fix lo and hi limits
//  bx_yz_lo_lo_2D[0] = dom_lo[0]-1;
//  bx_yz_lo_hi_2D[0] = dom_lo[0]-1;
//  bx_yz_hi_lo_2D[0] = dom_hi[0]+1;
//  bx_yz_hi_hi_2D[0] = dom_hi[0]+1;
//
//  bx_xz_lo_lo_2D[1] = dom_lo[1]-1;
//  bx_xz_lo_hi_2D[1] = dom_lo[1]-1;
//  bx_xz_hi_lo_2D[1] = dom_hi[1]+1;
//  bx_xz_hi_hi_2D[1] = dom_hi[1]+1;
//
//  bx_xy_lo_lo_2D[2] = dom_lo[2]-1;
//  bx_xy_lo_hi_2D[2] = dom_lo[2]-1;
//  bx_xy_hi_lo_2D[2] = dom_hi[2]+1;
//  bx_xy_hi_hi_2D[2] = dom_hi[2]+1;
//
//  // Create 2D boxes for CUDA loops
//  const Box bx_yz_lo_2D(bx_yz_lo_lo_2D, bx_yz_lo_hi_2D);
//  const Box bx_yz_hi_2D(bx_yz_hi_lo_2D, bx_yz_hi_hi_2D);
//
//  const Box bx_xz_lo_2D(bx_xz_lo_lo_2D, bx_xz_lo_hi_2D);
//  const Box bx_xz_hi_2D(bx_xz_hi_lo_2D, bx_xz_hi_hi_2D);
//
//  const Box bx_xy_lo_2D(bx_xy_lo_lo_2D, bx_xy_lo_hi_2D);
//  const Box bx_xy_hi_2D(bx_xy_hi_lo_2D, bx_xy_hi_hi_2D);
//
//  // Create InVects for following 3D Boxes
//  IntVect bx_yz_lo_hi_3D(vel_hi), bx_xz_lo_hi_3D(vel_hi), bx_xy_lo_hi_3D(vel_hi);
//  IntVect bx_yz_hi_lo_3D(vel_lo), bx_xz_hi_lo_3D(vel_lo), bx_xy_hi_lo_3D(vel_lo);
//
//  // Fix lo and hi limits
//  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
//  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;
//
//  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
//  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;
//
//  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
//  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;
//
//  // Create 3D boxes for CUDA loops
//  const Box bx_yz_lo_3D(vel_lo, bx_yz_lo_hi_3D);
//  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, vel_hi);
//
//  const Box bx_xz_lo_3D(vel_lo, bx_xz_lo_hi_3D);
//  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, vel_hi);
//
//  const Box bx_xy_lo_3D(vel_lo, bx_xy_lo_hi_3D);
//  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, vel_hi);
//
//  mfix_usr1(time);
//
//  if (nlft > 0)
//  {
//    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_yz_lo_3D, 3, i, j, k, n,
//    {
//      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
//      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
//
//      if((bct == bc_list.pinf) or (bct == bc_list.pout))
//        vel(i,j,k,n) = vel(dom_lo[0],j,k,n);
//      else if(bct == bc_list.minf)
//      {
//        if(n == 0)
//          vel(i,j,k,n) = get_bc_u_g(bcv);
//        else
//          vel(i,j,k,n) = 0;
//      }
//    });
//
//    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_yz_lo_2D, 3, i, j, k, n,
//    {
//      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
//
//      if(bct == bc_list.minf)
//        vel(i,j,k,n) = 2*vel(i,j,k,n) - vel(i+1,j,k,n);
//    });
//  }
//
//  if (nrgt > 0)
//  {
//    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_yz_hi_3D, 3, i, j, k, n,
//    {
//      const int bcv = bct_ilo(dom_hi[0]+1,j,k,1);
//      const int bct = bct_ilo(dom_hi[0]+1,j,k,0);
//
//      if((bct == bc_list.pinf) or (bct == bc_list.pout))
//        vel(i,j,k,n) = vel(dom_hi[0],j,k,n);
//      else if(bct == bc_list.minf)
//      {
//        if(n == 0)
//          vel(i,j,k,n) = get_bc_u_g(bcv);
//        else
//          vel(i,j,k,n) = 0;
//      }
//    });
//
//    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_yz_hi_2D, 3, i, j, k, n,
//    {
//      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
//
//      if(bct == bc_list.minf)
//        vel(i,j,k,n) = 2*vel(i,j,k,n) - vel(i-1,j,k,n);
//    });
//  }
//
//  if (nbot > 0)
//  {
//    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_xz_lo_3D, 3, i, j, k, n,
//    {
//      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
//      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
//
//      if((bct == bc_list.pinf) or (bct == bc_list.pout))
//        vel(i,j,k,n) = vel(i,dom_lo[1],k,n);
//      else if(bct == bc_list.minf)
//      {
//        if(n == 1)
//          vel(i,j,k,n) = get_bc_v_g(bcv);
//        else
//          vel(i,j,k,n) = 0;
//      }
//    });
//
//    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_xz_lo_2D, 3, i, j, k, n,
//    {
//      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
//
//      if(bct == bc_list.minf)
//        vel(i,j,k,n) = 2*vel(i,j,k,n) - vel(i,j+1,k,n);
//    });
//  }
//
//  if (ntop > 0)
//  {
//    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_xz_hi_3D, 3, i, j, k, n,
//    {
//      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
//      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
//
//      if((bct == bc_list.pinf) or (bct == bc_list.pout))
//        vel(i,j,k,n) = vel(i,dom_hi[1],k,n);
//      else if(bct == bc_list.minf)
//      {
//        if(n == 1)
//          vel(i,j,k,n) = get_bc_v_g(bcv);
//        else
//          vel(i,j,k,n) = 0;
//      }
//    });
//
//    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_xz_hi_2D, 3, i, j, k, n,
//    {
//      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
//
//      if(bct == bc_list.minf)
//        vel(i,j,k,n) = 2*vel(i,j,k,n) - vel(i,j-1,k,n);
//    });
//  }
//
//  if (ndwn > 0)
//  {
//    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_xy_lo_3D, 3, i, j, k, n,
//    {
//      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
//      const int bct = bct_klo(i,j,dom_lo[2]-1,0);
//
//      if((bct == bc_list.pinf) or (bct == bc_list.pout))
//        vel(i,j,k,n) = vel(i,j,dom_lo[2],n);
//      else if(bct == bc_list.minf)
//      {
//        if(n == 2)
//          vel(i,j,k,n) = get_bc_w_g(bcv);
//        else
//          vel(i,j,k,n) = 0;
//      }
//    });
//
//    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_xy_lo_2D, 3, i, j, k, n,
//    {
//      const int bct = bct_klo(i,j,dom_lo[2]-1,0);
//
//      if(bct == bc_list.minf)
//        vel(i,j,k,n) = 2*vel(i,j,k,n) - vel(i,j,k+1,n);
//    });
//  }
//
//  if (nup > 0)
//  {
//    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_xy_hi_3D, 3, i, j, k, n,
//    {
//      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
//      const int bct = bct_khi(i,j,dom_hi[2]+1,0);
//
//      if((bct == bc_list.pinf) or (bct == bc_list.pout))
//        vel(i,j,k,n) = vel(i,j,dom_hi[2],n);
//      else if(bct == bc_list.minf)
//      {
//        if(n == 2)
//          vel(i,j,k,n) = get_bc_w_g(bcv);
//        else
//          vel(i,j,k,n) = 0;
//      }
//    });
//
//    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_xy_hi_2D, 3, i, j, k, n,
//    {
//      const int bct = bct_khi(i,j,dom_hi[2]+1,0);
//
//      if(bct == bc_list.minf)
//        vel(i,j,k,n) = 2*vel(i,j,k,n) - vel(i,j,k-1,n);
//    });
//  }
//}


void
set_vec_bcs(const BcList& bc_list,
            FArrayBox& vec_fab,
            const IArrayBox& bct_ilo_fab,
            const IArrayBox& bct_ihi_fab,
            const IArrayBox& bct_jlo_fab,
            const IArrayBox& bct_jhi_fab,
            const IArrayBox& bct_klo_fab,
            const IArrayBox& bct_khi_fab,
            const Box& domain,
            const int* ng)
{
  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<Real> const& vec = vec_fab.array();

  IntVect vec_lo(vec_fab.loVect());
  IntVect vec_hi(vec_fab.hiVect());

  Array4<const int> const& bct_ilo = bct_ilo_fab.array();
  Array4<const int> const& bct_ihi = bct_ihi_fab.array();
  Array4<const int> const& bct_jlo = bct_jlo_fab.array();
  Array4<const int> const& bct_jhi = bct_jhi_fab.array();
  Array4<const int> const& bct_klo = bct_klo_fab.array();
  Array4<const int> const& bct_khi = bct_khi_fab.array();

  const int nlft = std::max(0, dom_lo[0]-vec_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-vec_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-vec_lo[2]);

  const int nrgt = std::max(0, vec_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, vec_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, vec_hi[2]-dom_hi[2]);

  // Create InVects for following 3D Boxes
  IntVect bx_yz_lo_hi_3D(vec_hi), bx_xz_lo_hi_3D(vec_hi), bx_xy_lo_hi_3D(vec_hi);
  IntVect bx_yz_hi_lo_3D(vec_lo), bx_xz_hi_lo_3D(vec_lo), bx_xy_hi_lo_3D(vec_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for CUDA loops
  const Box bx_yz_lo_3D(vec_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, vec_hi);

  const Box bx_xz_lo_3D(vec_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, vec_hi);

  const Box bx_xy_lo_3D(vec_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, vec_hi);

  if (nlft > 0)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_yz_lo_3D, 3, i, j, k, n,
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        vec(i,j,k,n) = vec(dom_lo[0],j,k,n);
      else if(bct == bc_list.minf)
        vec(i,j,k,n) = get_bc_ep_g(bcv) * get_bc_vel_g(bcv,n);
    });
  }

  if (nrgt > 0)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_yz_hi_3D, 3, i, j, k, n,
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        vec(i,j,k,n) = vec(dom_hi[0],j,k,n);
      else if(bct == bc_list.minf)
        vec(i,j,k,n) = get_bc_ep_g(bcv) * get_bc_vel_g(bcv,n);
    });
  }

  if (nbot > 0)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_xz_lo_3D, 3, i, j, k, n,
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        vec(i,j,k,n) = vec(i,dom_lo[1],k,n);
      else if(bct == bc_list.minf)
        vec(i,j,k,n) = get_bc_ep_g(bcv) * get_bc_vel_g(bcv,n);
    });
  }

  if (ntop > 0)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_xz_hi_3D, 3, i, j, k, n,
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        vec(i,j,k,n) = vec(i,dom_hi[1],k,n);
      else if(bct == bc_list.minf)
        vec(i,j,k,n) = get_bc_ep_g(bcv) * get_bc_vel_g(bcv,n);
    });
  }

  if (ndwn > 0)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_xy_lo_3D, 3, i, j, k, n,
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        vec(i,j,k,n) = vec(i,j,dom_lo[2],n);
      else if(bct == bc_list.minf)
        vec(i,j,k,n) = get_bc_ep_g(bcv) * get_bc_vel_g(bcv,n);
    });
  }

  if (nup > 0)
  {
    AMREX_CUDA_HOST_DEVICE_FOR_4D(bx_xy_hi_3D, 3, i, j, k, n,
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        vec(i,j,k,n) = vec(i,j,dom_hi[2],n);
      else if(bct == bc_list.minf)
        vec(i,j,k,n) = get_bc_ep_g(bcv) * get_bc_vel_g(bcv,n);
    });
  }
}
