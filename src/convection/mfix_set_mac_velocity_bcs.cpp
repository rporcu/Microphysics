// 
//  This subroutine sets the BCs for velocity components only.
//  
//  Author: Michele Rosso
// 
//  Date: December 20, 2017
// 

#include <mfix_set_mac_velocity_bcs.hpp>
#include <bc_mod_F.H>

void
set_mac_velocity_bcs(Real* time,
                     BcList& bc_list,
                     const Box& bx,
                     MFIter* mfi,
                     MultiFab& u_g_fab,
                     MultiFab& v_g_fab,
                     MultiFab& w_g_fab,
                     IArrayBox& bct_ilo_fab,
                     IArrayBox& bct_ihi_fab,
                     IArrayBox& bct_jlo_fab,
                     IArrayBox& bct_jhi_fab,
                     IArrayBox& bct_klo_fab,
                     IArrayBox& bct_khi_fab,
                     Box& domain,
                     Real** m_bc_vel_g,
                     const int* nghost)
{
  IntVect bx_lo(bx.loVect());
  IntVect bx_hi(bx.hiVect());

  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<Real> const& u_g = u_g_fab.array(*mfi);
  Array4<Real> const& v_g = v_g_fab.array(*mfi);
  Array4<Real> const& w_g = w_g_fab.array(*mfi);

  IntVect u_lo(u_g_fab[*mfi].loVect());
  IntVect u_hi(u_g_fab[*mfi].hiVect());

  IntVect v_lo(v_g_fab[*mfi].loVect());
  IntVect v_hi(v_g_fab[*mfi].hiVect());

  IntVect w_lo(w_g_fab[*mfi].loVect());
  IntVect w_hi(w_g_fab[*mfi].hiVect());

  Array4<int> const& bct_ilo = bct_ilo_fab.array();
  Array4<int> const& bct_ihi = bct_ihi_fab.array();
  Array4<int> const& bct_jlo = bct_jlo_fab.array();
  Array4<int> const& bct_jhi = bct_jhi_fab.array();
  Array4<int> const& bct_klo = bct_klo_fab.array();
  Array4<int> const& bct_khi = bct_khi_fab.array();

  const int nlft = std::max(0, dom_lo[0]-bx_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-bx_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-bx_lo[2]);

  const int nrgt = std::max(0, bx_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, bx_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, bx_hi[2]-dom_hi[2]);

  // Create InVects for following Boxes
  IntVect ulo_bx_yz_lo(bx_lo), vlo_bx_yz_lo(bx_lo), wlo_bx_yz_lo(bx_lo);
  IntVect ulo_bx_yz_hi(bx_hi), vlo_bx_yz_hi(bx_hi), wlo_bx_yz_hi(bx_hi);
  IntVect uhi_bx_yz_lo(bx_lo), vhi_bx_yz_lo(bx_lo), whi_bx_yz_lo(bx_lo);
  IntVect uhi_bx_yz_hi(bx_hi), vhi_bx_yz_hi(bx_hi), whi_bx_yz_hi(bx_hi);

  IntVect ulo_bx_xz_lo(bx_lo), vlo_bx_xz_lo(bx_lo), wlo_bx_xz_lo(bx_lo);
  IntVect ulo_bx_xz_hi(bx_hi), vlo_bx_xz_hi(bx_hi), wlo_bx_xz_hi(bx_hi);
  IntVect uhi_bx_xz_lo(bx_lo), vhi_bx_xz_lo(bx_lo), whi_bx_xz_lo(bx_lo);
  IntVect uhi_bx_xz_hi(bx_hi), vhi_bx_xz_hi(bx_hi), whi_bx_xz_hi(bx_hi);

  IntVect ulo_bx_xy_lo(bx_lo), vlo_bx_xy_lo(bx_lo), wlo_bx_xy_lo(bx_lo);
  IntVect ulo_bx_xy_hi(bx_hi), vlo_bx_xy_hi(bx_hi), wlo_bx_xy_hi(bx_hi);
  IntVect uhi_bx_xy_lo(bx_lo), vhi_bx_xy_lo(bx_lo), whi_bx_xy_lo(bx_lo);
  IntVect uhi_bx_xy_hi(bx_hi), vhi_bx_xy_hi(bx_hi), whi_bx_xy_hi(bx_hi);

  // Fix lo and hi limits
  // Box 'yz'
  ulo_bx_yz_lo[0] = u_lo[0];
  vlo_bx_yz_lo[0] = v_lo[0];
  wlo_bx_yz_lo[0] = w_lo[0];
  ulo_bx_yz_hi[0] = dom_lo[0];
  vlo_bx_yz_hi[0] = dom_lo[0]-1;
  wlo_bx_yz_hi[0] = dom_lo[0]-1;

  uhi_bx_yz_lo[0] = dom_hi[0]+1;
  vhi_bx_yz_lo[0] = dom_hi[0]+1;
  whi_bx_yz_lo[0] = dom_hi[0]+1;
  uhi_bx_yz_hi[0] = u_hi[0];
  vhi_bx_yz_hi[0] = v_hi[0];
  whi_bx_yz_hi[0] = w_hi[0];

  // Box 'xz'
  ulo_bx_xz_lo[1] = u_lo[1];
  vlo_bx_xz_lo[1] = v_lo[1];
  wlo_bx_xz_lo[1] = w_lo[1];
  ulo_bx_xz_hi[1] = dom_lo[1]-1;
  vlo_bx_xz_hi[1] = dom_lo[1];
  wlo_bx_xz_hi[1] = dom_lo[1]-1;

  uhi_bx_xz_lo[1] = dom_hi[1]+1;
  vhi_bx_xz_lo[1] = dom_hi[1]+1;
  whi_bx_xz_lo[1] = dom_hi[1]+1;
  uhi_bx_xz_hi[1] = u_hi[1];
  vhi_bx_xz_hi[1] = v_hi[1];
  whi_bx_xz_hi[1] = w_hi[1];

  // Box 'xy'
  ulo_bx_xy_lo[2] = u_lo[2];
  vlo_bx_xy_lo[2] = v_lo[2];
  wlo_bx_xy_lo[2] = w_lo[2];
  ulo_bx_xy_hi[2] = dom_lo[2]-1;
  vlo_bx_xy_hi[2] = dom_lo[2]-1;
  wlo_bx_xy_hi[2] = dom_lo[2];

  uhi_bx_xy_lo[2] = dom_hi[2]+1;
  vhi_bx_xy_lo[2] = dom_hi[2]+1;
  whi_bx_xy_lo[2] = dom_hi[2]+1;
  uhi_bx_xy_hi[2] = u_hi[2];
  vhi_bx_xy_hi[2] = v_hi[2];
  whi_bx_xy_hi[2] = w_hi[2];


  // Create 2D boxes for CUDA loops
  const Box ulo_bx_yz(ulo_bx_yz_lo, ulo_bx_yz_hi);
  const Box vlo_bx_yz(vlo_bx_yz_lo, vlo_bx_yz_hi);
  const Box wlo_bx_yz(wlo_bx_yz_lo, wlo_bx_yz_hi);

  const Box uhi_bx_yz(uhi_bx_yz_lo, uhi_bx_yz_hi);
  const Box vhi_bx_yz(vhi_bx_yz_lo, vhi_bx_yz_hi);
  const Box whi_bx_yz(whi_bx_yz_lo, whi_bx_yz_hi);

  const Box ulo_bx_xz(ulo_bx_xz_lo, ulo_bx_xz_hi);
  const Box vlo_bx_xz(vlo_bx_xz_lo, vlo_bx_xz_hi);
  const Box wlo_bx_xz(wlo_bx_xz_lo, wlo_bx_xz_hi);

  const Box uhi_bx_xz(uhi_bx_xz_lo, uhi_bx_xz_hi);
  const Box vhi_bx_xz(vhi_bx_xz_lo, vhi_bx_xz_hi);
  const Box whi_bx_xz(whi_bx_xz_lo, whi_bx_xz_hi);

  const Box ulo_bx_xy(ulo_bx_xy_lo, ulo_bx_xy_hi);
  const Box vlo_bx_xy(vlo_bx_xy_lo, vlo_bx_xy_hi);
  const Box wlo_bx_xy(wlo_bx_xy_lo, wlo_bx_xy_hi);

  const Box uhi_bx_xy(uhi_bx_xy_lo, uhi_bx_xy_hi);
  const Box vhi_bx_xy(vhi_bx_xy_lo, vhi_bx_xy_hi);
  const Box whi_bx_xy(whi_bx_xy_lo, whi_bx_xy_hi);

  mfix_usr1(time);

  if (nlft > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(ulo_bx_yz, i, j, k,
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if(((bct == bc_list.pinf) or (bct == bc_list.pout)) and (i != dom_lo[0]))
        u_g(i,j,k) = u_g(dom_lo[0],j,k);
      else if(bct == bc_list.minf)
        u_g(i,j,k) = m_bc_vel_g[bcv][0];
    });

    AMREX_HOST_DEVICE_FOR_3D(vlo_bx_yz, i, j, k,
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        v_g(i,j,k) = v_g(dom_lo[0],j,k);
      else if(bct == bc_list.minf)
        v_g(i,j,k) = 0;
    });

    AMREX_HOST_DEVICE_FOR_3D(wlo_bx_yz, i, j, k,
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        w_g(i,j,k) = w_g(dom_lo[0],j,k);
      else if(bct == bc_list.minf)
        w_g(i,j,k) = 0;
    });
  }

  if (nrgt > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(uhi_bx_yz, i, j, k,
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if(((bct == bc_list.pinf) or (bct == bc_list.pout)) and (i != dom_hi[0]+1))
        u_g(i,j,k) = u_g(dom_hi[0]+1,j,k);
      else if(bct == bc_list.minf)
        u_g(i,j,k) = m_bc_vel_g[bcv][0];
    });

    AMREX_HOST_DEVICE_FOR_3D(vhi_bx_yz, i, j, k,
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        v_g(i,j,k) = v_g(dom_hi[0],j,k);
      else if(bct == bc_list.minf)
        v_g(i,j,k) = 0;
    });

    AMREX_HOST_DEVICE_FOR_3D(whi_bx_yz, i, j, k,
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        w_g(i,j,k) = w_g(dom_hi[0],j,k);
      else if(bct == bc_list.minf)
        w_g(i,j,k) = 0;
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (nbot > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(ulo_bx_xz, i, j, k,
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        u_g(i,j,k) = u_g(i,dom_lo[1],k);
      else if(bct == bc_list.minf)
        u_g(i,j,k) = 0;
    });

    AMREX_HOST_DEVICE_FOR_3D(vlo_bx_xz, i, j, k,
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if(((bct == bc_list.pinf) or (bct == bc_list.pout)) and (j != dom_lo[1]))
        v_g(i,j,k) = v_g(i,dom_lo[1],k);
      else if(bct == bc_list.minf)
        v_g(i,j,k) = m_bc_vel_g[bcv][1];
    });

    AMREX_HOST_DEVICE_FOR_3D(wlo_bx_xz, i, j, k,
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        w_g(i,j,k) = w_g(i,dom_lo[1],k);
      else if(bct == bc_list.minf)
        w_g(i,j,k) = 0;
    });
  }

  if (ntop > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(uhi_bx_xz, i, j, k,
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        u_g(i,j,k) = u_g(i,dom_hi[1],k);
      else if(bct == bc_list.minf)
        u_g(i,j,k) = 0;
    });

    AMREX_HOST_DEVICE_FOR_3D(vhi_bx_xz, i, j, k,
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if(((bct == bc_list.pinf) or (bct == bc_list.pout)) and (j != dom_hi[1]+1))
        v_g(i,j,k) = v_g(i,dom_hi[1]+1,k);
      else if(bct == bc_list.minf)
        v_g(i,j,k) = m_bc_vel_g[bcv][1];
    });

    AMREX_HOST_DEVICE_FOR_3D(whi_bx_xz, i, j, k,
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        w_g(i,j,k) = w_g(i,dom_hi[1],k);
      else if(bct == bc_list.minf)
        w_g(i,j,k) = 0;
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (ndwn > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(ulo_bx_xy, i, j, k,
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        u_g(i,j,k) = u_g(i,j,dom_lo[2]);
      else if(bct == bc_list.minf)
        u_g(i,j,k) = 0;
    });

    AMREX_HOST_DEVICE_FOR_3D(vlo_bx_xy, i, j, k,
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        v_g(i,j,k) = v_g(i,j,dom_lo[2]);
      else if(bct == bc_list.minf)
        v_g(i,j,k) = 0;
    });

    AMREX_HOST_DEVICE_FOR_3D(wlo_bx_xy, i, j, k,
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if(((bct == bc_list.pinf) or (bct == bc_list.pout)) and (k != dom_lo[2]))
        w_g(i,j,k) = w_g(i,j,dom_lo[2]);
      else if(bct == bc_list.minf)
        w_g(i,j,k) = m_bc_vel_g[bcv][2];
    });
  }

  if (nup > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(uhi_bx_xy, i, j, k,
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        u_g(i,j,k) = u_g(i,j,dom_hi[2]);
      else if(bct == bc_list.minf)
        u_g(i,j,k) = 0;
    });

    AMREX_HOST_DEVICE_FOR_3D(vhi_bx_xy, i, j, k,
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
        v_g(i,j,k) = v_g(i,j,dom_hi[2]);
      else if(bct == bc_list.minf)
        v_g(i,j,k) = 0;
    });

    AMREX_HOST_DEVICE_FOR_3D(whi_bx_xy, i, j, k,
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if(((bct == bc_list.pinf) or (bct == bc_list.pout)) and (k != dom_hi[2]+1))
        w_g(i,j,k) = w_g(i,j,dom_hi[2]+1);
      else if(bct == bc_list.minf)
        w_g(i,j,k) = m_bc_vel_g[bcv][2];
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

}
