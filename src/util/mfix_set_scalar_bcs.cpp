//
//
//  This subroutine sets the BCs for all scalar variables involved in
//  the projection, EXCEPT the pressure and velocity components.
//  To set velocity BCs, use "set_velocity_bcs".
//
//  Author: Michele Rosso
//
//  Date: December 20, 2017
//
//

#include <mfix_set_scalar_bcs.hpp>
#include <bc_mod_F.H>
#include <eos_mod.hpp>
#include <fld_constants_mod_F.H>
#include <param_mod_F.H>

using namespace amrex;

void set_scalar_bcs(const BcList& bc_list,
                    FArrayBox& ep_g_fab,
                    FArrayBox& ro_g_fab,
                    FArrayBox& mu_g_fab,
                    const IArrayBox& bct_ilo_fab,
                    const IArrayBox& bct_ihi_fab,
                    const IArrayBox& bct_jlo_fab,
                    const IArrayBox& bct_jhi_fab,
                    const IArrayBox& bct_klo_fab,
                    const IArrayBox& bct_khi_fab,
                    const Box& domain,
                    Real* m_bc_ep_g,
                    Real* m_bc_t_g,
                    const int* ng)
{
  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<Real> const& ep_g = ep_g_fab.array();
  Array4<Real> const& ro_g = ro_g_fab.array();
  Array4<Real> const& mu_g = mu_g_fab.array();

  IntVect ep_g_lo(ep_g_fab.loVect());
  IntVect ep_g_hi(ep_g_fab.hiVect());

  Array4<const int> const& bct_ilo = bct_ilo_fab.array();
  Array4<const int> const& bct_ihi = bct_ihi_fab.array();
  Array4<const int> const& bct_jlo = bct_jlo_fab.array();
  Array4<const int> const& bct_jhi = bct_jhi_fab.array();
  Array4<const int> const& bct_klo = bct_klo_fab.array();
  Array4<const int> const& bct_khi = bct_khi_fab.array();

  const int nlft = std::max(0, dom_lo[0]-ep_g_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-ep_g_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-ep_g_lo[2]);

  const int nrgt = std::max(0, ep_g_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, ep_g_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, ep_g_hi[2]-dom_hi[2]);

  const Real ro_g0 = get_ro_g0();
  const Real mu_g0 = get_mu_g0();

  // Create InVects for following 2D Boxes
  IntVect bx_yz_lo_lo_2D(ep_g_lo), bx_yz_lo_hi_2D(ep_g_hi);
  IntVect bx_yz_hi_lo_2D(ep_g_lo), bx_yz_hi_hi_2D(ep_g_hi);
  IntVect bx_xz_lo_lo_2D(ep_g_lo), bx_xz_lo_hi_2D(ep_g_hi);
  IntVect bx_xz_hi_lo_2D(ep_g_lo), bx_xz_hi_hi_2D(ep_g_hi);
  IntVect bx_xy_lo_lo_2D(ep_g_lo), bx_xy_lo_hi_2D(ep_g_hi);
  IntVect bx_xy_hi_lo_2D(ep_g_lo), bx_xy_hi_hi_2D(ep_g_hi);

  // Fix lo and hi limits
  bx_yz_lo_lo_2D[0] = dom_lo[0]-1;
  bx_yz_lo_hi_2D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_2D[0] = dom_hi[0]+1;
  bx_yz_hi_hi_2D[0] = dom_hi[0]+1;

  bx_xz_lo_lo_2D[1] = dom_lo[1]-1;
  bx_xz_lo_hi_2D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_2D[1] = dom_hi[1]+1;
  bx_xz_hi_hi_2D[1] = dom_hi[1]+1;

  bx_xy_lo_lo_2D[2] = dom_lo[2]-1;
  bx_xy_lo_hi_2D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_2D[2] = dom_hi[2]+1;
  bx_xy_hi_hi_2D[2] = dom_hi[2]+1;

  // Create 2D boxes for CUDA loops
  const Box bx_yz_lo_2D(bx_yz_lo_lo_2D, bx_yz_lo_hi_2D);
  const Box bx_yz_hi_2D(bx_yz_hi_lo_2D, bx_yz_hi_hi_2D);

  const Box bx_xz_lo_2D(bx_xz_lo_lo_2D, bx_xz_lo_hi_2D);
  const Box bx_xz_hi_2D(bx_xz_hi_lo_2D, bx_xz_hi_hi_2D);

  const Box bx_xy_lo_2D(bx_xy_lo_lo_2D, bx_xy_lo_hi_2D);
  const Box bx_xy_hi_2D(bx_xy_hi_lo_2D, bx_xy_hi_hi_2D);

  // Create InVects for following 3D Boxes
  IntVect bx_yz_lo_hi_3D(ep_g_hi), bx_xz_lo_hi_3D(ep_g_hi), bx_xy_lo_hi_3D(ep_g_hi);
  IntVect bx_yz_hi_lo_3D(ep_g_lo), bx_xz_hi_lo_3D(ep_g_lo), bx_xy_hi_lo_3D(ep_g_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for CUDA loops
  const Box bx_yz_lo_3D(ep_g_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, ep_g_hi);

  const Box bx_xz_lo_3D(ep_g_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, ep_g_hi);

  const Box bx_xy_lo_3D(ep_g_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, ep_g_hi);

  const Real undefined = get_undefined();

  if (nlft > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_yz_lo_3D, i, j, k,
    {
      Real bc_ro_g(ro_g0);
      Real bc_mu_g(0);

      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        ep_g(i,j,k) = ep_g(dom_lo[0],j,k);
        ro_g(i,j,k) = ro_g(dom_lo[0],j,k);
        mu_g(i,j,k) = mu_g(dom_lo[0],j,k);
      }
      else if(bct == bc_list.minf)
      {
        if(is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(m_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        ro_g(i,j,k) = bc_ro_g;
        mu_g(i,j,k) = bc_mu_g;
      }
    });

    AMREX_HOST_DEVICE_FOR_3D(bx_yz_lo_2D, i, j, k,
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if(bct == bc_list.minf)
        ep_g(i,j,k) = 2*m_bc_ep_g[bcv] - ep_g(i+1,j,k);
    });
  }

  if (nrgt > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_yz_hi_3D, i, j, k,
    {
      Real bc_ro_g(ro_g0);
      Real bc_mu_g(0);

      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        ep_g(i,j,k) = ep_g(dom_hi[0],j,k);
        ro_g(i,j,k) = ro_g(dom_hi[0],j,k);
        mu_g(i,j,k) = mu_g(dom_hi[0],j,k);
      }
      else if(bct == bc_list.minf)
      {
        if(is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(m_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        ro_g(i,j,k) = bc_ro_g;
        mu_g(i,j,k) = bc_mu_g;
      }
    });

    AMREX_HOST_DEVICE_FOR_3D(bx_yz_hi_2D, i, j, k,
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if(bct == bc_list.minf)
        ep_g(i,j,k) = 2*m_bc_ep_g[bcv] - ep_g(i-1,j,k);
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (nbot > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xz_lo_3D, i, j, k,
    {
      Real bc_ro_g(ro_g0);
      Real bc_mu_g(0);

      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        ep_g(i,j,k) = ep_g(i,dom_lo[1],k);
        ro_g(i,j,k) = ro_g(i,dom_lo[1],k);
        mu_g(i,j,k) = mu_g(i,dom_lo[1],k);
      }
      else if(bct == bc_list.minf)
      {
        if(is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(m_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        ro_g(i,j,k) = bc_ro_g;
        mu_g(i,j,k) = bc_mu_g;
      }
    });

    AMREX_HOST_DEVICE_FOR_3D(bx_xz_lo_2D, i, j, k,
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if(bct == bc_list.minf)
        ep_g(i,j,k) = 2*m_bc_ep_g[bcv] - ep_g(i,j+1,k);
    });
  }

  if (ntop > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xz_hi_3D, i, j, k,
    {
      Real bc_ro_g(ro_g0);
      Real bc_mu_g(0);

      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        ep_g(i,j,k) = ep_g(i,dom_hi[1],k);
        ro_g(i,j,k) = ro_g(i,dom_hi[1],k);
        mu_g(i,j,k) = mu_g(i,dom_hi[1],k);
      }
      else if(bct == bc_list.minf)
      {
        if(is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(m_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        ro_g(i,j,k) = bc_ro_g;
        mu_g(i,j,k) = bc_mu_g;
      }
    });

    AMREX_HOST_DEVICE_FOR_3D(bx_xz_hi_2D, i, j, k,
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if(bct == bc_list.minf)
        ep_g(i,j,k) = 2*m_bc_ep_g[bcv] - ep_g(i,j-1,k);
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (ndwn > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xy_lo_3D, i, j, k,
    {
      Real bc_ro_g(ro_g0);
      Real bc_mu_g(0);

      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        ep_g(i,j,k) = ep_g(i,j,dom_lo[2]);
        ro_g(i,j,k) = ro_g(i,j,dom_lo[2]);
        mu_g(i,j,k) = mu_g(i,j,dom_lo[2]);
      }
      else if(bct == bc_list.minf)
      {
        if(is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(m_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        ro_g(i,j,k) = bc_ro_g;
        mu_g(i,j,k) = bc_mu_g;
      }
    });

    AMREX_HOST_DEVICE_FOR_3D(bx_xy_lo_2D, i, j, k,
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if(bct == bc_list.minf)
        ep_g(i,j,k) = 2*m_bc_ep_g[bcv] - ep_g(i,j,k+1);
    });
  }

  if (nup > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xy_hi_3D, i, j, k,
    {
      Real bc_ro_g(ro_g0);
      Real bc_mu_g(0);

      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout))
      {
        ep_g(i,j,k) = ep_g(i,j,dom_hi[2]);
        ro_g(i,j,k) = ro_g(i,j,dom_hi[2]);
        mu_g(i,j,k) = mu_g(i,j,dom_hi[2]);
      }
      else if(bct == bc_list.minf)
      {
        if(is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(m_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        ro_g(i,j,k) = bc_ro_g;
        mu_g(i,j,k) = bc_mu_g;
      }
    });

    AMREX_HOST_DEVICE_FOR_3D(bx_xy_hi_2D, i, j, k,
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if(bct == bc_list.minf)
        ep_g(i,j,k) = 2*m_bc_ep_g[bcv] - ep_g(i,j,k-1);
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif
}
