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

#include <mfix.H>
#include <bc_mod_F.H>
#include <eos_mod.hpp>
#include <fld_constants_mod_F.H>
#include <param_mod_F.H>

using namespace amrex;

void 
mfix::set_scalar_bcs(MFIter* mfi,
                     const int lev,
                     const Box& domain)
{
  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<amrex::Real> const& ep_g_a = ep_g[lev]->array(*mfi);
  Array4<amrex::Real> const& ro_g_a = ro_g[lev]->array(*mfi);
  Array4<amrex::Real> const& mu_g_a = mu_g[lev]->array(*mfi);

  IntVect ep_g_lo((*ep_g[lev])[*mfi].loVect());
  IntVect ep_g_hi((*ep_g[lev])[*mfi].hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

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

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  amrex::Real* p_bc_ep_g = m_bc_ep_g.data();
  amrex::Real* p_bc_t_g = m_bc_t_g.data();

  if (nlft > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_yz_lo_3D, i, j, k,
    {
      Real bc_ro_g(ro_g0);
      Real bc_mu_g(0);

      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == pinf) or (bct == pout))
      {
        ep_g_a(i,j,k) = ep_g_a(dom_lo[0],j,k);
        ro_g_a(i,j,k) = ro_g_a(dom_lo[0],j,k);
        mu_g_a(i,j,k) = mu_g_a(dom_lo[0],j,k);
      }
      else if(bct == minf)
      {
        if(is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(p_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        ro_g_a(i,j,k) = bc_ro_g;
        mu_g_a(i,j,k) = bc_mu_g;
      }
    });

    AMREX_HOST_DEVICE_FOR_3D(bx_yz_lo_2D, i, j, k,
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if(bct == minf)
        ep_g_a(i,j,k) = 2*p_bc_ep_g[bcv] - ep_g_a(i+1,j,k);
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

      if((bct == pinf) or (bct == pout))
      {
        ep_g_a(i,j,k) = ep_g_a(dom_hi[0],j,k);
        ro_g_a(i,j,k) = ro_g_a(dom_hi[0],j,k);
        mu_g_a(i,j,k) = mu_g_a(dom_hi[0],j,k);
      }
      else if(bct == minf)
      {
        if(is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(p_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        ro_g_a(i,j,k) = bc_ro_g;
        mu_g_a(i,j,k) = bc_mu_g;
      }
    });

    AMREX_HOST_DEVICE_FOR_3D(bx_yz_hi_2D, i, j, k,
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if(bct == minf)
        ep_g_a(i,j,k) = 2*p_bc_ep_g[bcv] - ep_g_a(i-1,j,k);
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

      if((bct == pinf) or (bct == pout))
      {
        ep_g_a(i,j,k) = ep_g_a(i,dom_lo[1],k);
        ro_g_a(i,j,k) = ro_g_a(i,dom_lo[1],k);
        mu_g_a(i,j,k) = mu_g_a(i,dom_lo[1],k);
      }
      else if(bct == minf)
      {
        if(is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(p_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        ro_g_a(i,j,k) = bc_ro_g;
        mu_g_a(i,j,k) = bc_mu_g;
      }
    });

    AMREX_HOST_DEVICE_FOR_3D(bx_xz_lo_2D, i, j, k,
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if(bct == minf)
        ep_g_a(i,j,k) = 2*p_bc_ep_g[bcv] - ep_g_a(i,j+1,k);
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

      if((bct == pinf) or (bct == pout))
      {
        ep_g_a(i,j,k) = ep_g_a(i,dom_hi[1],k);
        ro_g_a(i,j,k) = ro_g_a(i,dom_hi[1],k);
        mu_g_a(i,j,k) = mu_g_a(i,dom_hi[1],k);
      }
      else if(bct == minf)
      {
        if(is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(p_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        ro_g_a(i,j,k) = bc_ro_g;
        mu_g_a(i,j,k) = bc_mu_g;
      }
    });

    AMREX_HOST_DEVICE_FOR_3D(bx_xz_hi_2D, i, j, k,
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if(bct == minf)
        ep_g_a(i,j,k) = 2*p_bc_ep_g[bcv] - ep_g_a(i,j-1,k);
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

      if((bct == pinf) or (bct == pout))
      {
        ep_g_a(i,j,k) = ep_g_a(i,j,dom_lo[2]);
        ro_g_a(i,j,k) = ro_g_a(i,j,dom_lo[2]);
        mu_g_a(i,j,k) = mu_g_a(i,j,dom_lo[2]);
      }
      else if(bct == minf)
      {
        if(is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(p_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        ro_g_a(i,j,k) = bc_ro_g;
        mu_g_a(i,j,k) = bc_mu_g;
      }
    });

    AMREX_HOST_DEVICE_FOR_3D(bx_xy_lo_2D, i, j, k,
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if(bct == minf)
        ep_g_a(i,j,k) = 2*p_bc_ep_g[bcv] - ep_g_a(i,j,k+1);
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

      if((bct == pinf) or (bct == pout))
      {
        ep_g_a(i,j,k) = ep_g_a(i,j,dom_hi[2]);
        ro_g_a(i,j,k) = ro_g_a(i,j,dom_hi[2]);
        mu_g_a(i,j,k) = mu_g_a(i,j,dom_hi[2]);
      }
      else if(bct == minf)
      {
        if(is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(p_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        ro_g_a(i,j,k) = bc_ro_g;
        mu_g_a(i,j,k) = bc_mu_g;
      }
    });

    AMREX_HOST_DEVICE_FOR_3D(bx_xy_hi_2D, i, j, k,
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if(bct == minf)
        ep_g_a(i,j,k) = 2*p_bc_ep_g[bcv] - ep_g_a(i,j,k-1);
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif
}
