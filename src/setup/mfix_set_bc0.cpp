#include <mfix.H>
#include <bc_mod_F.H>
#include <eos_mod.hpp>
#include <fld_constants_mod_F.H>
#include <param_mod_F.H>

void 
mfix::set_bc0(const Box& sbx,
              MFIter* mfi,
              const int lev,
              const Box& domain)
{
  const Real ro_g0(get_ro_g0());
  const Real mu_g0(get_mu_g0());

  Array4<Real> const& a_ep_g = ep_g[lev]->array(*mfi);
  Array4<Real> const& a_ro_g = ro_g[lev]->array(*mfi);
  Array4<Real> const& a_mu_g = mu_g[lev]->array(*mfi);

  const IntVect sbx_lo(sbx.loVect());
  const IntVect sbx_hi(sbx.hiVect());

  const IntVect dom_lo(domain.loVect());
  const IntVect dom_hi(domain.hiVect());

  Array4<const int> const& a_bc_ilo = bc_ilo[lev]->array();
  Array4<const int> const& a_bc_ihi = bc_ihi[lev]->array();
  Array4<const int> const& a_bc_jlo = bc_jlo[lev]->array();
  Array4<const int> const& a_bc_jhi = bc_jhi[lev]->array();
  Array4<const int> const& a_bc_klo = bc_klo[lev]->array();
  Array4<const int> const& a_bc_khi = bc_khi[lev]->array();

  const int nlft = std::max(0,dom_lo[0]-sbx_lo[0]);
  const int nbot = std::max(0,dom_lo[1]-sbx_lo[1]);
  const int ndwn = std::max(0,dom_lo[2]-sbx_lo[2]);

  const int nrgt = std::max(0,sbx_hi[0]-dom_hi[0]);
  const int ntop = std::max(0,sbx_hi[1]-dom_hi[1]);
  const int nup  = std::max(0,sbx_hi[2]-dom_hi[2]);

  // Create InVects for following 3D Boxes
  IntVect bx_yz_lo_hi_3D(sbx_hi), bx_xz_lo_hi_3D(sbx_hi), bx_xy_lo_hi_3D(sbx_hi);
  IntVect bx_yz_hi_lo_3D(sbx_lo), bx_xz_hi_lo_3D(sbx_lo), bx_xy_hi_lo_3D(sbx_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for CUDA loops
  const Box bx_yz_lo_3D(sbx_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, sbx_hi);

  const Box bx_xz_lo_3D(sbx_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, sbx_hi);

  const Box bx_xy_lo_3D(sbx_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, sbx_hi);

  const Real undefined = get_undefined();

  if (nlft > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_yz_lo_3D, i, j, k,
    {
      const int bcv = a_bc_ilo(dom_lo[0]-1,j,k,1);
      const int bct = a_bc_ilo(dom_lo[0]-1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout) or (bct == bc_list.minf))
      {
        Real bc_ro_g(ro_g0);
        Real bc_mu_g(0);

        if (is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(m_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        a_ep_g(i,j,k) = m_bc_ep_g[bcv];
        a_ro_g(i,j,k) = bc_ro_g;
        a_mu_g(i,j,k) = bc_mu_g;
      }
    });
  }
  
  if (nrgt > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_yz_hi_3D, i, j, k,
    {
      const int bcv = a_bc_ihi(dom_hi[0]+1,j,k,1);
      const int bct = a_bc_ihi(dom_hi[0]+1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout) or (bct == bc_list.minf))
      {
        Real bc_ro_g(ro_g0);
        Real bc_mu_g(0);

        if (is_equal(mu_g0, undefined))
          bc_mu_g = sutherland(m_bc_t_g[bcv]);
        else
          bc_mu_g = mu_g0;

        a_ep_g(i,j,k) = m_bc_ep_g[bcv];
        a_ro_g(i,j,k) = bc_ro_g;
        a_mu_g(i,j,k) = bc_mu_g;
      }
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (nbot > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xz_lo_3D, i, j, k,
    {
      const int bcv = a_bc_jlo(i,dom_lo[1]-1,k,1);
      const int bct = a_bc_jlo(i,dom_lo[1]-1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout) or (bct == bc_list.minf))
      {
        Real bc_ro_g(ro_g0);
        Real bc_mu_g(0);

        if (is_equal(mu_g0, undefined))
           bc_mu_g = sutherland(m_bc_t_g[bcv]);
        else
           bc_mu_g = mu_g0;

        a_ep_g(i,j,k) = m_bc_ep_g[bcv];
        a_ro_g(i,j,k) = bc_ro_g;
        a_mu_g(i,j,k) = bc_mu_g;
      }
    });
  }

  if (ntop > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xz_hi_3D, i, j, k,
    {
      const int bcv = a_bc_jhi(i,dom_hi[1]+1,k,1);
      const int bct = a_bc_jhi(i,dom_hi[1]+1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout) or (bct == bc_list.minf))
      {
        Real bc_ro_g(ro_g0);
        Real bc_mu_g(0);

        if (is_equal(mu_g0, undefined))
           bc_mu_g = sutherland(m_bc_t_g[bcv]);
        else
           bc_mu_g = mu_g0;

        a_ep_g(i,j,k) = m_bc_ep_g[bcv];
        a_ro_g(i,j,k) = bc_ro_g;
        a_mu_g(i,j,k) = bc_mu_g;
      }
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (ndwn > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xy_lo_3D, i, j, k,
    {
      const int bcv = a_bc_klo(i,j,dom_lo[2]-1,1);
      const int bct = a_bc_klo(i,j,dom_lo[2]-1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout) or (bct == bc_list.minf))
      {
        Real bc_ro_g(ro_g0);
        Real bc_mu_g(0);

        if (is_equal(mu_g0, undefined))
           bc_mu_g = sutherland(m_bc_t_g[bcv]);
        else
           bc_mu_g = mu_g0;

        a_ep_g(i,j,k) = m_bc_ep_g[bcv];
        a_ro_g(i,j,k) = bc_ro_g;
        a_mu_g(i,j,k) = bc_mu_g;
      }
    });
  }

  if (nup > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xy_hi_3D, i, j, k,
    {
      const int bcv = a_bc_khi(i,j,dom_hi[2]+1,1);
      const int bct = a_bc_khi(i,j,dom_hi[2]+1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout) or (bct == bc_list.minf))
      {
        Real bc_ro_g(ro_g0);
        Real bc_mu_g(0);

        if (is_equal(mu_g0, undefined))
           bc_mu_g = sutherland(m_bc_t_g[bcv]);
        else
           bc_mu_g = mu_g0;

        a_ep_g(i,j,k) = m_bc_ep_g[bcv];
        a_ro_g(i,j,k) = bc_ro_g;
        a_mu_g(i,j,k) = bc_mu_g;
      }
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif
}
