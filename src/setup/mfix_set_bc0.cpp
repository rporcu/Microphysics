#include <mfix.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_mf_helpers.H>

void
mfix::set_bc0 (const Box& sbx,
               MFIter* mfi,
               const int lev,
               const Box& domain)
{
  const Real ro_g0  = FLUID::ro_g0;
  const Real mu_g0  = FLUID::mu_g0;
  const Real trac_0 = FLUID::trac_0;
  const Real MW_g0  = FLUID::MW_g0;

  const int nspecies_g = FLUID::nspecies;

  Gpu::ManagedVector< Real > MW_gk0_managed(nspecies_g);
  Gpu::ManagedVector< Real* > m_bc_X_gk_managed(nspecies_g);

  for (int n(0); n < nspecies_g; n++) {
    MW_gk0_managed[n] = FLUID::MW_gk0[n];
    m_bc_X_gk_managed[n] = m_bc_X_gk[n].data();
  }

  // Flag to understand if fluid is a mixture
  const int fluid_is_a_mixture = FLUID::is_a_mixture;

  Real* p_MW_gk0 = fluid_is_a_mixture ? MW_gk0_managed.data() : nullptr;
  Real** p_bc_X_gk = fluid_is_a_mixture ? m_bc_X_gk_managed.data() : nullptr;

  amrex::Real* p_bc_ep_g = m_bc_ep_g.data();

  Array4<Real> const& a_ep_g = m_leveldata[lev]->ep_g->array(*mfi);
  Array4<Real> const& a_ro_g = m_leveldata[lev]->ro_g->array(*mfi);
  Array4<Real> const& a_trac = m_leveldata[lev]->trac->array(*mfi);
  Array4<Real> const& a_mu_g = m_leveldata[lev]->mu_g->array(*mfi);
  Array4<Real> const& a_MW_g = m_leveldata[lev]->MW_g->array(*mfi);
 
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

  const int nlft = amrex::max(0,dom_lo[0]-sbx_lo[0]);
  const int nbot = amrex::max(0,dom_lo[1]-sbx_lo[1]);
  const int ndwn = amrex::max(0,dom_lo[2]-sbx_lo[2]);

  const int nrgt = amrex::max(0,sbx_hi[0]-dom_hi[0]);
  const int ntop = amrex::max(0,sbx_hi[1]-dom_hi[1]);
  const int nup  = amrex::max(0,sbx_hi[2]-dom_hi[2]);

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  if (nlft > 0)
  {
    IntVect bx_yz_lo_hi_3D(sbx_hi);

    // Fix lo and hi limits
    bx_yz_lo_hi_3D[0] = dom_lo[0]-1;

    const Box bx_yz_lo_3D(sbx_lo, bx_yz_lo_hi_3D);

    ParallelFor(bx_yz_lo_3D, [a_bc_ilo,minf,pinf,pout,a_ep_g,a_ro_g,a_trac,
        a_mu_g,p_bc_ep_g,ro_g0,trac_0,mu_g0,a_MW_g,MW_g0,p_bc_X_gk,p_MW_gk0,
        nspecies_g,dom_lo,fluid_is_a_mixture]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_ilo(dom_lo[0]-1,j,k,1);
      const int bct = a_bc_ilo(dom_lo[0]-1,j,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        a_ep_g(i,j,k) = p_bc_ep_g[bcv];
        a_ro_g(i,j,k) = ro_g0;
        a_trac(i,j,k) = trac_0;
        a_mu_g(i,j,k) = mu_g0;

        if (not fluid_is_a_mixture) {
          a_MW_g(i,j,k) = MW_g0;
        }
        else {
          Real MW_gk_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_gk_sum += p_bc_X_gk[n][bcv]/p_MW_gk0[n];
          }

          a_MW_g(i,j,k) = 1./MW_gk_sum;
        }
      }
    });
  }

  if (nrgt > 0)
  {
    IntVect bx_yz_hi_lo_3D(sbx_lo);

    // Fix lo and hi limits
    bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

    const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, sbx_hi);

    ParallelFor(bx_yz_hi_3D, [a_bc_ihi,minf,pinf,pout,a_ep_g,a_ro_g,a_trac,
        a_mu_g,p_bc_ep_g,ro_g0,trac_0,mu_g0,a_MW_g,MW_g0,p_bc_X_gk,p_MW_gk0,
        nspecies_g,dom_hi,fluid_is_a_mixture]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_ihi(dom_hi[0]+1,j,k,1);
      const int bct = a_bc_ihi(dom_hi[0]+1,j,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        a_ep_g(i,j,k) = p_bc_ep_g[bcv];
        a_ro_g(i,j,k) = ro_g0;
        a_trac(i,j,k) = trac_0;
        a_mu_g(i,j,k) = mu_g0;

        if (not fluid_is_a_mixture) {
          a_MW_g(i,j,k) = MW_g0;
        }
        else {
          Real MW_gk_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_gk_sum += p_bc_X_gk[n][bcv]/p_MW_gk0[n];
          }

          a_MW_g(i,j,k) = 1./MW_gk_sum;
        }
      }
    });
  }

  if (nbot > 0)
  {
    IntVect bx_xz_lo_hi_3D(sbx_hi);

    // Fix lo and hi limits
    bx_xz_lo_hi_3D[1] = dom_lo[1]-1;

    const Box bx_xz_lo_3D(sbx_lo, bx_xz_lo_hi_3D);

    ParallelFor(bx_xz_lo_3D, [a_bc_jlo,minf,pinf,pout,a_ep_g,a_ro_g,a_trac,
        a_mu_g,p_bc_ep_g,ro_g0,trac_0,mu_g0,a_MW_g,MW_g0,p_bc_X_gk,p_MW_gk0,
        nspecies_g,dom_lo,fluid_is_a_mixture]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_jlo(i,dom_lo[1]-1,k,1);
      const int bct = a_bc_jlo(i,dom_lo[1]-1,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        a_ep_g(i,j,k) = p_bc_ep_g[bcv];
        a_ro_g(i,j,k) = ro_g0;
        a_trac(i,j,k) = trac_0;
        a_mu_g(i,j,k) = mu_g0;

        if (not fluid_is_a_mixture) {
          a_MW_g(i,j,k) = MW_g0;
        }
        else {
          Real MW_gk_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_gk_sum += p_bc_X_gk[n][bcv]/p_MW_gk0[n];
          }

          a_MW_g(i,j,k) = 1./MW_gk_sum;
        }
      }
    });
  }

  if (ntop > 0)
  {
    IntVect bx_xz_hi_lo_3D(sbx_lo);

    // Fix lo and hi limits
    bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

    const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, sbx_hi);

    ParallelFor(bx_xz_hi_3D, [a_bc_jhi,minf,pinf,pout,a_ep_g,a_ro_g,a_trac,
        a_mu_g,p_bc_ep_g,ro_g0,trac_0,mu_g0,a_MW_g,MW_g0,p_bc_X_gk,p_MW_gk0,
        nspecies_g,dom_hi,fluid_is_a_mixture]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_jhi(i,dom_hi[1]+1,k,1);
      const int bct = a_bc_jhi(i,dom_hi[1]+1,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        a_ep_g(i,j,k) = p_bc_ep_g[bcv];
        a_ro_g(i,j,k) = ro_g0;
        a_trac(i,j,k) = trac_0;
        a_mu_g(i,j,k) = mu_g0;

        if (not fluid_is_a_mixture) {
          a_MW_g(i,j,k) = MW_g0;
        }
        else {
          Real MW_gk_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_gk_sum += p_bc_X_gk[n][bcv]/p_MW_gk0[n];
          }

          a_MW_g(i,j,k) = 1./MW_gk_sum;
        }
      }
    });
  }

  if (ndwn > 0)
  {
    IntVect bx_xy_lo_hi_3D(sbx_hi);

    // Fix lo and hi limits
    bx_xy_lo_hi_3D[2] = dom_lo[2]-1;

    const Box bx_xy_lo_3D(sbx_lo, bx_xy_lo_hi_3D);

    ParallelFor(bx_xy_lo_3D, [a_bc_klo,minf,pinf,pout,a_ep_g,a_ro_g,a_trac,
        a_mu_g,p_bc_ep_g,ro_g0,trac_0,mu_g0,a_MW_g,MW_g0,p_bc_X_gk,p_MW_gk0,
        nspecies_g,dom_lo,fluid_is_a_mixture]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_klo(i,j,dom_lo[2]-1,1);
      const int bct = a_bc_klo(i,j,dom_lo[2]-1,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        a_ep_g(i,j,k) = p_bc_ep_g[bcv];
        a_ro_g(i,j,k) = ro_g0;
        a_trac(i,j,k) = trac_0;
        a_mu_g(i,j,k) = mu_g0;

        if (not fluid_is_a_mixture) {
          a_MW_g(i,j,k) = MW_g0;
        }
        else {
          Real MW_gk_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_gk_sum += p_bc_X_gk[n][bcv]/p_MW_gk0[n];
          }

          a_MW_g(i,j,k) = 1./MW_gk_sum;
        }
      }
    });
  }

  if (nup > 0)
  {
    IntVect bx_xy_hi_lo_3D(sbx_lo);

    // Fix lo and hi limits
    bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

    const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, sbx_hi);

    ParallelFor(bx_xy_hi_3D, [a_bc_khi,minf,pinf,pout,a_ep_g,a_ro_g,a_trac,
        a_mu_g,p_bc_ep_g,ro_g0,trac_0,mu_g0,a_MW_g,MW_g0,p_bc_X_gk,p_MW_gk0,
        nspecies_g,dom_hi,fluid_is_a_mixture]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_khi(i,j,dom_hi[2]+1,1);
      const int bct = a_bc_khi(i,j,dom_hi[2]+1,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        a_ep_g(i,j,k) = p_bc_ep_g[bcv];
        a_ro_g(i,j,k) = ro_g0;
        a_trac(i,j,k) = trac_0;
        a_mu_g(i,j,k) = mu_g0;

        if (not fluid_is_a_mixture) {
          a_MW_g(i,j,k) = MW_g0;
        }
        else {
          Real MW_gk_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_gk_sum += p_bc_X_gk[n][bcv]/p_MW_gk0[n];
          }

          a_MW_g(i,j,k) = 1./MW_gk_sum;
        }
      }
    });
  }

}
