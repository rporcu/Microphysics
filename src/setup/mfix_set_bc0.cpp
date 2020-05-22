#include <mfix.H>

#include <MFIX_FLUID_Parms.H>
#include <MFIX_SPECIES_Parms.H>
#include <MFIX_MFHelpers.H>

void
mfix::set_bc0 (const Box& sbx,
               MFIter* mfi,
               const int lev,
               const Box& domain)
{
  // TODO the following won't be constant but computed as functions of the other
  // variables
  //const Real T_g0   = FLUID::T_g0;
  const Real ro_g0  = FLUID::ro_g0;
  const Real cp_g0  = FLUID::cp_g0;
  const Real k_g0   = FLUID::k_g0;
  const Real mu_g0  = FLUID::mu_g0;
  const Real trac_0 = FLUID::trac_0;

  //const Real h_g0   = FLUID::T_g0 * FLUID::cp_g0;

  const int nspecies_g = FLUID::nspecies_g;

  Gpu::ManagedVector< Real > D_g0_managed(nspecies_g, 0.);

  for (int n(0); n < nspecies_g; n++) {
    D_g0_managed[n] = FLUID::D_g0[n];
  }

  Real* D_g0 = D_g0_managed.data();

  Array4<Real> const& a_ep_g = m_leveldata[lev]->ep_g->array(*mfi);
  Array4<Real> const& a_h_g  = m_leveldata[lev]->h_g->array(*mfi);
  Array4<Real> const& a_T_g  = m_leveldata[lev]->T_g->array(*mfi);
  Array4<Real> const& a_ro_g = m_leveldata[lev]->ro_g->array(*mfi);
  Array4<Real> const& a_trac = m_leveldata[lev]->trac->array(*mfi);
  Array4<Real> const& a_cp_g = m_leveldata[lev]->cp_g->array(*mfi);
  Array4<Real> const& a_k_g  = m_leveldata[lev]->k_g->array(*mfi);
  Array4<Real> const& a_mu_g = m_leveldata[lev]->mu_g->array(*mfi);
 
  // The following auxiliary MultiFabs are helpful when we do not solve for
  // species mass fractions
  MultiFab* X_g_aux;
  MultiFab* D_g_aux;
  if (FLUID::nspecies_g == 0) {
    X_g_aux = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0).release();
    D_g_aux = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0).release();
  }
  else {
    X_g_aux = m_leveldata[lev]->X_g;
    D_g_aux = m_leveldata[lev]->D_g;
  }

  Array4<Real> const& a_X_g  = X_g_aux->array(*mfi);
  Array4<Real> const& a_D_g  = D_g_aux->array(*mfi);

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

  amrex::Real* p_bc_ep_g = m_bc_ep_g.data();
  amrex::Real* p_bc_t_g = m_bc_t_g.data();
  
  Gpu::ManagedVector< Real* > m_bc_X_g_managed(nspecies_g);
  
  if (advect_fluid_species) {
    for (int n(0); n < nspecies_g; n++)
      m_bc_X_g_managed[n] = m_bc_X_g[n].data();
  }

  Real** p_bc_X_g = m_bc_X_g_managed.data();

  if (nlft > 0)
  {
    IntVect bx_yz_lo_hi_3D(sbx_hi);

    // Fix lo and hi limits
    bx_yz_lo_hi_3D[0] = dom_lo[0]-1;

    const Box bx_yz_lo_3D(sbx_lo, bx_yz_lo_hi_3D);

    ParallelFor(bx_yz_lo_3D, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_ilo(dom_lo[0]-1,j,k,1);
      const int bct = a_bc_ilo(dom_lo[0]-1,j,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        Real bc_ro_g(ro_g0);
        Real bc_trac(trac_0);
        Real bc_cp_g(cp_g0);
        Real bc_k_g(k_g0);
        Real bc_mu_g(mu_g0);

        a_ep_g(i,j,k) = p_bc_ep_g[bcv];
        a_ro_g(i,j,k) = bc_ro_g;
        a_trac(i,j,k) = bc_trac;
        a_mu_g(i,j,k) = bc_mu_g;
        a_T_g(i,j,k)  = p_bc_t_g[bcv];
        a_cp_g(i,j,k) = bc_cp_g;
        a_k_g(i,j,k)  = bc_k_g;
        a_h_g(i,j,k)  = a_cp_g(i,j,k)*p_bc_t_g[bcv];

        for (int n(0); n < nspecies_g; n++) {
          a_D_g(i,j,k,n) = D_g0[n];
          a_X_g(i,j,k,n) = p_bc_X_g[n][bcv];
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

    ParallelFor(bx_yz_hi_3D, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_ihi(dom_hi[0]+1,j,k,1);
      const int bct = a_bc_ihi(dom_hi[0]+1,j,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        Real bc_ro_g(ro_g0);
        Real bc_trac(trac_0);
        Real bc_cp_g(cp_g0);
        Real bc_k_g(k_g0);
        Real bc_mu_g(mu_g0);

        a_ep_g(i,j,k) = p_bc_ep_g[bcv];
        a_ro_g(i,j,k) = bc_ro_g;
        a_trac(i,j,k) = bc_trac;
        a_mu_g(i,j,k) = bc_mu_g;
        a_T_g(i,j,k)  = p_bc_t_g[bcv];
        a_cp_g(i,j,k) = bc_cp_g;
        a_k_g(i,j,k)  = bc_k_g;
        a_h_g(i,j,k)  = a_cp_g(i,j,k)*p_bc_t_g[bcv];

        for (int n(0); n < nspecies_g; n++) {
          a_D_g(i,j,k,n) = D_g0[n];
          a_X_g(i,j,k,n) = p_bc_X_g[n][bcv];
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

    ParallelFor(bx_xz_lo_3D, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_jlo(i,dom_lo[1]-1,k,1);
      const int bct = a_bc_jlo(i,dom_lo[1]-1,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        Real bc_ro_g(ro_g0);
        Real bc_trac(trac_0);
        Real bc_cp_g(cp_g0);
        Real bc_k_g(k_g0);
        Real bc_mu_g(mu_g0);

        a_ep_g(i,j,k) = p_bc_ep_g[bcv];
        a_ro_g(i,j,k) = bc_ro_g;
        a_trac(i,j,k) = bc_trac;
        a_mu_g(i,j,k) = bc_mu_g;
        a_T_g(i,j,k)  = p_bc_t_g[bcv];
        a_cp_g(i,j,k) = bc_cp_g;
        a_k_g(i,j,k)  = bc_k_g;
        a_h_g(i,j,k)  = a_cp_g(i,j,k)*p_bc_t_g[bcv];

        for (int n(0); n < nspecies_g; n++) {
          a_D_g(i,j,k,n) = D_g0[n];
          a_X_g(i,j,k,n) = p_bc_X_g[n][bcv];
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

    ParallelFor(bx_xz_hi_3D, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_jhi(i,dom_hi[1]+1,k,1);
      const int bct = a_bc_jhi(i,dom_hi[1]+1,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        Real bc_ro_g(ro_g0);
        Real bc_trac(trac_0);
        Real bc_cp_g(cp_g0);
        Real bc_k_g(k_g0);
        Real bc_mu_g(mu_g0);

        a_ep_g(i,j,k) = p_bc_ep_g[bcv];
        a_ro_g(i,j,k) = bc_ro_g;
        a_trac(i,j,k) = bc_trac;
        a_mu_g(i,j,k) = bc_mu_g;
        a_T_g(i,j,k)  = p_bc_t_g[bcv];
        a_cp_g(i,j,k) = bc_cp_g;
        a_k_g(i,j,k)  = bc_k_g;
        a_h_g(i,j,k)  = a_cp_g(i,j,k)*p_bc_t_g[bcv];

        for (int n(0); n < nspecies_g; n++) {
          a_D_g(i,j,k,n) = D_g0[n];
          a_X_g(i,j,k,n) = p_bc_X_g[n][bcv];
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

    ParallelFor(bx_xy_lo_3D, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_klo(i,j,dom_lo[2]-1,1);
      const int bct = a_bc_klo(i,j,dom_lo[2]-1,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        Real bc_ro_g(ro_g0);
        Real bc_trac(trac_0);
        Real bc_cp_g(cp_g0);
        Real bc_k_g(k_g0);
        Real bc_mu_g(mu_g0);

        a_ep_g(i,j,k) = p_bc_ep_g[bcv];
        a_ro_g(i,j,k) = bc_ro_g;
        a_trac(i,j,k) = bc_trac;
        a_mu_g(i,j,k) = bc_mu_g;
        a_T_g(i,j,k)  = p_bc_t_g[bcv];
        a_cp_g(i,j,k) = bc_cp_g;
        a_k_g(i,j,k)  = bc_k_g;
        a_h_g(i,j,k)  = a_cp_g(i,j,k)*p_bc_t_g[bcv];

        for (int n(0); n < nspecies_g; n++) {
          a_D_g(i,j,k,n) = D_g0[n];
          a_X_g(i,j,k,n) = p_bc_X_g[n][bcv];
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

    ParallelFor(bx_xy_hi_3D, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_khi(i,j,dom_hi[2]+1,1);
      const int bct = a_bc_khi(i,j,dom_hi[2]+1,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        Real bc_ro_g(ro_g0);
        Real bc_trac(trac_0);
        Real bc_cp_g(cp_g0);
        Real bc_k_g(k_g0);
        Real bc_mu_g(mu_g0);

        a_ep_g(i,j,k) = p_bc_ep_g[bcv];
        a_ro_g(i,j,k) = bc_ro_g;
        a_trac(i,j,k) = bc_trac;
        a_mu_g(i,j,k) = bc_mu_g;
        a_T_g(i,j,k)  = p_bc_t_g[bcv];
        a_cp_g(i,j,k) = bc_cp_g;
        a_k_g(i,j,k)  = bc_k_g;
        a_h_g(i,j,k)  = a_cp_g(i,j,k)*p_bc_t_g[bcv];

        for (int n(0); n < nspecies_g; n++) {
          a_D_g(i,j,k,n) = D_g0[n];
          a_X_g(i,j,k,n) = p_bc_X_g[n][bcv];
        }
      }
    });
  }

  Gpu::synchronize();

  if (FLUID::nspecies_g == 0) {
    delete X_g_aux;
    delete D_g_aux;
  }

  if (advect_fluid_species) {
    delete[] D_g0;
  }
}
