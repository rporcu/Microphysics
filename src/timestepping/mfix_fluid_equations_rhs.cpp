#include <mfix.H>

#include <mfix_mf_helpers.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>


void
mfix::mfix_density_rhs (Vector< MultiFab* > const& rhs,
                        Vector< MultiFab* > const& ro_gk_txfr)
{
  if (solve_reactions) {
    for (int lev = 0; lev <= finest_level; lev++) {
      for (int n_g(0); n_g < FLUID::nspecies; n_g++) {
        MultiFab::Add(*rhs[lev], *ro_gk_txfr[lev], n_g, 0, 1, rhs[lev]->nGrow());
      }
    }

    for (int lev = 0; lev <= finest_level; lev++)
      EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
  }
}


void
mfix::mfix_enthalpy_rhs (const bool explicit_diffusion,
                         Vector< MultiFab* > const& rhs,
                         Vector< MultiFab* > const& T_g,
                         Vector< MultiFab* > const& ep_g,
                         Vector< MultiFab* > const& ro_g,
                         Vector< MultiFab* > const& k_g,
                         Vector< MultiFab* > const& T_g_on_eb,
                         Vector< MultiFab* > const& k_g_on_eb,
                         Vector< MultiFab* > const& X_gk,
                         Vector< MultiFab* > const& D_gk,
                         Vector< MultiFab* > const& h_gk)
{
  if (explicit_diffusion) {
    diffusion_op->ComputeLapT(rhs, T_g, ep_g, k_g, T_g_on_eb, k_g_on_eb);

    for (int lev = 0; lev <= finest_level; lev++)
      EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
  }

  if (FLUID::is_a_mixture)
  {
    // Temporary for computing other terms of RHS
    Vector< MultiFab* > h_gk_D_gk(nlev, nullptr);
    Vector< MultiFab* > auxiliary(nlev, nullptr);

    // Allocate memory for computing fluid species contribution
    for (int lev(0); lev <= finest_level; lev++) {
      h_gk_D_gk[lev] = MFHelpers::createFrom(*D_gk[lev]).release();
      auxiliary[lev] = MFHelpers::createFrom(*X_gk[lev], 0.).release();
    }

    // Transform h_gk_D_gk into h_gk * D_gk
    for (int lev(0); lev <= finest_level; lev++) {
      MultiFab::Multiply(*h_gk_D_gk[lev], *h_gk[lev], 0, 0, FLUID::nspecies,
          h_gk_D_gk[lev]->nGrow());
    }

    // Compute the mixed enthalpy/species term
    diffusion_op->ComputeLapX(auxiliary, X_gk, ro_g, ep_g, h_gk_D_gk);

    for (int lev(0); lev <= finest_level; lev++) {
      // Add the contribution due to the nth specie
      for (int n(0); n < FLUID::nspecies; n++) {
        MultiFab::Add(*rhs[lev], *auxiliary[lev], n, 0, 1, rhs[lev]->nGrow());
      }
    }

    for (int lev = 0; lev <= finest_level; lev++)
      EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
  }
}


void
mfix::mfix_scalar_rhs (const bool explicit_diffusion,
                       Vector< MultiFab* > const& rhs,
                       Vector< MultiFab* > const& trac,
                       Vector< MultiFab* > const& ep_g,
                       Vector< MultiFab* > const& ro_g,
                       const Vector<Real>& mu_s_in)
{
  if (explicit_diffusion) {
    diffusion_op->ComputeLapS(rhs, trac, ro_g, ep_g, mu_s_in);

    for (int lev = 0; lev <= finest_level; lev++)
      EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
  }
}


void
mfix::mfix_species_X_rhs (const bool explicit_diffusion,
                          Vector< MultiFab* > const& rhs,
                          Vector< MultiFab* > const& X_gk,
                          Vector< MultiFab* > const& ep_g,
                          Vector< MultiFab* > const& ro_g,
                          Vector< MultiFab* > const& D_gk,
                          Vector< MultiFab* > const& ro_gk_txfr)
{
  if (explicit_diffusion) {
    diffusion_op->ComputeLapX(rhs, X_gk, ro_g, ep_g, D_gk);

    for (int lev = 0; lev <= finest_level; lev++)
      EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
  }

  if (solve_reactions) {
    for (int lev = 0; lev <= finest_level; lev++) {
      rhs[lev]->plus(*ro_gk_txfr[lev], 0, FLUID::nspecies, rhs[lev]->nGrow());
    }

    for (int lev = 0; lev <= finest_level; lev++)
      EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
  }
}
