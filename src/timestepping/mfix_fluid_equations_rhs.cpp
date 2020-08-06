#include <mfix.H>

#include <mfix_mf_helpers.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>


void
mfix::mfix_momentum_rhs (Vector< MultiFab* > const& rhs,
                         Vector< MultiFab* > const& vel_g)
{
  diffusion_op->ComputeDivTau(rhs, vel_g, get_ro_g(), get_ep_g(), get_mu_g());

  for (int lev = 0; lev <= finest_level; lev++)
    EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
}

void
mfix::mfix_enthalpy_rhs (Vector< MultiFab* > const& rhs,
                         Vector< MultiFab* > const& T_g,
                         Vector< MultiFab* > const& X_gk)
{
  diffusion_op->ComputeLapT(rhs, T_g, get_ep_g(), get_k_g(), get_T_g_on_eb(),
      get_k_g_on_eb());

  if (FLUID::is_a_mixture)
  {
    Vector< MultiFab* > D_gk = get_D_gk();
    Vector< MultiFab* > h_gk = get_h_gk();

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
      MultiFab::Multiply(*h_gk_D_gk[lev], *h_gk[lev], 0, 0, FLUID::nspecies_g,
          h_gk_D_gk[lev]->nGrow());
    }

    // Compute the mixed enthalpy/species term
    diffusion_op->ComputeLapX(auxiliary, X_gk, get_ro_g(), get_ep_g(), h_gk_D_gk);

    for (int lev(0); lev <= finest_level; lev++) {
      // Add the contribution due to the nth specie
      for (int n(0); n < FLUID::nspecies_g; n++) {
        MultiFab::Add(*rhs[lev], *auxiliary[lev], n, 0, 1, rhs[lev]->nGrow());
      }
    }
  }

  // TODO: Add term due to interphase enthalpy transfer?  right now it happens
  // after predictor and corrector, implicitly. We might consider to move it
  // here.

  for (int lev = 0; lev <= finest_level; lev++)
    EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
}

void
mfix::mfix_scalar_rhs (Vector< MultiFab* > const& rhs,
                       Vector< MultiFab* > const& sca)
{
  diffusion_op->ComputeLapS(rhs, get_trac_old(), get_ro_g(), get_ep_g(),
      mu_s);

  for (int lev = 0; lev <= finest_level; lev++)
    EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
}

void
mfix::mfix_species_X_rhs (Vector< MultiFab* > const& rhs,
                          Vector< MultiFab* > const& X_gk)
{
  diffusion_op->ComputeLapX(rhs, X_gk, get_ro_g(), get_ep_g(), get_D_gk());

  for (int lev = 0; lev <= finest_level; lev++)
    EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
}
