#include <mfix.H>

#include <mfix_mf_helpers.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>


void
mfix::mfix_density_rhs (Vector< MultiFab*      > const& rhs,
                        Vector< MultiFab const*> const& ro_gk_txfr)
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
mfix::mfix_enthalpy_rhs (Vector< MultiFab*      > const& rhs,
                         Vector< MultiFab const*> const& ep_g,
                         Vector< MultiFab const*> const& ro_g,
                         Vector< MultiFab*      > const& X_gk,
                         Vector< MultiFab const*> const& D_gk,
                         Vector< MultiFab const*> const& h_gk)
{
  for (int lev = 0; lev <= finest_level; lev++)
    rhs[lev]->setVal(0.);

  if (FLUID::is_a_mixture)
  {
    const int nspecies_g = FLUID::nspecies;

    // Temporary for computing other terms of RHS
    Vector<MultiFab*> lap_hX_gk(nlev, nullptr);

    // Allocate memory for computing fluid species contributio
    for (int lev(0); lev <= finest_level; lev++) {
      lap_hX_gk[lev] = new MultiFab(grids[lev], dmap[lev], FLUID::nspecies,
                                    nghost_state(), MFInfo(), *ebfactory[lev]);
    }

    // Compute the mixed enthalpy/species term
    diffusion_op->ComputeLaphX(lap_hX_gk, X_gk, ro_g, ep_g, D_gk, h_gk);

    for (int lev(0); lev <= finest_level; lev++) {
      // Add the contribution due to the nth specie
      for (int n(0); n < nspecies_g; n++) {
        MultiFab::Add(*rhs[lev], *lap_hX_gk[lev], n, 0, 1, rhs[lev]->nGrow());
      }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
      delete lap_hX_gk[lev];
    }
  }
}


void
mfix::mfix_scalar_rhs (const bool /*explicit_diffusion*/,
                       Vector< MultiFab* > const& /*lap_trac*/,
                       Vector< MultiFab* > const& /*trac*/,
                       Vector< MultiFab* > const& /*ep_g*/,
                       Vector< MultiFab* > const& /*ro_g*/,
                       const Vector<Real>& /*mu_s_in*/)
{
}


void
mfix::mfix_species_X_rhs (Vector< MultiFab*      > const& rhs,
                          Vector< MultiFab const*> const& ro_gk_txfr)
{

  if (solve_reactions) {

    for (int lev = 0; lev <= finest_level; lev++) {
      rhs[lev]->plus(*ro_gk_txfr[lev], 0, FLUID::nspecies, rhs[lev]->nGrow());
      EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
    }

  } else {

    for (int lev = 0; lev <= finest_level; lev++)
      rhs[lev]->setVal(0.);

  }
}
