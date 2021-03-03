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
                         Vector< MultiFab const*> const& /*X_gk*/,
                         Vector< MultiFab const*> const& D_gk,
                         Vector< MultiFab const*> const& h_gk)
{
  for (int lev = 0; lev <= finest_level; lev++)
    rhs[lev]->setVal(0.);

  if (FLUID::is_a_mixture)
  {

    // Temporary for computing other terms of RHS
    Vector< MultiFab* > X_gk_tmp(nlev, nullptr);
    Vector< MultiFab* > h_gk_D_gk(nlev, nullptr);
    Vector< MultiFab* > auxiliary(nlev, nullptr);
#if 0
    // Allocate memory for computing fluid species contribution
    for (int lev(0); lev <= finest_level; lev++) {
      h_gk_D_gk[lev] = MFHelpers::createFrom(*D_gk[lev]).release();
      auxiliary[lev] = MFHelpers::createFrom(*X_gk[lev], 0.).release();
    }
#endif

    // Allocate memory for computing fluid species contribution
    for (int lev(0); lev <= finest_level; lev++) {

      X_gk_tmp[lev] = new MultiFab(grids[lev], dmap[lev], FLUID::nspecies,
                                   nghost_state(), MFInfo(), *ebfactory[lev]);

      h_gk_D_gk[lev] = new MultiFab(grids[lev], dmap[lev], FLUID::nspecies,
                                    nghost_state(), MFInfo(), *ebfactory[lev]);

      auxiliary[lev] = new MultiFab(grids[lev], dmap[lev], FLUID::nspecies,
                                    nghost_state(), MFInfo(), *ebfactory[lev]);

    }

    // Transform h_gk_D_gk into h_gk * D_gk
    for (int lev(0); lev <= finest_level; lev++) {

      MultiFab::Copy(*X_gk_tmp[lev], *D_gk[lev], 0, 0, FLUID::nspecies,
                     X_gk_tmp[lev]->nGrow());

      MultiFab::Copy(*h_gk_D_gk[lev], *D_gk[lev], 0, 0, FLUID::nspecies,
                     h_gk_D_gk[lev]->nGrow());

      MultiFab::Multiply(*h_gk_D_gk[lev], *h_gk[lev], 0, 0, FLUID::nspecies,
                         h_gk_D_gk[lev]->nGrow());
    }

    // Compute the mixed enthalpy/species term
    diffusion_op->ComputeLapX(auxiliary, X_gk_tmp, ro_g, ep_g, GetVecOfConstPtrs(h_gk_D_gk));

    for (int lev(0); lev <= finest_level; lev++) {
      // Add the contribution due to the nth specie
      for (int n(0); n < FLUID::nspecies; n++) {
        MultiFab::Add(*rhs[lev], *auxiliary[lev], n, 0, 1, rhs[lev]->nGrow());
      }
    }

    for (int lev = 0; lev <= finest_level; lev++)
      EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);

    for (int lev = 0; lev <= finest_level; lev++) {
      delete X_gk_tmp[lev];
      delete h_gk_D_gk[lev];
      delete auxiliary[lev];
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
