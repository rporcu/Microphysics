#include <mfix.H>

#include <mfix_mf_helpers.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_reactions_parms.H>


void
mfix::mfix_density_rhs (Vector< MultiFab*      > const& rhs,
                        Vector< MultiFab const*> const& chem_txfr)
{
  for (int lev = 0; lev <= finest_level; lev++)
    rhs[lev]->setVal(0.);

  if (reactions.solve) {
    ChemTransfer chem_txfr_idxs(fluid.nspecies, reactions.nreactions);

    for (int lev = 0; lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Tilebox
        Box bx = mfi.tilebox();

        const int nspecies_g = fluid.nspecies;
        const int start_idx  = chem_txfr_idxs.ro_gk_txfr;

        Array4<Real      > const& rhs_arr        = rhs[lev]->array(mfi);
        Array4<Real const> const& ro_gk_txfr_arr = chem_txfr[lev]->const_array(mfi,start_idx);

        amrex::ParallelFor(bx, [nspecies_g,rhs_arr,ro_gk_txfr_arr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real rhs_loc(0);

          for (int n_g(0); n_g < nspecies_g; ++n_g)
            rhs_loc += ro_gk_txfr_arr(i,j,k,n_g);

          rhs_arr(i,j,k) = rhs_loc;
        });
      }
    }
  }

  for (int lev = 0; lev <= finest_level; lev++)
    EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
}


void
mfix::mfix_enthalpy_rhs (Vector< MultiFab*      > const& rhs,
                         Vector< MultiFab const*> const& ep_g,
                         Vector< MultiFab const*> const& ro_g,
                         Vector< MultiFab*      > const& X_gk,
                         Vector< MultiFab const*> const& T_g,
                         Vector< MultiFab const*> const& chem_txfr)
{
  for (int lev = 0; lev <= finest_level; lev++)
    rhs[lev]->setVal(0.);

  if (fluid.is_a_mixture) {
    const int nspecies_g = fluid.nspecies;

    // Temporary for computing other terms of RHS
    Vector<MultiFab*> lap_hX_gk(nlev, nullptr);

    // Allocate memory for computing fluid species contributio
    for (int lev(0); lev <= finest_level; lev++) {
      lap_hX_gk[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies,
                                    nghost_state(), MFInfo(), *ebfactory[lev]);

      lap_hX_gk[lev]->setVal(0.);
    }

    // Compute the mixed enthalpy/species term
    diffusion_op->ComputeLaphX(lap_hX_gk, X_gk, ro_g, ep_g, T_g);

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

  if (reactions.solve) {
    ChemTransfer chem_txfr_idxs(fluid.nspecies, reactions.nreactions);

    const int start_idx = chem_txfr_idxs.h_g_txfr;

    for (int lev(0); lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*ep_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Tilebox
        Box bx = mfi.tilebox();

        Array4<Real      > const& rhs_arr      = rhs[lev]->array(mfi);
        Array4<Real const> const& h_g_txfr_arr = chem_txfr[lev]->const_array(mfi,start_idx);

        amrex::ParallelFor(bx, [rhs_arr,h_g_txfr_arr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          rhs_arr(i,j,k) += h_g_txfr_arr(i,j,k);
        });
      }
    }
  }

  for (int lev = 0; lev <= finest_level; lev++)
    EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
}


void
mfix::mfix_scalar_rhs (/*Vector< MultiFab* > const& rhs,*/
                       Vector< MultiFab const* > const& /*trac*/,
                       Vector< MultiFab const* > const& /*ep_g*/,
                       Vector< MultiFab const* > const& /*ro_g*/,
                       Vector<Real> const& /*mu_s_in*/)
{
}


void
mfix::mfix_species_X_rhs (Vector< MultiFab*      > const& rhs,
                          Vector< MultiFab const*> const& chem_txfr)
{
  for (int lev = 0; lev <= finest_level; lev++)
    rhs[lev]->setVal(0.);

  if (reactions.solve) {
    const int nspecies_g = fluid.nspecies;

    ChemTransfer chem_txfr_idxs(nspecies_g, reactions.nreactions);

    const int start_idx = chem_txfr_idxs.ro_gk_txfr;

    for (int lev = 0; lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Tilebox
        Box bx = mfi.tilebox ();

        Array4<Real      > const& rhs_arr          = rhs[lev]->array(mfi);
        Array4<Real const> const& species_txfr_arr = chem_txfr[lev]->const_array(mfi,start_idx);

        amrex::ParallelFor(bx, [rhs_arr,species_txfr_arr,nspecies_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          for (int n_g(0); n_g < nspecies_g; ++n_g) {
            rhs_arr(i,j,k,n_g) = species_txfr_arr(i,j,k,n_g);
          }
        });
      }
    }
  }

  for (int lev = 0; lev <= finest_level; lev++)
    EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
}


void
mfix::mfix_momentum_rhs (Vector< MultiFab* > const& rhs,
                         Vector< MultiFab const* > const& ep_g,
                         Vector< MultiFab const* > const& chem_txfr)
{
  for (int lev = 0; lev <= finest_level; lev++)
    rhs[lev]->setVal(0.);

  if (reactions.solve) {
    ChemTransfer chem_txfr_idxs(fluid.nspecies, reactions.nreactions);

    const int start_idx = chem_txfr_idxs.vel_g_txfr;

    for (int lev = 0; lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Tilebox
        Box bx = mfi.tilebox ();

        Array4<Real      > const& rhs_arr           = rhs[lev]->array(mfi);
        Array4<Real const> const& epg_arr           = ep_g[lev]->const_array(mfi);
        Array4<Real const> const& momentum_txfr_arr = chem_txfr[lev]->const_array(mfi,start_idx);

        amrex::ParallelFor(bx, [epg_arr,rhs_arr,momentum_txfr_arr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          rhs_arr(i,j,k,0) += momentum_txfr_arr(i,j,k,0);
          rhs_arr(i,j,k,1) += momentum_txfr_arr(i,j,k,1);
          rhs_arr(i,j,k,2) += momentum_txfr_arr(i,j,k,2);
        });
      }
    }
  }

  for (int lev = 0; lev <= finest_level; lev++)
    EB_set_covered(*rhs[lev], 0, rhs[lev]->nComp(), rhs[lev]->nGrow(), 0.);
}
