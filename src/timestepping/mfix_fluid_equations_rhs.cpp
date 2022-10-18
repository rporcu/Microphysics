#include <mfix.H>

#include <mfix_mf_helpers.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_fluid.H>
#include <mfix_species.H>
#include <mfix_reactions.H>


void
mfix::mfix_density_rhs (Vector< MultiFab*      > const& rhs,
                        Vector< MultiFab const*> const& txfr)
{
  for (int lev = 0; lev <= finest_level; lev++)
    rhs[lev]->setVal(0.);

  if (reactions.solve()) {
    InterphaseTxfrIndexes txfr_idxs(fluid.nspecies(), reactions.nreactions());

    for (int lev = 0; lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Tilebox
        Box bx = mfi.tilebox();

        const int nspecies_g = fluid.nspecies();
        const int start_idx  = txfr_idxs.chem_ro_gk;

        Array4<Real      > const& rhs_arr        = rhs[lev]->array(mfi);
        Array4<Real const> const& ro_gk_txfr_arr = txfr[lev]->const_array(mfi,start_idx);

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
                         Vector< MultiFab const*> const& /*ro_g*/,
                         Vector< MultiFab*      > const& /*X_gk*/,
                         Vector< MultiFab const*> const& /*T_g*/,
                         Vector< MultiFab const*> const& txfr)
{
  for (int lev = 0; lev <= finest_level; lev++)
    rhs[lev]->setVal(0.);

  if (reactions.solve()) {
    InterphaseTxfrIndexes txfr_idxs(fluid.nspecies(), reactions.nreactions());

    const int start_idx = txfr_idxs.chem_h;

    for (int lev(0); lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*ep_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Tilebox
        Box bx = mfi.tilebox();

        Array4<Real      > const& rhs_arr      = rhs[lev]->array(mfi);
        Array4<Real const> const& h_g_txfr_arr = txfr[lev]->const_array(mfi,start_idx);

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
mfix::mfix_species_X_rhs (Vector< MultiFab*      > const& rhs,
                          Vector< MultiFab const*> const& txfr)
{
  for (int lev = 0; lev <= finest_level; lev++)
    rhs[lev]->setVal(0.);

  if (reactions.solve()) {
    const int nspecies_g = fluid.nspecies();

    InterphaseTxfrIndexes txfr_idxs(nspecies_g, reactions.nreactions());

    const int ro_gk_txfr_idx = txfr_idxs.chem_ro_gk;

    for (int lev = 0; lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Tilebox
        Box bx = mfi.tilebox ();

        Array4<Real      > const& rhs_arr        = rhs[lev]->array(mfi);
        Array4<Real const> const& ro_gk_txfr_arr = txfr[lev]->const_array(mfi,ro_gk_txfr_idx);

        amrex::ParallelFor(bx, [rhs_arr,ro_gk_txfr_arr,nspecies_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          for (int n_g(0); n_g < nspecies_g; ++n_g) {
            rhs_arr(i,j,k,n_g) = ro_gk_txfr_arr(i,j,k,n_g);
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
                         Vector< MultiFab const* > const& txfr)
{
  for (int lev = 0; lev <= finest_level; lev++)
    rhs[lev]->setVal(0.);

  if (reactions.solve()) {
    InterphaseTxfrIndexes txfr_idxs(fluid.nspecies(), reactions.nreactions());

    const int start_idx = txfr_idxs.chem_vel;

    for (int lev = 0; lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Tilebox
        Box bx = mfi.tilebox ();

        Array4<Real      > const& rhs_arr           = rhs[lev]->array(mfi);
        Array4<Real const> const& epg_arr           = ep_g[lev]->const_array(mfi);
        Array4<Real const> const& momentum_txfr_arr = txfr[lev]->const_array(mfi,start_idx);

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
