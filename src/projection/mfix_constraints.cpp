#include <mfix.H>
#include <mfix_bc_parms.H>
#include <mfix_mf_helpers.H>
#include <mfix_utils.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLNodeLaplacian.H>


//
// Compute Incompressible Fluid constraint
//
void
mfix::mfix_incompressible_fluid_rhs (Vector< MultiFab* > const& rhs)
{
  for (int lev(0); lev <= finest_level; ++lev)
    rhs[lev]->setVal(0.);

  for (int lev(0); lev <= finest_level; lev++) {
    EB_set_covered(*rhs[lev], 0, 1, 0, 0.0);
  }
}


//
// Compute Ideal Gas constraint for an open system
//
void
mfix::mfix_idealgas_opensystem_rhs (Vector< MultiFab*       > const& rhs,
                                    Vector< MultiFab const* > const& enthalpy_rhs,
                                    Vector< MultiFab const* > const& species_rhs,
                                    Vector< MultiFab const* > const& ro_g,
                                    Vector< MultiFab const* > const& T_g,
                                    Vector< MultiFab const* > const& X_gk)
{
  const int solve_enthalpy = fluid.solve_enthalpy;
  const int solve_species = fluid.solve_species;
  const int fluid_is_a_mixture = fluid.is_a_mixture;
  const int nspecies_g = fluid.nspecies;

  auto& fluid_parms = *fluid.parameters;

  if (solve_enthalpy || solve_species) {

    for (int lev(0); lev <= finest_level; ++lev)
      rhs[lev]->setVal(0.);

    for (int lev(0); lev <= finest_level; lev++) {
      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ro_g[lev]->Factory());
      const auto& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        Array4<const Real> empty_arr;

        auto const& rhs_arr   = rhs[lev]->array(mfi);
        auto const& h_RHS_arr = solve_enthalpy ? enthalpy_rhs[lev]->const_array(mfi) : empty_arr;
        auto const& ro_g_arr  = ro_g[lev]->const_array(mfi);
        auto const& T_g_arr   = T_g[lev]->const_array(mfi);
        auto const& X_gk_arr  = fluid_is_a_mixture ? X_gk[lev]->const_array(mfi) : empty_arr;
        auto const& X_RHS_arr = solve_species ? species_rhs[lev]->const_array(mfi) : empty_arr;

        auto const& flags_arr = flags.const_array(mfi);

        amrex::ParallelFor(bx, [rhs_arr,h_RHS_arr,ro_g_arr,T_g_arr,X_gk_arr,
            X_RHS_arr,flags_arr,solve_enthalpy,fluid_is_a_mixture,nspecies_g,
            fluid_parms,solve_species]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

          Real rhs_value(0);

          if (!flags_arr(i,j,k).isCovered()) {

            const Real rog_loc = ro_g_arr(i,j,k);
            const Real Tg_loc  = T_g_arr(i,j,k);

            Real MW_g_loc(0);
            Real cp_g_loc(0);

            // set initial fluid molecular weight
            if (fluid_is_a_mixture) {
              for (int n(0); n < nspecies_g; n++) {
                const Real MW_gk = fluid_parms.get_MW_gk<run_on>(n);

                const Real cp_gk = fluid_parms.calc_cp_gk<run_on>(Tg_loc,n);

                MW_g_loc += X_gk_arr(i,j,k,n) / MW_gk;
                cp_g_loc += X_gk_arr(i,j,k,n) * cp_gk;
              }
              MW_g_loc = 1. / MW_g_loc;

            } else {

              cp_g_loc = fluid_parms.calc_cp_g<run_on>(Tg_loc);
            }

            if (solve_enthalpy) {
              rhs_value += h_RHS_arr(i,j,k) / (rog_loc*cp_g_loc*Tg_loc);
            }

            if (fluid_is_a_mixture) {
              for (int n(0); n < nspecies_g; ++n) {
                const Real MW_gk = fluid_parms.get_MW_gk<run_on>(n);

                Real coeff = MW_g_loc / MW_gk;

                const Real h_gk = fluid_parms.calc_h_gk<run_on>(Tg_loc, n, cell_is_covered);

                coeff -= h_gk / (cp_g_loc*Tg_loc);

                rhs_value += (coeff / rog_loc) * X_RHS_arr(i,j,k,n);
              }
            } else if (solve_species) {

              Real coeff = 1.;

              const Real h_g_loc = fluid_parms.calc_h_g<run_on>(Tg_loc, cell_is_covered);

              coeff -= h_g_loc / (cp_g_loc*Tg_loc);

              rhs_value += (coeff / rog_loc) * X_RHS_arr(i,j,k,0);
            }
          }

          rhs_arr(i,j,k) = rhs_value;
        });
      }
    }

    for (int lev(0); lev <= finest_level; lev++) {
      EB_set_covered(*rhs[lev], 0, 1, 0, 0.0);
    }

  } else { // no solve_enthalpy and no solve_fluid_species

    this->mfix_incompressible_fluid_rhs(rhs);
  }
}


//
// Compute Ideal Gas constraint for a closed system
//
void
mfix::mfix_idealgas_closedsystem_rhs (Vector< MultiFab*       > const& rhs,
                                      Vector< MultiFab const* > const& enthalpy_rhs,
                                      Vector< MultiFab const* > const& species_rhs,
                                      Vector< MultiFab const* > const& ep_g,
                                      Vector< MultiFab const* > const& ro_g,
                                      Vector< MultiFab const* > const& T_g,
                                      Vector< MultiFab const* > const& X_gk,
                                      Vector< MultiFab const* > const& pressure_g,
                                      Vector< Real >& avgSigma,
                                      Vector< Real >& avgTheta)
{
  Vector< MultiFab* > Sigma(finest_level+1);
  Vector< MultiFab* > Theta(finest_level+1);

  for (int lev = 0; lev <= finest_level; ++lev) {
    Sigma[lev] = new MultiFab(grids[lev], dmap[lev], 1, rhs[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    Theta[lev] = new MultiFab(grids[lev], dmap[lev], 1, rhs[lev]->nGrow(), MFInfo(), *ebfactory[lev]);

    Sigma[lev]->setVal(0.);
    Theta[lev]->setVal(0.);
  }

  mfix_idealgas_opensystem_rhs(Sigma, enthalpy_rhs, species_rhs, ro_g, T_g, X_gk);

  const int solve_enthalpy = fluid.solve_enthalpy;
  const int nspecies_g = fluid.nspecies;
  const int fluid_is_a_mixture = fluid.is_a_mixture;

  auto& fluid_parms = *fluid.parameters;

  // Compute Theta
  for (int lev(0); lev <= finest_level; ++lev) {

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ep_g[lev]->Factory());
    const auto& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*Theta[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

      Array4< Real const > dummy_arr;

      Array4< Real       > const& theta_arr  = Theta[lev]->array(mfi);
      Array4< Real const > const& ep_g_arr   = ep_g[lev]->const_array(mfi);
      Array4< Real const > const& ro_g_arr   = ro_g[lev]->const_array(mfi);
      Array4< Real const > const& X_gk_arr   = fluid_is_a_mixture ? X_gk[lev]->const_array(mfi) : dummy_arr;
      Array4< Real const > const& T_g_arr    = T_g[lev]->const_array(mfi);
      Array4< Real const > const& pres_g_arr = pressure_g[lev]->const_array(mfi);

      auto const& flags_arr = flags.const_array(mfi);

      ParallelFor(bx, [theta_arr,ep_g_arr,T_g_arr,X_gk_arr,pres_g_arr,ro_g_arr,
          flags_arr,fluid_is_a_mixture,nspecies_g,fluid_parms,solve_enthalpy]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (!flags_arr(i,j,k).isCovered()) {

          const Real Tg_loc = T_g_arr(i,j,k);

          Real cp_g_loc(0);

          // set initial fluid molecular weight
          if (fluid_is_a_mixture) {
            for (int n(0); n < nspecies_g; n++) {

              const Real cp_gk = fluid_parms.calc_cp_gk<run_on>(Tg_loc,n);

              cp_g_loc += X_gk_arr(i,j,k,n) * cp_gk;
            }

          } else {

            cp_g_loc = fluid_parms.calc_cp_g<run_on>(Tg_loc);
          }

          theta_arr(i,j,k) = ep_g_arr(i,j,k) * ((1./pres_g_arr(i,j,k)) -
              (1./(cp_g_loc*ro_g_arr(i,j,k)*Tg_loc)));
        }
      });
    }
  }

  for (int lev(0); lev <= finest_level; ++lev) {

    const amrex::MultiFab& volfrac = ebfactory[lev]->getVolFrac();

    avgSigma[lev] = Utils::volEpsWgtAvg(*Sigma[lev], 0, *ep_g[lev], volfrac);
    avgTheta[lev] = Utils::volEpsWgtAvg(*Theta[lev], 0, *ep_g[lev], volfrac);
  }

  // Compute rhs
  for (int lev(0); lev <= finest_level; ++lev) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

      Array4< Real       > const& rhs_arr    = rhs[lev]->array(mfi);
      Array4< Real const > const& sigma_arr  = Sigma[lev]->const_array(mfi);
      Array4< Real const > const& theta_arr  = Theta[lev]->const_array(mfi);

      const Real avg_sigma = avgSigma[lev];
      const Real avg_theta = avgTheta[lev];

      amrex::ParallelFor(bx, [rhs_arr,sigma_arr,theta_arr,avg_sigma,avg_theta,
          solve_enthalpy]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const Real d_sigma = sigma_arr(i,j,k) - avg_sigma;
        const Real d_theta = theta_arr(i,j,k) - avg_theta;

        rhs_arr(i,j,k) = d_sigma - (avg_sigma/avg_theta)*d_theta;
      });
    }
  }

  for (int lev = 0; lev <= finest_level; ++lev) {
    delete Sigma[lev];
    delete Theta[lev];
  }
}
