#include <mfix.H>
#include <mfix_bc_parms.H>
#include <mfix_mf_helpers.H>

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
// Compute incompressibility constraint div(ep_g * u) for a closed system
//
void
mfix::mfix_closed_system_rhs (Vector< MultiFab*       > const& rhs,
                              Vector< MultiFab const* > const& lap_T,
                              Vector< MultiFab const* > const& enthalpy_rhs,
                              Vector< MultiFab const* > const& lap_X,
                              Vector< MultiFab const* > const& species_rhs,
                              Vector< MultiFab const* > const& ep_g,
                              Vector< MultiFab const* > const& ro_g,
                              Vector< MultiFab const* > const& T_g,
                              Vector< MultiFab*       > const& X_gk,
                              Vector< MultiFab const* > const& txfr,
                              Vector< MultiFab const* > const& chem_txfr,
                              Vector< MultiFab const* > const& pressure_g,
                              Vector< Real >& avgSigma,
                              Vector< Real >& avgTheta)
{
  Vector< MultiFab* > Sigma(finest_level+1);
  Vector< MultiFab* > Theta(finest_level+1);

  for (int lev = 0; lev <= finest_level; ++lev) {
    Sigma[lev] = new MultiFab(grids[lev], dmap[lev], 1, nghost_state(), MFInfo(), *ebfactory[lev]);
    Theta[lev] = new MultiFab(grids[lev], dmap[lev], 1, nghost_state(), MFInfo(), *ebfactory[lev]);

    Sigma[lev]->setVal(0.);
    Theta[lev]->setVal(0.);
  }

  mfix_open_system_rhs(Sigma, lap_T, enthalpy_rhs, lap_X, species_rhs, ro_g,
                       T_g, X_gk, txfr, chem_txfr);

  const int nspecies_g = fluid.nspecies;
  const int is_mixture = fluid.is_a_mixture;

  // IC values for MW_g
  const Real MW_g0 = fluid.MW_g0;

  Gpu::DeviceVector<Real> MW_gk0_d(nspecies_g);
  if (is_mixture)
    Gpu::copyAsync(Gpu::hostToDevice, fluid.MW_gk0.begin(), fluid.MW_gk0.end(), MW_gk0_d.begin());
  Real* p_MW_gk0 = is_mixture ? MW_gk0_d.data() : nullptr;

  auto& fluid_parms = *fluid.parameters;

  // Compute Theta
  for (int lev(0); lev <= finest_level; ++lev)
  {
    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ep_g[lev]->Factory());
    const auto& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*Theta[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      Array4< Real       > const& theta_arr  = Theta[lev]->array(mfi);
      Array4< Real const > const& ep_g_arr   = ep_g[lev]->const_array(mfi);
      Array4< Real const > const& X_gk_arr   = X_gk[lev]->const_array(mfi);
      Array4< Real const > const& T_g_arr    = T_g[lev]->const_array(mfi);
      Array4< Real const > const& pres_g_arr = pressure_g[lev]->const_array(mfi);

      auto const& flags_arr = flags.const_array(mfi);

      const Real R = fluid.R;

      ParallelFor(bx, [theta_arr,ep_g_arr,T_g_arr,X_gk_arr,pres_g_arr,R,
          flags_arr,p_MW_gk0,MW_g0,is_mixture,nspecies_g,fluid_parms]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const Real Tg_loc = T_g_arr(i,j,k);

        Real MW_g(0);
        Real cp_g(0);

        // set initial fluid molecular weight
        if (is_mixture) {
          for (int n(0); n < nspecies_g; n++) {
            MW_g += X_gk_arr(i,j,k,n) / p_MW_gk0[n];
            cp_g += X_gk_arr(i,j,k,n) * fluid_parms.calc_cp_gk(Tg_loc,n);
          }

          MW_g = 1. / MW_g;
        }
        else {
          MW_g = MW_g0;
          cp_g = fluid_parms.calc_cp_g(Tg_loc);
        }

        if (!flags_arr(i,j,k).isCovered()) {
          const Real coeff = ep_g_arr(i,j,k) / pres_g_arr(i,j,k);
          theta_arr(i,j,k) = coeff * (1 - R/(MW_g * cp_g));
        } else {
          theta_arr(i,j,k) = 0.;
        }
      });
    }
  }

  for (int lev(0); lev <= finest_level; ++lev)
  {
    Box domain(geom[lev].Domain());

    bool local = true;

    avgSigma[lev] = volWgtSumBox(lev, *Sigma[lev], 0, domain, local);
    avgTheta[lev] = volWgtSumBox(lev, *Theta[lev], 0, domain, local);

    // Compute parallel reductions
    ParallelDescriptor::ReduceRealSum(avgSigma[lev]);
    ParallelDescriptor::ReduceRealSum(avgTheta[lev]);
  }

  // Compute rhs
  for (int lev(0); lev <= finest_level; ++lev)
  {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      Array4< Real       > const& rhs_arr    = rhs[lev]->array(mfi);
      Array4< Real const > const& sigma_arr  = Sigma[lev]->const_array(mfi);
      Array4< Real const > const& theta_arr  = Theta[lev]->const_array(mfi);

      const Real avg_sigma = avgSigma[lev];
      const Real avg_theta = avgTheta[lev];

      amrex::ParallelFor(bx, [rhs_arr,sigma_arr,theta_arr,avg_sigma,avg_theta]
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


//
// Compute incompressibility constraint div(ep_g * u) for an open system
//
void
mfix::mfix_open_system_rhs (Vector< MultiFab*      > const& rhs,
                            Vector< MultiFab const*> const& lap_T,
                            Vector< MultiFab const*> const& enthalpy_rhs,
                            Vector< MultiFab const*> const& lap_X,
                            Vector< MultiFab const*> const& species_rhs,
                            Vector< MultiFab const*> const& ro_g,
                            Vector< MultiFab const*> const& T_g,
                            Vector< MultiFab*      > const& X_gk,
                            Vector< MultiFab const*> const& /*txfr*/,
                            Vector< MultiFab const*> const& /*chem_txfr*/)
{
  for (int lev(0); lev <= finest_level; ++lev)
    rhs[lev]->setVal(0.);

  const int adv_enthalpy = advect_enthalpy;
  const int is_mixture = fluid.is_a_mixture;
  const int nspecies_g = fluid.nspecies;

  Gpu::DeviceVector< Real > MW_gk_d(nspecies_g);
  if (is_mixture)
    Gpu::copyAsync(Gpu::hostToDevice, fluid.MW_gk0.begin(), fluid.MW_gk0.end(), MW_gk_d.begin());
  Real* p_MW_gk = is_mixture ? MW_gk_d.data() : nullptr;

  const Real MW_g0 = fluid.MW_g0;

  auto& fluid_parms = *fluid.parameters;

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
      auto const& lap_T_arr = adv_enthalpy ? lap_T[lev]->const_array(mfi) : empty_arr;
      auto const& h_RHS_arr = adv_enthalpy ? enthalpy_rhs[lev]->const_array(mfi) : empty_arr;
      auto const& ro_g_arr  = adv_enthalpy ? ro_g[lev]->const_array(mfi) : empty_arr;
      auto const& T_g_arr   = adv_enthalpy ? T_g[lev]->const_array(mfi) : empty_arr;
      auto const& X_gk_arr  = is_mixture ? X_gk[lev]->const_array(mfi) : empty_arr;
      auto const& lap_X_arr = is_mixture ? lap_X[lev]->const_array(mfi) : empty_arr;
      auto const& X_RHS_arr = is_mixture ? species_rhs[lev]->const_array(mfi) : empty_arr;

      auto const& flags_arr = flags.const_array(mfi);

      amrex::ParallelFor(bx, [rhs_arr,lap_T_arr,h_RHS_arr,ro_g_arr,T_g_arr,
          X_gk_arr,lap_X_arr,X_RHS_arr,flags_arr,adv_enthalpy,is_mixture,
          nspecies_g,p_MW_gk,MW_g0,fluid_parms]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real rhs_value(0);

        const Real rog_loc = ro_g_arr(i,j,k);
        const Real Tg_loc  = adv_enthalpy ? T_g_arr(i,j,k) : 0.;

        Real MW_g(0);
        Real cp_g(0);

        // set initial fluid molecular weight
        if (is_mixture) {
          for (int n(0); n < nspecies_g; n++) {
            MW_g += X_gk_arr(i,j,k,n) / p_MW_gk[n];
            cp_g += X_gk_arr(i,j,k,n) * fluid_parms.calc_cp_gk(Tg_loc,n);
          }
          MW_g = 1. / MW_g;
        }
        else {
          MW_g = MW_g0;
          cp_g = fluid_parms.calc_cp_g(Tg_loc);
        }

        if (!flags_arr(i,j,k).isCovered()) {
          if (adv_enthalpy) {
            rhs_value += (h_RHS_arr(i,j,k) + lap_T_arr(i,j,k)) / (rog_loc*cp_g*Tg_loc);
          }

          if (is_mixture) {
            for (int n(0); n < nspecies_g; ++n) {
              Real coeff = MW_g / p_MW_gk[n];

              if (adv_enthalpy) {
                const Real h_gk = fluid_parms.calc_h_gk(Tg_loc,n);
                coeff -= h_gk / (cp_g*Tg_loc);
              }

              rhs_value += (coeff / rog_loc) * (X_RHS_arr(i,j,k,n) + lap_X_arr(i,j,k,n));
            }
          }
        }

        rhs_arr(i,j,k) = rhs_value;
      });
    }
  }

  for (int lev(0); lev <= finest_level; lev++) {
    EB_set_covered(*rhs[lev], 0, 1, 0, 0.0);
  }

  for (int lev(0); lev <= finest_level; lev++) {
    rhs[lev]->FillBoundary(geom[lev].periodicity());
  }
}
