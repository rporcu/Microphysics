#include <mfix.H>
#include <mfix_bc_parms.H>
#include <mfix_fluid_parms.H>
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
// Compute incompressibility constraint div(ep_g * u) for an open system
//
void
mfix::mfix_open_system_rhs (Vector< MultiFab*      > const& rhs,
                            Vector< MultiFab*      > const& lap_T_star,
                            Vector< MultiFab*      > const& lap_X_star,
                            Vector< MultiFab const*> const& ep_g,
                            Vector< MultiFab const*> const& ro_g,
                            Vector< MultiFab const*> const& MW_g,
                            Vector< MultiFab const*> const& T_g,
                            Vector< MultiFab const*> const& cp_g,
                            Vector< MultiFab const*> const& X_gk,
                            Vector< MultiFab const*> const& D_gk,
                            Vector< MultiFab const*> const& h_gk,
                            Vector< MultiFab const*> const& txfr,
                            Vector< MultiFab const*> const& ro_gk_txfr)
{
  Vector< MultiFab* > S_h(nlev, nullptr);
  Vector< MultiFab* > S_sk(nlev, nullptr);

  const int adv_enthalpy = advect_enthalpy;
  const int fluid_is_mixture = FLUID::is_a_mixture;

  if (adv_enthalpy)
  {
    for (int lev(0); lev <= finest_level; lev++) {
      S_h[lev] = MFHelpers::createFrom(*rhs[lev], 0., 1).release();
    }

    // Compute S_h, aka enthalpy RHS
    mfix_enthalpy_rhs(S_h, ep_g, ro_g, X_gk, D_gk, h_gk);

    for (int lev(0); lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        Array4<      Real> const& S_h_arr   = S_h[lev]->array(mfi);
        Array4<const Real> const& lap_T_arr = lap_T_star[lev]->const_array(mfi);
        Array4<const Real> const& txfr_arr  = txfr[lev]->const_array(mfi);
        Array4<const Real> const& T_g_arr   = T_g[lev]->const_array(mfi);

        amrex::ParallelFor(bx, [S_h_arr,txfr_arr,T_g_arr,lap_T_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          S_h_arr(i,j,k) += lap_T_arr(i,j,k)
            + txfr_arr(i,j,k,Transfer::gammaTp)
            - txfr_arr(i,j,k,Transfer::gamma)*T_g_arr(i,j,k);
        });
      }
    }

    // Divide S_h by ( ro_g c_p T_g )
    for (int lev(0); lev <= finest_level; lev++) {
      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(S_h[lev]->Factory());
      const auto& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        Array4< Real > const& S_h_arr  = S_h[lev]->array(mfi);
        Array4< const Real > const& ro_g_arr = ro_g[lev]->const_array(mfi);
        Array4< const Real > const& T_g_arr  = T_g[lev]->const_array(mfi);
        Array4< const Real > const& cp_g_arr = cp_g[lev]->const_array(mfi);
        auto const& flags_arr = flags.const_array(mfi);

        amrex::ParallelFor(bx, [flags_arr,S_h_arr,ro_g_arr,cp_g_arr,T_g_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (not flags_arr(i,j,k).isCovered())
            S_h_arr(i,j,k) /= (ro_g_arr(i,j,k)*cp_g_arr(i,j,k)*T_g_arr(i,j,k));
          else
            S_h_arr(i,j,k) = 0.0;
        });
      }
    }
  }

  const int nspecies_g = FLUID::nspecies;

  if (fluid_is_mixture)
  {
    for (int lev(0); lev <= finest_level; lev++) {
      S_sk[lev] = MFHelpers::createFrom(*rhs[lev], 0., 1, nspecies_g).release();
    }

    Gpu::DeviceVector< Real > MW_gk_d(nspecies_g);
    Gpu::copyAsync(Gpu::hostToDevice, FLUID::MW_gk0.begin(), FLUID::MW_gk0.end(),
                   MW_gk_d.begin());
    Gpu::synchronize();
    Real* p_MW_gk = MW_gk_d.data();

    // compute S_sk
    mfix_species_X_rhs(S_sk, ro_gk_txfr);

    for (int lev(0); lev <= finest_level; lev++)
    {
      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(rhs[lev]->Factory());
      const auto& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.tilebox();

        Array4<       Real > const& S_sk_arr  = S_sk[lev]->array(mfi);
        Array4< const Real > const& lap_X_arr = lap_X_star[lev]->const_array(mfi);
        Array4< const Real > const& ro_g_arr  = ro_g[lev]->const_array(mfi);
        Array4< const Real > const& MW_g_arr  = MW_g[lev]->const_array(mfi);

        Array4< const Real > const& T_g_arr  = adv_enthalpy ?
          T_g[lev]->const_array(mfi) : Array4<const Real>();
        Array4< const Real > const& cp_g_arr = adv_enthalpy ?
          cp_g[lev]->const_array(mfi) : Array4<const Real>();
        Array4< const Real > const& h_gk_arr = adv_enthalpy ?
          h_gk[lev]->const_array(mfi) : Array4<const Real>();
        auto const& flags_arr = flags.const_array(mfi);

        amrex::ParallelFor(bx, nspecies_g, [flags_arr,S_sk_arr,MW_g_arr,p_MW_gk,
            h_gk_arr,cp_g_arr,T_g_arr,ro_g_arr,adv_enthalpy,lap_X_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          if (not flags_arr(i,j,k).isCovered()) {
            S_sk_arr(i,j,k,n) += lap_X_arr(i,j,k,n);

            amrex::Real coeff = MW_g_arr(i,j,k)/p_MW_gk[n];

            if (adv_enthalpy)
              coeff -= h_gk_arr(i,j,k,n)/(cp_g_arr(i,j,k)*T_g_arr(i,j,k));

            S_sk_arr(i,j,k,n) *= coeff / ro_g_arr(i,j,k);
          }
          else {
            S_sk_arr(i,j,k,n) = 0.;
          }
        });
      }
    }
  }

  for (int lev(0); lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      Array4<Real> const& rhs_arr  = rhs[lev]->array(mfi);
      Array4<const Real> const& S_h_arr = adv_enthalpy ?
        S_h[lev]->const_array(mfi) : Array4<const Real>();
      Array4<const Real> const& S_sk_arr = fluid_is_mixture ?
        S_sk[lev]->const_array(mfi) : Array4<const Real>();

      amrex::ParallelFor(bx, [S_h_arr,S_sk_arr,nspecies_g,rhs_arr,adv_enthalpy,
          fluid_is_mixture]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real rhs_sum(0);

        if (adv_enthalpy) {
          rhs_sum += S_h_arr(i,j,k);
        }

        if (fluid_is_mixture) {
          for (int n(0); n < nspecies_g; n++)
            rhs_sum += S_sk_arr(i,j,k,n);
        }

        rhs_arr(i,j,k) = rhs_sum;
      });
    }
  }

  for (int lev(0); lev <= finest_level; lev++) {
    EB_set_covered(*rhs[lev], 0, 1, 0, 0.0);
  }

  for (int lev(0); lev <= finest_level; lev++) {
    rhs[lev]->FillBoundary(geom[lev].periodicity());
  }

  for (int lev(0); lev <= finest_level; lev++)
  {
    if (adv_enthalpy) {
      delete S_h[lev];
    }

    if (fluid_is_mixture) {
      delete S_sk[lev];
    }
  }
}


//
// Compute div(ep_g * u)
//
void
mfix::mfix_compute_diveu (Real time)
{
  // Note that the solver imposes the boundary conditions with the right scalings so we don't
  //      fill any ghost cells here.
  Vector< MultiFab* > epu(nlev, nullptr);

  for (int lev = 0; lev < nlev; lev++)
    {
      MultiFab& vel_g = *(m_leveldata[lev]->vel_g);
      // We only need one ghost cell here -- so no need to make it bigger
      epu[lev] = new MultiFab(vel_g.boxArray(), vel_g.DistributionMap(),
                              vel_g.nComp(), 1 , MFInfo(), *ebfactory[lev]);

      epu[lev]->setVal(1.e200);

      Box domain(geom[lev].Domain());

      MultiFab::Copy(*epu[lev], vel_g, 0, 0, 3, epu[lev]->nGrow());

      for (int n = 0; n < 3; n++)
        MultiFab::Multiply(*epu[lev], *(m_leveldata[lev]->ep_g), 0, n, 1, epu[lev]->nGrow());

      epu[lev]->FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      // Extrapolate Dirichlet values to ghost cells -- but do it differently in that
      //  no-slip walls are treated exactly like slip walls --
      // Note that this routine is essential to impose the correct inflow bc's on
      //  the product ep_g * vel_g
      for (MFIter mfi(*epu[lev], false); mfi.isValid(); ++mfi)
      {
        set_vec_bcs(lev, (*epu[lev])[mfi], domain);
      }

      epu[lev]->FillBoundary(geom[lev].periodicity());

      // We set these to zero because if the values in the covered cells are undefined,
      //   even though they are multiplied by zero in the divu computation, we can still get NaNs
      EB_set_covered(*epu[lev], 0, epu[lev]->nComp(), 1, 0.0);
    }

  // Define the operator in order to compute the multi-level divergence
  //
  //        (del dot b sigma grad)) phi
  //
  LPInfo info;
  MLNodeLaplacian matrix(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));

  // Set domain BCs for Poisson's solver
  // The domain BCs refer to level 0 only

  matrix.setDomainBC(BC::ppe_lobc, BC::ppe_hibc);

  matrix.compDivergence(get_diveu(), epu);

  for(int lev(0); lev < nlev; lev++)
    delete epu[lev];

  // Restore velocities to carry Dirichlet values on faces
  int extrap_dir_bcs = 0;
  mfix_set_velocity_bcs(time, get_vel_g(), extrap_dir_bcs);
}
