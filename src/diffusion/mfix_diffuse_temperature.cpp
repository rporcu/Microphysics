#include <mfix_diffusion_op.H>

#include <mfix.H>
#include <mfix_eb_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_solvers.H>

#include <AMReX_EB_utils.H>

using namespace amrex;
using namespace Solvers;

//
// Implicit solve for scalar diffusion
//
void DiffusionOp::diffuse_temperature (const Vector< MultiFab* >& T_g,
                                       const Vector< MultiFab* >& ep_g,
                                       const Vector< MultiFab* >& ro_g,
                                       const Vector< MultiFab* >& h_g,
                                       const Vector< MultiFab* >& X_gk,
                                       const Vector< MultiFab* >& T_g_on_eb,
                                       Real dt,
                                       const Real abstol,
                                       const Real reltol,
                                       const int maxiter)
{
  BL_PROFILE("DiffusionOp::diffuse_temperature");

  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

  int finest_level = amrcore->finestLevel();

  // Update the coefficients of the matrix going into the solve based on the
  // current state of the simulation. Recall that the relevant matrix is
  //
  //      alpha a - beta div ( b grad )   <--->   rho - dt div ( k_g grad )
  //
  // So the constants and variable coefficients are:
  //
  //      alpha: 0
  //      beta: dt
  //      A: ro_g ep_g
  //      B: ep_g k_g

  if(verbose > 0)
    amrex::Print() << "Diffusing temperature ..." << std::endl;

  auto& fluid_parms = *fluid.parameters;
  const int fluid_is_a_mixture = fluid.is_a_mixture;
  const int nspecies_g = fluid.nspecies;

  DampedNewton::ResidueMF R = [&] (const Vector< MultiFab* >& residue,
                                   const Vector< MultiFab* >& Tg_arg) -> void
  {
    for(int lev = 0; lev <= finest_level; lev++) {
      // Zero these out just to have a clean start
      residue[lev]->setVal(0.0);
    }

    Vector< MultiFab* > residue_aux(finest_level+1, nullptr);

    for(int lev = 0; lev <= finest_level; lev++) {
      residue_aux[lev] = new MultiFab(residue[lev]->boxArray(), residue[lev]->DistributionMap(),
                                      residue[lev]->nComp(), residue[lev]->nGrow(),
                                      MFInfo(), residue[lev]->Factory());
      residue_aux[lev]->setVal(0.0);
    }

    // Set alpha and beta
    temperature_matrix->setScalars(0, dt);

    Vector<BCRec> bcs_dummy; // This is just to satisfy the call to EBterp...
    bcs_dummy.resize(3);

    for(int lev = 0; lev <= finest_level; lev++)
    {
      MultiFab ep_k_g(ep_g[lev]->boxArray(), ep_g[lev]->DistributionMap(),
                      ep_g[lev]->nComp(), 1, MFInfo(), ep_g[lev]->Factory());

      // Initialize to 0
      ep_k_g.setVal(0.);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.growntilebox(IntVect(1,1,1));

        if (bx.ok())
        {
          Array4<Real      > const& ep_k_g_array = ep_k_g.array(mfi);
          Array4<Real const> const& ep_g_array   = ep_g[lev]->const_array(mfi);
          Array4<Real const> const& T_g_array    = Tg_arg[lev]->const_array(mfi);

          amrex::ParallelFor(bx, [ep_g_array,T_g_array,ep_k_g_array,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const Real Tg = T_g_array(i,j,k);
            ep_k_g_array(i,j,k) = ep_g_array(i,j,k)*fluid_parms.calc_k_g(Tg);
          });
        }
      }

//      ep_k_g.FillBoundary(geom[lev].periodicity());

      EB_interp_CellCentroid_to_FaceCentroid (ep_k_g, GetArrOfPtrs(b[lev]), 0, 0, 1, geom[lev], bcs_dummy);

//      b[lev][0]->FillBoundary(geom[lev].periodicity());
//      b[lev][1]->FillBoundary(geom[lev].periodicity());
//      b[lev][2]->FillBoundary(geom[lev].periodicity());

      if (EB::fix_temperature) {
        // The following is a WIP in AMReX
        //temperature_matrix->setPhiOnCentroid();

        MultiFab k_g_on_eb(T_g_on_eb[lev]->boxArray(), T_g_on_eb[lev]->DistributionMap(),
                           T_g_on_eb[lev]->nComp(), T_g_on_eb[lev]->nGrow(), MFInfo(),
                           T_g_on_eb[lev]->Factory());

        k_g_on_eb.setVal(0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*T_g_on_eb[lev]); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.growntilebox(IntVect(1,1,1));

          if (bx.ok()) {
            Array4<Real      > const& k_g_on_eb_array = k_g_on_eb.array(mfi);
            Array4<Real const> const& T_g_on_eb_array = T_g_on_eb[lev]->const_array(mfi);

            amrex::ParallelFor(bx, [k_g_on_eb_array,T_g_on_eb_array,fluid_parms]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              if (T_g_on_eb_array(i,j,k) > 0)
                k_g_on_eb_array(i,j,k) = fluid_parms.calc_k_g(T_g_on_eb_array(i,j,k));
            });
          }
        }

        k_g_on_eb.FillBoundary(geom[lev].periodicity());

        temperature_matrix->setEBDirichlet(lev, *T_g_on_eb[lev], k_g_on_eb);
      }

      temperature_matrix->setBCoeffs(lev, GetArrOfConstPtrs(b[lev]),
                                     MLMG::Location::FaceCentroid);

      temperature_matrix->setLevelBC(lev, GetVecOfConstPtrs(Tg_arg)[lev]);
    }

    MLMG solver(*temperature_matrix);

    solver.apply(residue_aux, Tg_arg);

    for(int lev = 0; lev <= finest_level; lev++) {
      MultiFab::Copy(*residue[lev], *residue_aux[lev], 0, 0, 1, 1);
      EB_set_covered(*residue[lev], 0, residue[lev]->nComp(), residue[lev]->nGrow(), 0.);
    }

    for(int lev = 0; lev <= finest_level; lev++) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.growntilebox(IntVect(1,1,1));

        const EBFArrayBox& epg_fab = static_cast<EBFArrayBox const&>((*ep_g[lev])[mfi]);
        const EBCellFlagFab& flags = epg_fab.getEBCellFlagFab();

        if (bx.ok()) {

          Array4<Real const> dummy_arr;

          Array4<Real      > const& residue_array = residue[lev]->array(mfi);
          Array4<Real const> const& ep_g_array    = ep_g[lev]->const_array(mfi);
          Array4<Real const> const& ro_g_array    = ro_g[lev]->const_array(mfi);
          Array4<Real const> const& T_g_array     = Tg_arg[lev]->const_array(mfi);
          Array4<Real const> const& h_g_array     = h_g[lev]->const_array(mfi);
          Array4<Real const> const& X_gk_array    = fluid_is_a_mixture ? X_gk[lev]->const_array(mfi) : dummy_arr;

          auto const& flags_arr = flags.const_array();

          amrex::ParallelFor(bx, [residue_array,ep_g_array,T_g_array,ro_g_array,
              h_g_array,X_gk_array,fluid_parms,fluid_is_a_mixture,nspecies_g,
              flags_arr,run_on_device]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

            const Real ep_ro_g = ep_g_array(i,j,k)*ro_g_array(i,j,k);

            const Real Tg = T_g_array(i,j,k);
            Real hg(0.);

            if (!fluid_is_a_mixture) {
              hg = run_on_device ?
                fluid_parms.calc_h_g<RunOn::Device>(Tg, cell_is_covered) :
                fluid_parms.calc_h_g<RunOn::Host>(Tg, cell_is_covered);
            } else {
              for (int n_g(0); n_g < nspecies_g; ++n_g) {
                const Real Xgk = X_gk_array(i,j,k,n_g);

                hg += run_on_device ?
                  Xgk*fluid_parms.calc_h_gk<RunOn::Device>(Tg, n_g, cell_is_covered) :
                  Xgk*fluid_parms.calc_h_gk<RunOn::Host>(Tg, n_g, cell_is_covered);
              }
            }

            residue_array(i,j,k) += ep_ro_g * hg;

            // subtract the RHS (ep_g*ro_g*tilde_h_g)
            residue_array(i,j,k) -= ep_ro_g * h_g_array(i,j,k);
          });
        }
      }

      residue[lev]->FillBoundary(geom[lev].periodicity());
      EB_set_covered(*residue[lev], 0, residue[lev]->nComp(), residue[lev]->nGrow(), 0.);

      delete residue_aux[lev];
    }

    return;
  };

  DampedNewton::GradientMF partial_R = [&] (const Vector< MultiFab* >& gradient,
                                            const Vector< MultiFab* >& Tg_arg) -> void
  {
    for(int lev = 0; lev <= finest_level; lev++) {
      // Zero these out just to have a clean start
      gradient[lev]->setVal(0.0);
    }

    Vector< MultiFab* > gradient_aux(finest_level+1, nullptr);

    for(int lev = 0; lev <= finest_level; lev++) {
      gradient_aux[lev] = new MultiFab(gradient[lev]->boxArray(), gradient[lev]->DistributionMap(),
                                       gradient[lev]->nComp(), gradient[lev]->nGrow(),
                                       MFInfo(), gradient[lev]->Factory());
      gradient_aux[lev]->setVal(0.0);
    }

    // Set alpha and beta
    temperature_matrix->setScalars(0, dt);

    Vector<BCRec> bcs_dummy; // This is just to satisfy the call to EBterp...
    bcs_dummy.resize(3);

    for(int lev = 0; lev <= finest_level; lev++) {

      MultiFab ep_dk_g(ep_g[lev]->boxArray(), ep_g[lev]->DistributionMap(),
                       ep_g[lev]->nComp(), 1, MFInfo(), ep_g[lev]->Factory());

      // Initialize to 0
      ep_dk_g.setVal(0.);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.growntilebox(IntVect(1,1,1));

        if (bx.ok())
        {
          Array4<Real      > const& ep_dk_g_array = ep_dk_g.array(mfi);
          Array4<Real const> const& ep_g_array    = ep_g[lev]->const_array(mfi);
          Array4<Real const> const& T_g_array     = Tg_arg[lev]->const_array(mfi);

          amrex::ParallelFor(bx, [ep_g_array,T_g_array,ep_dk_g_array,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const Real Tg = T_g_array(i,j,k);
            ep_dk_g_array(i,j,k) = ep_g_array(i,j,k)*fluid_parms.calc_partial_k_g(Tg);
          });
        }
      }

//      ep_dk_g.FillBoundary(geom[lev].periodicity());

      EB_interp_CellCentroid_to_FaceCentroid (ep_dk_g, GetArrOfPtrs(b[lev]), 0, 0, 1, geom[lev], bcs_dummy);

//      b[lev][0]->FillBoundary(geom[lev].periodicity());
//      b[lev][1]->FillBoundary(geom[lev].periodicity());
//      b[lev][2]->FillBoundary(geom[lev].periodicity());

      if (EB::fix_temperature) {
        // The following is a WIP in AMReX
        //temperature_matrix->setPhiOnCentroid();

        MultiFab dk_g_on_eb(T_g_on_eb[lev]->boxArray(), T_g_on_eb[lev]->DistributionMap(),
                            T_g_on_eb[lev]->nComp(), T_g_on_eb[lev]->nGrow(), MFInfo(),
                            T_g_on_eb[lev]->Factory());

        dk_g_on_eb.setVal(0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*T_g_on_eb[lev]); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.growntilebox(IntVect(1,1,1));

          if (bx.ok()) {
            Array4<Real      > const& dk_g_on_eb_array = dk_g_on_eb.array(mfi);
            Array4<Real const> const& T_g_on_eb_array  = T_g_on_eb[lev]->const_array(mfi);

            amrex::ParallelFor(bx, [dk_g_on_eb_array,T_g_on_eb_array,fluid_parms]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              if (T_g_on_eb_array(i,j,k) > 0)
                dk_g_on_eb_array(i,j,k) = fluid_parms.calc_partial_k_g(T_g_on_eb_array(i,j,k));
            });
          }
        }

        dk_g_on_eb.FillBoundary(geom[lev].periodicity());

        temperature_matrix->setEBDirichlet(lev, *T_g_on_eb[lev], dk_g_on_eb);
      }

      temperature_matrix->setBCoeffs(lev, GetArrOfConstPtrs(b[lev]),
                                     MLMG::Location::FaceCentroid);

      temperature_matrix->setLevelBC(lev, GetVecOfConstPtrs(Tg_arg)[lev]);
    }

    MLMG solver(*temperature_matrix);

    solver.apply(gradient_aux, Tg_arg);

    for(int lev = 0; lev <= finest_level; lev++) {
      MultiFab::Copy(*gradient[lev], *gradient_aux[lev], 0, 0, 1, 1);
      EB_set_covered(*gradient[lev], 0, gradient[lev]->nComp(), gradient[lev]->nGrow(), 0.);
    }

    for (int lev = 0; lev <= finest_level; lev++) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.growntilebox(IntVect(1,1,1));

        if (bx.ok()) {

          Array4<Real const> dummy_arr;

          Array4<Real      > const& gradient_array = gradient[lev]->array(mfi);
          Array4<Real const> const& ep_g_array     = ep_g[lev]->const_array(mfi);
          Array4<Real const> const& ro_g_array     = ro_g[lev]->const_array(mfi);
          Array4<Real const> const& T_g_array      = Tg_arg[lev]->const_array(mfi);
          Array4<Real const> const& X_gk_array     = fluid_is_a_mixture ? X_gk[lev]->const_array(mfi) : dummy_arr;

          amrex::ParallelFor(bx, [gradient_array,ep_g_array,T_g_array,ro_g_array,
              X_gk_array,fluid_is_a_mixture,nspecies_g,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const Real ep_ro_g = ep_g_array(i,j,k)*ro_g_array(i,j,k);
            const Real Tg = T_g_array(i,j,k);

            Real partial_hg(0.);

            if (!fluid_is_a_mixture) {
              partial_hg = fluid_parms.calc_partial_h_g<RunOn::Gpu>(Tg);
            } else {
              for (int n_g(0); n_g < nspecies_g; ++n_g) {
                partial_hg += X_gk_array(i,j,k,n_g)*fluid_parms.calc_partial_h_gk<RunOn::Gpu>(Tg,n_g);
              }
            }

            gradient_array(i,j,k) += ep_ro_g * partial_hg;
          });
        }
      }

      gradient[lev]->FillBoundary(geom[lev].periodicity());
      EB_set_covered(*gradient[lev], 0, gradient[lev]->nComp(), gradient[lev]->nGrow(), 0.);
      delete gradient_aux[lev];
    }

    return;
  };

  // **************************************************************************
  // Solve the nonlinear equation
  // **************************************************************************

  DampedNewton::NormMF norm0 = [&] (const Vector<MultiFab*>& vec_of_MFs) -> Real
  {
    Vector<Real> vec_of_norms(vec_of_MFs.size(), 0.);
    std::transform(vec_of_MFs.begin(), vec_of_MFs.end(), vec_of_norms.begin(),
                   [](MultiFab* MF) { return MF->norm0(); });

    Real norm(0);
    for (auto value: vec_of_norms)
      norm = amrex::max(norm, value);

    return norm;
  };

  // **************************************************************************
  // **************************************************************************
  // **************************************************************************

  // Damped Newton solution
  try {
    DampedNewton::DampingFactor damping_factor(0., .25);
    DampedNewton::solve(T_g, R, partial_R, norm0, damping_factor,
                        abstol, reltol, maxiter);

  } catch (std::exception& first_exc) {

    first_exc.what();

    try {
      DampedNewton::DampingFactor damping_factor(0., .5);
      DampedNewton::solve(T_g, R, partial_R, norm0, damping_factor,
                          10*abstol, 10*reltol, maxiter);

    } catch (std::exception& second_exc) {

      second_exc.what();

      try {
        DampedNewton::DampingFactor damping_factor(1., .5);
        DampedNewton::solve(T_g, R, partial_R, norm0, damping_factor,
                            100*abstol, 100*reltol, maxiter);
      } catch (std::exception& third_exc) {

        if (ParallelDescriptor::IOProcessor()) {
          third_exc.what();
        }

        amrex::Abort("Damped-Newton solver did not converge");
      }
    }
  }

  // **************************************************************************
  // **************************************************************************
  // **************************************************************************

  for(int lev = 0; lev <= finest_level; lev++) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi) {

      const EBFArrayBox& epg_fab = static_cast<EBFArrayBox const&>((*ep_g[lev])[mfi]);
      const EBCellFlagFab& flags = epg_fab.getEBCellFlagFab();

      Box const& bx = mfi.growntilebox({1,1,1});

      Array4<Real const> dummy_arr;
      Array4<Real      > const& h_g_array  = h_g[lev]->array(mfi);
      Array4<Real const> const& T_g_array  = T_g[lev]->const_array(mfi);
      Array4<Real const> const& X_gk_array = fluid_is_a_mixture ? X_gk[lev]->const_array(mfi) : dummy_arr;

      auto const& flags_arr = flags.const_array();

      amrex::ParallelFor(bx, [h_g_array,T_g_array,X_gk_array,fluid_parms,
          fluid_is_a_mixture,nspecies_g,run_on_device,flags_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

        const Real Tg = T_g_array(i,j,k);

        Real hg(0.);

        if (!fluid_is_a_mixture) {
          hg = run_on_device ?
            fluid_parms.calc_h_g<RunOn::Device>(Tg, cell_is_covered) :
            fluid_parms.calc_h_g<RunOn::Host>(Tg, cell_is_covered);
        } else {
          for (int n_g(0); n_g < nspecies_g; ++n_g) {
            const Real Xgk = X_gk_array(i,j,k,n_g);
            hg += run_on_device ?
              Xgk*fluid_parms.calc_h_gk<RunOn::Device>(Tg, n_g, cell_is_covered) :
              Xgk*fluid_parms.calc_h_gk<RunOn::Host>(Tg, n_g, cell_is_covered);
          }
        }

        h_g_array(i,j,k) = hg;
      });
    }

    h_g[lev]->FillBoundary(geom[lev].periodicity());
  }
}
