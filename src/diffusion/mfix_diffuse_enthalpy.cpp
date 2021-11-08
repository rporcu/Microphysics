#include <mfix_diffusion_op.H>

#include <mfix_eb_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_solvers.H>

#include <AMReX_EBFArrayBox.H>

using namespace amrex;
using namespace Solvers;


//
// Implicit solve for scalar diffusion
//
void DiffusionOp::diffuse_enthalpy (const Vector< MultiFab* >& h_g,
                                    const Vector< MultiFab* >& ep_g,
                                    const Vector< MultiFab* >& ro_g,
                                    const Vector< MultiFab* >& T_g,
                                    const Vector< MultiFab* >& X_gk,
                                    const Vector< MultiFab* >& T_g_on_eb,
                                    Real dt)
{
  BL_PROFILE("DiffusionOp::diffuse_enthalpy");

  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

  int finest_level = amrcore->finestLevel();

  // Update the coefficients of the matrix going into the solve based on the
  // current state of the simulation. Recall that the relevant matrix is
  //
  //      alpha a - beta div ( b grad )   <--->   rho - dt div ( k_g grad )
  //
  // So the constants and variable coefficients are:
  //
  //      alpha: 1
  //      beta: dt
  //      a: ro_g ep_g cp_g
  //      b: ep_g (k_g / cp_g)

  if(verbose > 0)
    amrex::Print() << "Diffusing enthalpy ..." << std::endl;

  // Set alpha and beta
  temperature_matrix->setScalars(1.0, dt);

  auto& fluid_parms = *fluid.parameters;

  const auto fluid_is_a_mixture = fluid.is_a_mixture;

  Vector<BCRec> bcs_dummy; // This is just to satisfy the call to EB_interp...
  bcs_dummy.resize(1);

  for(int lev = 0; lev <= finest_level; lev++)
  {
    MultiFab ep_k_Icp_g(ep_g[lev]->boxArray(), ep_g[lev]->DistributionMap(),
                    ep_g[lev]->nComp(), 1, //ep_g[lev]->nGrow(),
                    MFInfo(), ep_g[lev]->Factory());

    // Initialize to 0
    ep_k_Icp_g.setVal(0.);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.growntilebox(IntVect(1,1,1));

      if (bx.ok())
      {
        Array4<Real const> dummy_arr;

        Array4<Real      > const& ep_k_Icp_g_array = ep_k_Icp_g.array(mfi);
        Array4<Real const> const& ep_g_array       = ep_g[lev]->const_array(mfi);
        Array4<Real const> const& T_g_array        = T_g[lev]->const_array(mfi);
        Array4<Real const> const& X_gk_array       = fluid_is_a_mixture ?
          X_gk[lev]->const_array(mfi) : dummy_arr;

        amrex::ParallelFor(bx, [ep_g_array,X_gk_array,T_g_array,ep_k_Icp_g_array,
            fluid_parms,fluid_is_a_mixture,run_on_device]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real cp_g(0);

          if (!fluid_is_a_mixture)
            cp_g = run_on_device ?
              fluid_parms.calc_cp_g<RunOn::Device>(T_g_array(i,j,k)) :
              fluid_parms.calc_cp_g<RunOn::Host>(T_g_array(i,j,k));
          else
            for(int n(0); n<fluid_parms.m_nspecies; ++n)
              cp_g += run_on_device ?
                X_gk_array(i,j,k,n)*fluid_parms.calc_cp_gk<RunOn::Device>(T_g_array(i,j,k), n) :
                X_gk_array(i,j,k,n)*fluid_parms.calc_cp_gk<RunOn::Host>(T_g_array(i,j,k), n);

          ep_k_Icp_g_array(i,j,k) = ep_g_array(i,j,k)*fluid_parms.calc_k_g(T_g_array(i,j,k)) / cp_g;
        });
      }
    }

    EB_interp_CellCentroid_to_FaceCentroid (ep_k_Icp_g, GetArrOfPtrs(b[lev]), 0, 0, 1, geom[lev], bcs_dummy);

    // Turn "ep_g" into (rho * ep_g * cp_g)
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.tilebox();

      Array4<Real      > const& ep_g_array  = ep_g[lev]->array(mfi);
      Array4<Real const> const& ro_g_array  = ro_g[lev]->const_array(mfi);
      Array4<Real const> const& T_g_array   = T_g[lev]->const_array(mfi);

      amrex::ParallelFor(bx, [ep_g_array,T_g_array,ro_g_array,fluid_parms,fluid_is_a_mixture]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { 
        ep_g_array(i,j,k) *= ro_g_array(i,j,k);
      });
    }

    // This sets the coefficients
    temperature_matrix->setACoeffs (lev, (*ep_g[lev]));

    temperature_matrix->setBCoeffs (lev, GetArrOfConstPtrs(b[lev]),
        MLMG::Location::FaceCentroid);

    // Zero these out just to have a clean start because they have 3 components
    //      (due to re-use with velocity solve)
    phi[lev]->setVal(0.0);
    rhs[lev]->setVal(0.0);

    // Set the right hand side to equal rhs
    MultiFab::Copy((*rhs[lev]), (*h_g[lev]), 0, 0, 1, 0);

    // Multiply rhs by (rho * ep_g * cp_g) -- we are solving
    //
    //   rho ep_g c_p T_star = rho ep_g c_p T_old + dt ( - div (rho u ep_g h)
    //   + rho div (ep_g k_g (grad T)) )
    //
    MultiFab::Multiply((*rhs[lev]), (*ep_g[lev]), 0, 0, 1, 0);

    // Turn "ep_g" back into (ep_g)
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.tilebox();

      Array4<Real      > const& ep_g_array = ep_g[lev]->array(mfi);
      Array4<Real const> const& ro_g_array = ro_g[lev]->const_array(mfi);
      Array4<Real const> const& T_g_array  = T_g[lev]->const_array(mfi);

      amrex::ParallelFor(bx, [ep_g_array,T_g_array,ro_g_array,fluid_parms]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        ep_g_array(i,j,k) /= ro_g_array(i,j,k);
      });
    }

    if (EB::fix_temperature) {
      // The following is a WIP in AMReX
      //temperature_matrix->setPhiOnCentroid();

      MultiFab k_Icp_g(T_g_on_eb[lev]->boxArray(), T_g_on_eb[lev]->DistributionMap(),
                       T_g_on_eb[lev]->nComp(), T_g_on_eb[lev]->nGrow(), MFInfo(),
                       T_g_on_eb[lev]->Factory());

      k_Icp_g.setVal(0);

      MultiFab h_g_on_eb(T_g_on_eb[lev]->boxArray(), T_g_on_eb[lev]->DistributionMap(),
                         T_g_on_eb[lev]->nComp(), T_g_on_eb[lev]->nGrow(), MFInfo(),
                         T_g_on_eb[lev]->Factory());

      h_g_on_eb.setVal(0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*T_g_on_eb[lev]); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.growntilebox(T_g_on_eb[lev]->nGrowVect());

        const EBFArrayBox& T_g_on_eb_fab = static_cast<EBFArrayBox const&>((*T_g_on_eb[lev])[mfi]);
        const EBCellFlagFab& flags = T_g_on_eb_fab.getEBCellFlagFab();

        if (bx.ok()) {
          Array4<Real const> dummy_arr;

          Array4<Real      > const& k_Icp_g_array   = k_Icp_g.array(mfi);
          Array4<Real      > const& h_g_on_eb_array = h_g_on_eb.array(mfi);
          Array4<Real const> const& X_gk_array      = fluid_is_a_mixture ?
            X_gk[lev]->const_array(mfi) : dummy_arr;
          Array4<Real const> const& T_g_on_eb_array = T_g_on_eb[lev]->const_array(mfi);

          auto const& flags_arr = flags.const_array();

          amrex::ParallelFor(bx, [k_Icp_g_array,h_g_on_eb_array,X_gk_array,
              T_g_on_eb_array,fluid_parms,fluid_is_a_mixture,run_on_device,
              flags_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

            const Real Tg = T_g_on_eb_array(i,j,k);

            if (Tg > 0) {
              Real cp_g(0);
              Real h_g(0);

              if (!fluid_is_a_mixture) {

                cp_g = run_on_device ?
                  fluid_parms.calc_cp_g<RunOn::Device>(Tg) :
                  fluid_parms.calc_cp_g<RunOn::Host>(Tg);

                h_g = run_on_device ?
                  fluid_parms.calc_h_g<RunOn::Device>(Tg, cell_is_covered) :
                  fluid_parms.calc_h_g<RunOn::Host>(Tg, cell_is_covered);

              } else {

                for(int n(0); n < fluid_parms.m_nspecies; ++n) {

                  cp_g += run_on_device ?
                    X_gk_array(i,j,k,n)*fluid_parms.calc_cp_gk<RunOn::Device>(Tg, n) :
                    X_gk_array(i,j,k,n)*fluid_parms.calc_cp_gk<RunOn::Host>(Tg, n);

                  h_g += run_on_device ?
                    X_gk_array(i,j,k,n)*fluid_parms.calc_h_gk<RunOn::Device>(Tg, n, cell_is_covered) :
                    X_gk_array(i,j,k,n)*fluid_parms.calc_h_gk<RunOn::Host>(Tg, n, cell_is_covered);

                }
              }

              k_Icp_g_array(i,j,k) = fluid_parms.calc_k_g(Tg) / cp_g;
              h_g_on_eb_array(i,j,k) = h_g;
            }
          });
        }
      }

      k_Icp_g.FillBoundary(geom[lev].periodicity());

      temperature_matrix->setEBDirichlet(lev, h_g_on_eb, k_Icp_g);
    }

    MultiFab::Copy(*phi[lev], *h_g[lev], 0, 0, 1, 1);

    temperature_matrix->setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);
  }

  MLMG solver(*temperature_matrix);
  setSolverSettings(solver);

  // This ensures that ghost cells of sol are correctly filled when returned
  // from the solver
  solver.setFinalFillBC(true);

  solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

  for(int lev = 0; lev <= finest_level; lev++)
  {
    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(T_g[lev]->Factory());
    const auto& flags = factory.getMultiEBCellFlagFab();
    const auto& volfrac = factory.getVolFrac();

    phi[lev]->FillBoundary(geom[lev].periodicity());
    MultiFab::Copy(*h_g[lev], *phi[lev], 0, 0, 1, 1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.growntilebox({1,1,1});

      Array4<Real const> dummy_arr;

      Array4<Real const> const& h_g_array  = h_g[lev]->const_array(mfi);
      Array4<Real      > const& T_g_array  = T_g[lev]->array(mfi);
      Array4<Real const> const& X_gk_array = fluid_is_a_mixture ?
        X_gk[lev]->const_array(mfi) : dummy_arr;
      Array4<Real const> const& ep_g_array = ep_g[lev]->const_array(mfi);

      auto const& flags_arr = flags.const_array(mfi);
      auto const& volfrac_arr = volfrac.const_array(mfi);

      amrex::ParallelFor(bx, [h_g_array,T_g_array,X_gk_array,fluid_parms,
          fluid_is_a_mixture,run_on_device,ep_g_array,flags_arr,volfrac_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

        if (!cell_is_covered) {

          const Real epg_loc = ep_g_array(i,j,k);
          const Real vfrac = volfrac_arr(i,j,k);

          // Residual computation
          auto R = [&] AMREX_GPU_DEVICE (Real Tg_arg) {
            
            Real hg_loc(0);

            if (!fluid_is_a_mixture) {

              hg_loc = run_on_device ?
                fluid_parms.calc_h_g<RunOn::Device>(Tg_arg, cell_is_covered) :
                fluid_parms.calc_h_g<RunOn::Host>(Tg_arg, cell_is_covered);

            } else {

              for (int n(0); n < fluid_parms.m_nspecies; ++n) {
                const Real h_gk = run_on_device ?
                  fluid_parms.calc_h_gk<RunOn::Device>(Tg_arg, n, cell_is_covered) :
                  fluid_parms.calc_h_gk<RunOn::Host>(Tg_arg, n, cell_is_covered);

                hg_loc += X_gk_array(i,j,k,n)*h_gk;
              }
            }

            return hg_loc - h_g_array(i,j,k);
          };

          // Partial derivative computation
          auto partial_R = [&] AMREX_GPU_DEVICE (Real Tg_arg) {
            
            Real gradient(0);

            if (!fluid_is_a_mixture) {

              gradient = run_on_device ?
                fluid_parms.calc_partial_h_g<RunOn::Device>(Tg_arg) :
                fluid_parms.calc_partial_h_g<RunOn::Host>(Tg_arg);

            } else {

              for (int n(0); n < fluid_parms.m_nspecies; ++n) {
                const Real h_gk = run_on_device ?
                  fluid_parms.calc_partial_h_gk<RunOn::Device>(Tg_arg,n) :
                  fluid_parms.calc_partial_h_gk<RunOn::Host>(Tg_arg,n);

                gradient += X_gk_array(i,j,k,n)*h_gk;
              }
            }

            return gradient;
          };

          Real Tg(T_g_array(i,j,k));

          int solver_iterations(0);

          {
            DampedNewton::DampingFactor damping_factor(0., 0.);
            solver_iterations =
              DampedNewton::solve(Tg, R, partial_R, damping_factor(epg_loc, vfrac),
                  1.e-8, 1.e-8, 500);

          } if (solver_iterations == 500) {

            DampedNewton::DampingFactor damping_factor(1., 0.);
            solver_iterations =
              DampedNewton::solve(Tg, R, partial_R, damping_factor(epg_loc, vfrac),
                  1.e-7, 1.e-7, 500);

          } if (solver_iterations == 500) {

            DampedNewton::DampingFactor damping_factor(1., 1.);
            solver_iterations =
              DampedNewton::solve(Tg, R, partial_R, damping_factor(epg_loc, vfrac),
                  1.e-6, 1.e-6, 500);

          } if (solver_iterations == 500) {
            amrex::Abort("Damped-Newton solver did not converge");
          }

          T_g_array(i,j,k) = Tg;
        }
      });
    }
  }
}
