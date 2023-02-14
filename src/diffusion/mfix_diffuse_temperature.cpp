#include <mfix_diffusion_op.H>

#include <mfix.H>
#include <mfix_eb.H>
#include <mfix_fluid.H>
#include <mfix_solvers.H>

#include <AMReX_EB_utils.H>

using namespace amrex;
using namespace MFIXSolvers;

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
  //      A: ro_g ep_g cp_g
  //      B: ep_g k_g

  if(verbose > 0)
    amrex::Print() << "Diffusing temperature ..." << std::endl;

  const auto& fluid_parms = fluid.parameters();
  const int fluid_is_a_mixture = fluid.isMixture();
  const int nspecies_g = fluid.nspecies();

  amrex::Vector<amrex::MultiFab*> A(finest_level+1, nullptr);

  for (int lev(0); lev <= finest_level; ++lev) {
    A[lev] = new amrex::MultiFab(grids[lev], dmap[lev], 1, 1, MFInfo(),
                                 *ebfactory[lev]);
    A[lev]->setVal(0.);
  }

  // **************************************************************************
  // Define the norm function
  // **************************************************************************

  auto norm0 = [&] (const Vector<MultiFab*>& vec_of_MFs) -> Real
  {
    Vector<Real> vec_of_norms(vec_of_MFs.size(), 0.);
    std::transform(vec_of_MFs.begin(), vec_of_MFs.end(), vec_of_norms.begin(),
      [](MultiFab* MF) { return MF->norm0(0, 0, false, true); });

    Real norm(0);
    for (auto value: vec_of_norms)
      norm = amrex::max(norm, value);

    ParallelDescriptor::ReduceRealMax(norm);

    return norm;
  };

  // **************************************************************************
  // Do the Newton iterations
  // **************************************************************************

  amrex::Vector<amrex::MultiFab*> update(finest_level+1, nullptr);
//  amrex::Vector<amrex::MultiFab*> residue(finest_level+1, nullptr);

  for (int lev(0); lev <= finest_level; ++lev) {
    update[lev] = new amrex::MultiFab(grids[lev], dmap[lev], 1, nghost,
                                      MFInfo(), *ebfactory[lev]);
    update[lev]->setVal(0.);

//    residue[lev] = new amrex::MultiFab(grids[lev], dmap[lev], 1, nghost,
//                                       MFInfo(), *ebfactory[lev]);
//    residue[lev]->setVal(0.);
  }

  int iter(0);
  const amrex::Real update_tol = reltol*norm0(T_g) + abstol;
  //const amrex::Real residue_tol = reltol*norm0(residue) + abstol;

  Vector<MultiFab*> k_g_on_eb(finest_level+1, nullptr);

  if (m_embedded_boundaries.fix_temperature()) {

    for (int lev = 0; lev <= finest_level; ++lev) {

      // The following is a WIP in AMReX
      //temperature_matrix->setPhiOnCentroid();

      k_g_on_eb[lev] = new MultiFab(T_g_on_eb[lev]->boxArray(),
                                    T_g_on_eb[lev]->DistributionMap(),
                                    T_g_on_eb[lev]->nComp(), T_g_on_eb[lev]->nGrow(),
                                    MFInfo(), T_g_on_eb[lev]->Factory());

      k_g_on_eb[lev]->setVal(0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*T_g_on_eb[lev]); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.growntilebox(IntVect(1,1,1));

        if (bx.ok()) {
          Array4<Real      > const& k_g_on_eb_array = k_g_on_eb[lev]->array(mfi);
          Array4<Real const> const& T_g_on_eb_array = T_g_on_eb[lev]->const_array(mfi);

          amrex::ParallelFor(bx, [k_g_on_eb_array,T_g_on_eb_array,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            if (T_g_on_eb_array(i,j,k) > 0)
              k_g_on_eb_array(i,j,k) = fluid_parms.calc_k_g(T_g_on_eb_array(i,j,k));
          });
        }
      }
    }
  }

  do {

    // Set alpha and beta
    temperature_matrix->setScalars(1.0, dt);

    Vector<BCRec> bcs_dummy; // This is just to satisfy the call to EB_interp...
    bcs_dummy.resize(3);

    for(int lev = 0; lev <= finest_level; lev++) {

      A[lev]->setVal(0.);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.tilebox();

        const EBFArrayBox& epg_fab = static_cast<EBFArrayBox const&>((*ep_g[lev])[mfi]);
        const EBCellFlagFab& flags = epg_fab.getEBCellFlagFab();

        if (bx.ok())
        {
          Array4<Real const> dummy_arr;

          Array4<Real      > const& A_array       = A[lev]->array(mfi);
          Array4<Real const> const& ep_g_array    = ep_g[lev]->const_array(mfi);
          Array4<Real const> const& ro_g_array    = ro_g[lev]->const_array(mfi);
          Array4<Real const> const& T_g_array     = T_g[lev]->const_array(mfi);
          Array4<Real const> const& X_gk_array    = fluid_is_a_mixture ?
            X_gk[lev]->const_array(mfi) : dummy_arr;

          auto const& flags_arr = flags.const_array();

          amrex::ParallelFor(bx, [ep_g_array,T_g_array,ro_g_array,fluid_parms,
              A_array,X_gk_array,dt,fluid_is_a_mixture,nspecies_g,flags_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

            if (!cell_is_covered) {

              const Real T_g_loc = T_g_array(i,j,k);
              
              Real cp_g_loc(0);

              if (!fluid_is_a_mixture) {

                cp_g_loc = fluid_parms.calc_cp_g<run_on>(T_g_loc);

              } else {

                for (int n(0); n < nspecies_g; ++n) {
                  
                  const Real cp_gk = fluid_parms.calc_cp_gk<run_on>(T_g_loc, n);

                  cp_g_loc += X_gk_array(i,j,k,n)*cp_gk;
                }
              }

              A_array(i,j,k) = ep_g_array(i,j,k)*ro_g_array(i,j,k)*cp_g_loc;
            }
          });
        }
      }

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
          Array4<Real      > const& ep_k_g_array  = ep_k_g.array(mfi);
          Array4<Real const> const& ep_g_array    = ep_g[lev]->const_array(mfi);
          Array4<Real const> const& T_g_array     = T_g[lev]->const_array(mfi);

          amrex::ParallelFor(bx, [ep_g_array,T_g_array,ep_k_g_array,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const Real T_g_loc = T_g_array(i,j,k);

            ep_k_g_array(i,j,k) = ep_g_array(i,j,k)*fluid_parms.calc_k_g(T_g_loc);
          });
        }
      }

      EB_interp_CellCentroid_to_FaceCentroid (ep_k_g, GetArrOfPtrs(b[lev]), 0,
                                              0, 1, geom[lev], bcs_dummy);

      if (m_embedded_boundaries.fix_temperature()) {
        temperature_matrix->setEBDirichlet(lev, *T_g_on_eb[lev], *k_g_on_eb[lev]);
      }

      temperature_matrix->setACoeffs(lev, *A[lev]);
      temperature_matrix->setBCoeffs(lev, GetArrOfConstPtrs(b[lev]),
                                     MLMG::Location::FaceCentroid);

      phi[lev]->setVal(0.);
      rhs[lev]->setVal(0.);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.tilebox();

        const EBFArrayBox& epg_fab = static_cast<EBFArrayBox const&>((*ep_g[lev])[mfi]);
        const EBCellFlagFab& flags = epg_fab.getEBCellFlagFab();

        if (bx.ok())
        {
          Array4<Real const> dummy_arr;

          Array4<Real      > const& rhs_array     = rhs[lev]->array(mfi);
          Array4<Real const> const& ep_g_array    = ep_g[lev]->const_array(mfi);
          Array4<Real const> const& ro_g_array    = ro_g[lev]->const_array(mfi);
          Array4<Real const> const& T_g_array     = T_g[lev]->const_array(mfi);
          Array4<Real const> const& h_g_array     = h_g[lev]->const_array(mfi);
          Array4<Real const> const& X_gk_array    = fluid_is_a_mixture ?
            X_gk[lev]->const_array(mfi) : dummy_arr;

          auto const& flags_arr = flags.const_array();

          amrex::ParallelFor(bx, [ep_g_array,T_g_array,ro_g_array,fluid_parms,
              X_gk_array,h_g_array,dt,fluid_is_a_mixture,nspecies_g,rhs_array,flags_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

            if (!cell_is_covered) {
              const Real T_g_loc = T_g_array(i,j,k);
              const Real ep_ro_g_loc = ep_g_array(i,j,k)*ro_g_array(i,j,k);
              
              Real cp_g_loc(0);
              Real h_g_loc(0);

              if (!fluid_is_a_mixture) {

                cp_g_loc = fluid_parms.calc_cp_g<run_on>(T_g_loc);
                h_g_loc = fluid_parms.calc_h_g<run_on>(T_g_loc);

              } else {

                for (int n(0); n < nspecies_g; ++n) {
                  
                  const Real cp_gk = fluid_parms.calc_cp_gk<run_on>(T_g_loc, n);
                  const Real h_gk = fluid_parms.calc_h_gk<run_on>(T_g_loc, n);

                  cp_g_loc += X_gk_array(i,j,k,n)*cp_gk;
                  h_g_loc += X_gk_array(i,j,k,n)*h_gk;
                }
              }

              rhs_array(i,j,k) = ep_ro_g_loc*cp_g_loc*T_g_loc -
                                 ep_ro_g_loc*h_g_loc +
                                 ep_ro_g_loc*h_g_array(i,j,k);
            }
          });
        }
      }

      MultiFab::Copy(*phi[lev], *T_g[lev], 0, 0, 1, 1);
      temperature_matrix->setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);
    } // end of loop on lev

    MLMG solver(*temperature_matrix);
    setSolverSettings(solver);

    // This ensures that ghost cells of sol are correctly filled when returned
    // from the solver
    solver.setFinalFillBC(true);

    solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

    // TODO
    // Here compute update = phi - T_g

    for(int lev = 0; lev <= finest_level; lev++) {
      phi[lev]->FillBoundary(geom[lev].periodicity());
    }

    for(int lev = 0; lev <= finest_level; lev++) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        const EBFArrayBox& epg_fab = static_cast<EBFArrayBox const&>((*ep_g[lev])[mfi]);
        const EBCellFlagFab& flags = epg_fab.getEBCellFlagFab();

        if (bx.ok())
        {
          Array4<Real      > const& update_array = update[lev]->array(mfi);
          Array4<Real      > const& T_g_array    = T_g[lev]->array(mfi);
          Array4<Real const> const& phi_array    = phi[lev]->const_array(mfi);

          auto const& flags_arr = flags.const_array();

          amrex::ParallelFor(bx, [update_array,T_g_array,phi_array,flags_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

            if (!cell_is_covered) {
              update_array(i,j,k) = phi_array(i,j,k) - T_g_array(i,j,k);
              T_g_array(i,j,k) = phi_array(i,j,k);
            }
          });
        }
      }
    } // end of loop on lev

    iter++;
    if (iter > maxiter) {
      Print() << "Newton solver iterations = " << iter << "\n";
      Print() << "Newton solver update = " << norm0(update) << "\n";
      amrex::Abort("Newton solver did not converge");
    }

  } while (//(norm0(residue) > residue_tol) ||
           (norm0(update) > update_tol));

  for (int lev(0); lev <= finest_level; ++lev) {
    delete A[lev];
    delete update[lev];
//    delete residue[lev];
  }

  if (m_embedded_boundaries.fix_temperature()) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      delete k_g_on_eb[lev];
    }
  }

  // **************************************************************************
  // **************************************************************************
  // **************************************************************************

  // TODO: update this during the Newton iterations
  for(int lev = 0; lev <= finest_level; lev++) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi) {

      const EBFArrayBox& epg_fab = static_cast<EBFArrayBox const&>((*ep_g[lev])[mfi]);
      const EBCellFlagFab& flags = epg_fab.getEBCellFlagFab();

      Box const& bx = mfi.growntilebox(IntVect(1,1,1));

      Array4<Real const> dummy_arr;

      Array4<Real      > const& h_g_array  = h_g[lev]->array(mfi);
      Array4<Real const> const& T_g_array  = T_g[lev]->const_array(mfi);
      Array4<Real const> const& X_gk_array = fluid_is_a_mixture ?
        X_gk[lev]->const_array(mfi) : dummy_arr;

      auto const& flags_arr = flags.const_array();

      amrex::ParallelFor(bx, [h_g_array,T_g_array,X_gk_array,fluid_parms,
          fluid_is_a_mixture,nspecies_g,flags_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

        const Real Tg = T_g_array(i,j,k);

        Real hg(0.);

        if (!fluid_is_a_mixture) {
          hg = fluid_parms.calc_h_g<run_on>(Tg, cell_is_covered);
        } else {
          for (int n_g(0); n_g < nspecies_g; ++n_g) {
            const Real Xgk = X_gk_array(i,j,k,n_g);
            hg += Xgk*fluid_parms.calc_h_gk<run_on>(Tg, n_g, cell_is_covered);
          }
        }

        h_g_array(i,j,k) = hg;
      });
    }
  }
}
