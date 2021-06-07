#include <mfix_diffusion_op.H>

#include <mfix_eb_parms.H>
#include <mfix_fluid_parms.H>

using namespace amrex;

//
// Implicit solve for scalar diffusion
//
void DiffusionOp::diffuse_temperature (const Vector< MultiFab* >& T_g,
                                       const Vector< MultiFab* >& ep_g,
                                       const Vector< MultiFab* >& ro_g,
                                       const Vector< MultiFab* >& h_g,
                                       const Vector< MultiFab* >& T_g_on_eb,
                                       Real dt)
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
  //      alpha: 1
  //      beta: dt
  //      a: ro_g ep_g cp_g
  //      b: ep_g k_g

  if(verbose > 0)
    amrex::Print() << "Diffusing temperature ..." << std::endl;

  // Set alpha and beta
  temperature_matrix->setScalars(1.0, dt);

  auto& fluid_parms = *fluid.parameters;

  Vector<BCRec> bcs_dummy; // This is just to satisfy the call to EB_interp...
  bcs_dummy.resize(1);

  for(int lev = 0; lev <= finest_level; lev++)
  {
    MultiFab ep_k_g(ep_g[lev]->boxArray(), ep_g[lev]->DistributionMap(),
                    ep_g[lev]->nComp(), 1, //ep_g[lev]->nGrow(),
                    MFInfo(), ep_g[lev]->Factory());

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
        Array4<Real const> const& T_g_array    = T_g[lev]->const_array(mfi);

        amrex::ParallelFor(bx, [ep_g_array,T_g_array,ep_k_g_array,fluid_parms]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          ep_k_g_array(i,j,k) = ep_g_array(i,j,k)*fluid_parms.calc_k_g(T_g_array(i,j,k));
        });
      }
    }

    EB_interp_CellCentroid_to_FaceCentroid (ep_k_g, GetArrOfPtrs(b[lev]), 0, 0, 1, geom[lev], bcs_dummy);

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

      amrex::ParallelFor(bx, [ep_g_array,T_g_array,ro_g_array,fluid_parms]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        ep_g_array(i,j,k) *= ro_g_array(i,j,k)*fluid_parms.calc_cp_g(T_g_array(i,j,k));
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
    MultiFab::Copy((*rhs[lev]), (*T_g[lev]), 0, 0, 1, 0);

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

      Array4<Real      > const& ep_g_array  = ep_g[lev]->array(mfi);
      Array4<Real const> const& ro_g_array  = ro_g[lev]->const_array(mfi);
      Array4<Real const> const& T_g_array   = T_g[lev]->const_array(mfi);

      amrex::ParallelFor(bx, [ep_g_array,T_g_array,ro_g_array,fluid_parms]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        ep_g_array(i,j,k) /= ro_g_array(i,j,k)*fluid_parms.calc_cp_g(T_g_array(i,j,k));
      });
    }

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
        Box const& bx = mfi.growntilebox(T_g_on_eb[lev]->nGrowVect());

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

    MultiFab::Copy(*phi[lev], *T_g[lev], 0, 0, 1, 1);

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
    phi[lev]->FillBoundary(geom[lev].periodicity());
    MultiFab::Copy(*T_g[lev], *phi[lev], 0, 0, 1, 1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.growntilebox({1,1,1});

      Array4<Real      > const& h_g_array  = h_g[lev]->array(mfi);
      Array4<Real const> const& T_g_array  = T_g[lev]->const_array(mfi);

      amrex::ParallelFor(bx, [h_g_array,T_g_array,fluid_parms]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        h_g_array(i,j,k) = fluid_parms.calc_h_g(T_g_array(i,j,k));
      });
    }
  }
}
