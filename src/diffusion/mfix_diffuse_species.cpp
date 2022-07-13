#include <mfix_diffusion_op.H>
#include <mfix_fluid.H>

using namespace amrex;

//
// Implicit solve for species mass fraction
//
void DiffusionOp::diffuse_species (const Vector< MultiFab* >&    X_gk_in,
                                   const Vector< Array< MultiFab*, AMREX_SPACEDIM>>& J_gk,
                                   const Vector< MultiFab* >& ep_ro_g_in,
                                   const Vector< MultiFab* >& /*T_g_in*/,
                                         Vector<BCRec> const& species_h_bcrec,
                                   Real dt)
{
    BL_PROFILE("DiffusionOp::diffuse_species");

    int finest_level = amrcore->finestLevel();

    // Update the coefficients of the matrix going into the solve based on the current state of the
    // simulation. Recall that the relevant matrix is
    //
    //      alpha a - beta div ( b grad )   <--->   rho - dt div ( mu_s grad )
    //
    // So the constants and variable coefficients are:
    //
    //      alpha: 1
    //      beta: dt
    //      a: ro_g ep_g
    //      b: ro_g ep_g D_gk

    if(verbose > 0)
      amrex::Print() << "Diffusing species mass fractions ..." << std::endl;

    // Set alpha and beta
    species_matrix->setScalars(1.0, dt);

    // Number of fluid species
    const int nspecies_g = fluid.nspecies();

    for(int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab ep_ro_D_gk(ep_ro_g_in[lev]->boxArray(),
            ep_ro_g_in[lev]->DistributionMap(), nspecies_g, 1, MFInfo(),
            ep_ro_g_in[lev]->Factory());

        ep_ro_D_gk.setVal(0.);

        const auto& fluid_parms = fluid.parameters();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*ep_ro_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.growntilebox(IntVect(1,1,1));

          if (bx.ok())
          {
            Array4<Real      > const& ep_ro_D_gk_arr = ep_ro_D_gk.array(mfi);
            Array4<Real const> const& ep_ro_g_arr    = ep_ro_g_in[lev]->const_array(mfi);

            amrex::ParallelFor(bx, [ep_ro_g_arr,ep_ro_D_gk_arr,nspecies_g,
                fluid_parms]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              const Real ep_ro_g = ep_ro_g_arr(i,j,k);
              const Real val = ep_ro_g*fluid_parms.get_D_g();

              for (int n(0); n < nspecies_g; ++n) {
                ep_ro_D_gk_arr(i,j,k,n) = val;
              }
            });
          }
        }

        EB_interp_CellCentroid_to_FaceCentroid (ep_ro_D_gk, GetArrOfPtrs(species_b[lev]),
                                                0, 0, nspecies_g, geom[lev], species_h_bcrec);

        // This sets the coefficients
        species_matrix->setACoeffs (lev, (*ep_ro_g_in[lev]));
        species_matrix->setBCoeffs (lev, GetArrOfConstPtrs(species_b[lev]),
                                    MLMG::Location::FaceCentroid);

        // Zero these out just to have a clean start because they have 3 components
        //      (due to re-use with velocity solve)
        species_phi[lev]->setVal(0.0);
        species_rhs[lev]->setVal(0.0);

        // Set rhs equal to X_gk and
        // Multiply rhs by (rho * ep_g * D_gk) -- we are solving
        //
        //      rho ep_g X_star = rho ep_g T_old + dt ( - div (rho u ep_g X_gk) + rho div (ep_g D_gk (grad X_gk)) )
        //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*X_gk_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.growntilebox(IntVect(0));

          if (bx.ok())
          {
            Array4<Real const> const& ep_ro_g_arr = ep_ro_g_in[lev]->const_array(mfi);
            Array4<Real const> const& X_gk_arr    = X_gk_in[lev]->const_array(mfi);
            Array4<Real      > const& rhs_arr     = species_rhs[lev]->array(mfi);

            amrex::ParallelFor(bx, [ep_ro_g_arr,X_gk_arr,rhs_arr,nspecies_g]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              const Real ep_ro_g = ep_ro_g_arr(i,j,k);

              for (int n(0); n < nspecies_g; ++n)
                rhs_arr(i,j,k,n) = X_gk_arr(i,j,k,n)*ep_ro_g;
            });
          }
        }

        MultiFab::Copy(*species_phi[lev], *X_gk_in[lev], 0, 0, nspecies_g, 1);
        species_matrix->setLevelBC(lev, GetVecOfConstPtrs(species_phi)[lev]);
    }

    MLMG solver(*species_matrix);
    setSolverSettings(solver);

    // This ensures that ghost cells of sol are correctly filled when returned from the solver
    solver.setFinalFillBC(true);

    solver.solve(GetVecOfPtrs(species_phi), GetVecOfConstPtrs(species_rhs), mg_rtol, mg_atol);

    solver.getFluxes(J_gk, GetVecOfPtrs(species_phi), MLMG::Location::FaceCentroid);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab::Copy(*X_gk_in[lev], *species_phi[lev], 0, 0, nspecies_g, 1);
    }
}
