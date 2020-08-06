#include <mfix_diffusion_op.H>

using namespace amrex;

//
// Implicit solve for species mass fraction
//
void DiffusionOp::diffuse_species (      Vector< MultiFab* >     X_gk_in,
                                   const Vector< MultiFab* > ep_ro_g_in,
                                   const Vector< MultiFab* >     D_gk_in,
                                   Real dt)
{
    BL_PROFILE("DiffusionOp::diffuse_species");

    int finest_level = amrcore->finestLevel();

    // Number of fluid species
    const int nspecies_g = X_gk_in[0]->nComp();

    for (int n(0); n < nspecies_g; n++)
    {
      Vector< MultiFab* > X_gk(finest_level+1);

      for(int lev = 0; lev <= finest_level; lev++)
      {
        X_gk[lev] = new MultiFab(X_gk_in[lev]->boxArray(), X_gk_in[lev]->DistributionMap(),
            1, 1, MFInfo(), X_gk_in[lev]->Factory());

        X_gk[lev]->setVal(0.0);

        MultiFab::Copy(*X_gk[lev], *X_gk_in[lev], n, 0, 1, 1);
      }

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

      // Set alpha and beta
      species_matrix->setScalars(1.0, dt);

      Vector<BCRec> bcs_X; // This is just to satisfy the call to EB_interp...
      bcs_X.resize(3);

      for(int lev = 0; lev <= finest_level; lev++)
      {
          MultiFab ep_ro_D_gk(ep_ro_g_in[lev]->boxArray(),
              ep_ro_g_in[lev]->DistributionMap(), 1, 1, MFInfo(),
              ep_ro_g_in[lev]->Factory());

          ep_ro_D_gk.setVal(0.);

          MultiFab::Copy(ep_ro_D_gk, *ep_ro_g_in[lev], 0, 0, 1, 1);
          MultiFab::Multiply(ep_ro_D_gk, *D_gk_in[lev], n, 0, 1, 1);

          EB_interp_CellCentroid_to_FaceCentroid (ep_ro_D_gk, GetArrOfPtrs(b[lev]), 0, 0, 1, geom[lev], bcs_X);

          // This sets the coefficients
          species_matrix->setACoeffs (lev, (*ep_ro_g_in[lev]));
          species_matrix->setBCoeffs (lev, GetArrOfConstPtrs(b[lev]), MLMG::Location::FaceCentroid);
      }

      if(verbose > 0)
        amrex::Print() << "Diffusing fluid mass fractions one at a time ..."
                       << std::endl;

      for(int lev = 0; lev <= finest_level; lev++)
      {
          // Zero these out just to have a clean start because they have 3 components
          //      (due to re-use with velocity solve)
          phi[lev]->setVal(0.0);
          rhs[lev]->setVal(0.0);

          // Set the right hand side to equal rhs
          MultiFab::Copy((*rhs[lev]), (*X_gk[lev]), 0, 0, 1, 0);

          // Multiply rhs by (rho * ep_g * D_gk) -- we are solving
          //
          //      rho ep_g X_star = rho ep_g T_old + dt ( - div (rho u ep_g X_gk) + rho div (ep_g D_gk (grad X_gk)) )
          //
          MultiFab::Multiply((*rhs[lev]), (*ep_ro_g_in[lev]), 0, 0, 1, 0);

          MultiFab::Copy(*phi[lev], *X_gk[lev], 0, 0, 1, 1);
          species_matrix->setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);
      }

      MLMG solver(*species_matrix);
      setSolverSettings(solver);

      // This ensures that ghost cells of sol are correctly filled when returned from the solver
      solver.setFinalFillBC(true);

      solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

      for(int lev = 0; lev <= finest_level; lev++)
      {
          phi[lev]->FillBoundary(geom[lev].periodicity());
          MultiFab::Copy(*X_gk_in[lev], *phi[lev], 0, n, 1, 1);
      }

      for (int lev(0); lev <= finest_level; lev++) {
        delete X_gk[lev];
      }
    }
}
