#include <DiffusionOp.H>

#include <MFIX_EB_Parms.H>
#include <MFIX_FLUID_Parms.H>

using namespace amrex;

//
// Implicit solve for scalar diffusion
//
void DiffusionOp::diffuse_temperature (Vector< MultiFab* > T_g,
                                       const Vector< MultiFab* > ep_g,
                                       const Vector< MultiFab* > ro_g,
                                       const Vector< MultiFab* > h_g,
                                       const Vector< MultiFab* > cp_g,
                                       const Vector< MultiFab* > k_g,
                                       const Vector< MultiFab* > T_g_on_eb,
                                       const Vector< MultiFab* > k_g_on_eb,
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

  Vector<BCRec> bcs_s; // This is just to satisfy the call to EBterp...
  bcs_s.resize(3);

  for(int lev = 0; lev <= finest_level; lev++)
  {
    MultiFab ep_g_k_g(ep_g[lev]->boxArray(), ep_g[lev]->DistributionMap(), 1, 1,
        MFInfo(), ep_g[lev]->Factory());

    // Initialize to 0
    ep_g_k_g.setVal(0.);

    MultiFab::Copy(ep_g_k_g, *ep_g[lev], 0, 0, 1, 1);
    MultiFab::Multiply(ep_g_k_g, *k_g[lev], 0, 0, 1, 1);

    EB_interp_CellCentroid_to_FaceCentroid (ep_g_k_g, GetArrOfPtrs(b[lev]), 0,
        0, 1, geom[lev], bcs_s);

    // Turn "ep_g" into (rho * ep_g * cp_g)
    MultiFab::Multiply((*ep_g[lev]), (*ro_g[lev]), 0, 0, 1, 0);
    MultiFab::Multiply((*ep_g[lev]), (*cp_g[lev]), 0, 0, 1, 0);

    // This sets the coefficients
    temperature_matrix->setACoeffs (lev, (*ep_g[lev]));

    temperature_matrix->setBCoeffs (lev, GetArrOfConstPtrs(b[lev]),
        MLMG::Location::FaceCentroid);

    // Zero these out just to have a clean start because they have 3 components
    //      (due to re-use with velocity solve)
    phi[lev]->setVal(0.0);
    rhs[lev]->setVal(0.0);

    // Set the right hand side to equal rhs
    MultiFab::Copy((*rhs[lev]),(*T_g[lev]), 0, 0, 1, 0);

    // Multiply rhs by (rho * ep_g * cp_g) -- we are solving
    //
    //   rho ep_g c_p T_star = rho ep_g c_p T_old + dt ( - div (rho u ep_g h)
    //   + rho div (ep_g k_g (grad T)) )
    //
    MultiFab::Multiply((*rhs[lev]), (*ep_g[lev]), 0, 0, 1, 0);

    // Turn "ep_g" back into (ep_g)
    MultiFab::Divide((*ep_g[lev]), (*ro_g[lev]), 0, 0, 1, 0);
    MultiFab::Divide((*ep_g[lev]), (*cp_g[lev]), 0, 0, 1, 0);

    if (EB::fix_temperature) {
      // The following is a WIP in AMReX
      //temperature_matrix->setPhiOnCentroid();
      temperature_matrix->setEBDirichlet(lev, *T_g_on_eb[lev], *k_g_on_eb[lev]);
    }

    MultiFab::Copy(*phi[lev],*T_g[lev], 0, 0, 1, 1);

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
    MultiFab::Copy(*h_g[lev], *phi[lev], 0, 0, 1, 1);
    MultiFab::Multiply(*h_g[lev], *cp_g[lev], 0, 0, 1, 1);
  }
}
