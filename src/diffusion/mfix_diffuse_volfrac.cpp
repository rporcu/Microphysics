#include <DiffusionOp.H>

using namespace amrex;

//
// Implicit solve for scalar diffusion
//
void DiffusionOp::diffuse_volfrac (Vector< MultiFab* >& scal_in,
                                   Real dcoeff)
{
    BL_PROFILE("DiffusionOp::diffuse_volfrac");

    int finest_level = amrcore->finestLevel();

    // Update the coefficients of the matrix going into the solve based on the current state of the
    // simulation. Recall that the relevant matrix is
    //
    //      alpha a - beta div ( b grad )   <--->   rho - dt div ( mu_s grad )
    //
    // So the constants and variable coefficients are:
    //
    //      alpha: 1
    //      beta: 1
    //      a: 1
    //      b: dcoeff

    // Set alpha and beta
    scal_matrix->setScalars(1.0, 1.0);

    int ntrac = 1;

    for(int lev = 0; lev <= finest_level; lev++)
    {
        for(int dir = 0; dir < 3; dir++)
           for(int n = 0; n < ntrac; n++)
             b[lev][dir]->setVal(dcoeff,1);

        for(int dir = 0; dir < 3; dir++)
            b[lev][dir]->FillBoundary(geom[lev].periodicity());

        // This sets the coefficients
        scal_matrix->setACoeffs (lev, 1.0);
        scal_matrix->setBCoeffs (lev, GetArrOfConstPtrs(b[lev]));
    }

    if(verbose > 0)
        amrex::Print() << "Diffusing volfrac ..." << std::endl; 

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Zero these out just to have a clean start because they have 3 components
        //      (due to re-use with velocity solve)
        phi[lev]->setVal(0.0);
        rhs[lev]->setVal(0.0);

        // Set the right hand side to equal rhs
        MultiFab::Copy((*rhs[lev]),(*scal_in[lev]), 0, 0, ntrac, 0);

        MultiFab::Copy(*phi[lev],*scal_in[lev], 0, 0, ntrac, 1);
        scal_matrix->setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);
    }

    MLMG solver(*scal_matrix);
    setSolverSettings(solver);

    // This ensures that ghost cells of sol are correctly filled when returned from the solver
    solver.setFinalFillBC(true);

    solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        phi[lev]->FillBoundary(geom[lev].periodicity());
        MultiFab::Copy(*scal_in[lev], *phi[lev], 0, 0, ntrac, 1);
    }
}
