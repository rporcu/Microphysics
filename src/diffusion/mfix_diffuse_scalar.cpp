#include <DiffusionOp.H>

using namespace amrex;

//
// Implicit solve for scalar diffusion
//
void DiffusionOp::diffuse_scalar (Vector< MultiFab* > scal_in,
                                  const Vector< MultiFab* > ep_ro_in,
                                  const Vector< Real > mu_s,
                                  Real dt)
{
    BL_PROFILE("DiffusionOp::diffuse_scalar");

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
    //      a: ro
    //      b: eta

    // Set alpha and beta
    scal_matrix->setScalars(1.0, dt);

    int ntrac = scal_in[0]->nComp();

    for(int lev = 0; lev <= finest_level; lev++)
    {
        for(int dir = 0; dir < 3; dir++)
           for(int n = 0; n < ntrac; n++)
             b[lev][dir]->setVal(mu_s[n],n,1);

        // This sets the coefficients
        scal_matrix->setACoeffs (lev, (*ep_ro_in[lev]));
        scal_matrix->setBCoeffs (lev, GetArrOfConstPtrs(b[lev]),MLMG::Location::FaceCentroid);
    }

    if(verbose > 0)
        amrex::Print() << "Diffusing tracers one at a time ..." << std::endl; 

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Zero these out just to have a clean start because they have 3 components
        //      (due to re-use with velocity solve)
        phi[lev]->setVal(0.0);
        rhs[lev]->setVal(0.0);

        // Set the right hand side to equal rhs
        MultiFab::Copy((*rhs[lev]),(*scal_in[lev]), 0, 0, ntrac, 0);

        // Multiply rhs by rho  -- we are solving
        //
        //      rho s_star = rho s_old + dt ( - div (rho u s) + rho div (mu (grad s)) )
        //
        for(int n = 0; n < ntrac; n++)
           MultiFab::Multiply((*rhs[lev]), (*ep_ro_in[lev]), 0, n, 1, 0);

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
