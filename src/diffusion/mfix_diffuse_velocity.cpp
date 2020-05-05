#include <AMReX_MultiFabUtil.H>
#include <DiffusionOp.H>

using namespace amrex;

//
// Implicit tensor solve for velocity diffusion
//
void DiffusionOp::diffuse_velocity (Vector< MultiFab* > vel_in,
                                    const Vector< MultiFab* > ep_ro_in,
                                    const Vector< MultiFab* > eta_in,
                                    Real dt)
{
    BL_PROFILE("DiffusionOp::diffuse_velocity");

    int finest_level = amrcore->finestLevel();

    // Update the coefficients of the matrix going into the solve based on the
    // current state of the simulation. Recall that the relevant matrix is
    //
    //      alpha a - beta div ( b grad )   <--->   rho - dt div ( eta grad )
    //
    // So the constants and variable coefficients are:
    //
    //      alpha: 1
    //      beta: dt
    //      a: ro
    //      b: eta

    // Set alpha and beta
    vel_matrix->setScalars(1.0, dt);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Compute the spatially varying b coefficients (on faces) to equal the
        // apparent viscosity
        average_cellcenter_to_face(GetArrOfPtrs(b[lev]), *eta_in[lev], geom[lev]);
        
        // This sets the coefficients
        vel_matrix->setACoeffs(lev, (*ep_ro_in[lev]));
        vel_matrix->setShearViscosity  (lev, GetArrOfConstPtrs(b[lev]));
        vel_matrix->setEBShearViscosity(lev, (*eta_in[lev]));
    }

    if(verbose > 0)
        amrex::Print() << "Diffusing velocity components all together..." << std::endl; 

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Set the right hand side to equal rho
        MultiFab::Copy((*rhs[lev]),(*vel_in[lev]), 0, 0, 3, 0);

        // Multiply rhs by rho to get momentum
        // Note that vel holds the updated velocity:
        //
        //      u_old + dt ( - u grad u + div ( eta (grad u)^T ) / rho - grad p / rho + gravity )
        //
        for (int i = 0; i < 3; i++)
           MultiFab::Multiply((*rhs[lev]), (*ep_ro_in[lev]), 0, i, 1, 0);

        // By this point we must have filled the Dirichlet values of phi stored in ghost cells
        MultiFab::Copy(*phi[lev],*vel_in[lev], 0, 0, 3, 1);
        phi[lev]->FillBoundary(geom[lev].periodicity());
        vel_matrix->setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);

        // matrix->setEBHomogDirichlet(lev, *eta_in[lev]);
    }

    MLMG solver(*vel_matrix);
    setSolverSettings(solver);

    // This ensures that ghost cells of sol are correctly filled when returned from the solver
    solver.setFinalFillBC(true);

    solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        phi[lev]->FillBoundary(geom[lev].periodicity());
        MultiFab::Copy(*vel_in[lev], *phi[lev], 0, 0, AMREX_SPACEDIM, 1);
    }

    if(verbose > 0)
        amrex::Print() << " Done diffusing all velocity components" << std::endl;
}
