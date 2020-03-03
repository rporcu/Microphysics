#include <AMReX_MultiFabUtil.H>
#include <DiffusionOp.H>

using namespace amrex;

void DiffusionOp::diffuse_drag (Vector< MultiFab* > drag_in,
                                Real dcoeff)
{
    BL_PROFILE("DiffusionOp::diffuse_drag");

    int finest_level = amrcore->finestLevel();

    // Set alpha and beta
    vel_matrix->setScalars(1.0, 1.0);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Note that b doesn't need any ghost cells
        for(int dir = 0; dir < 3; dir++)
             b[lev][dir]->setVal(dcoeff);
        
        // This sets the coefficients
        vel_matrix->setACoeffs(lev, 1.0);
        vel_matrix->setShearViscosity  (lev, dcoeff);
        vel_matrix->setEBShearViscosity(lev, dcoeff);
    }

    if(verbose > 0)
        amrex::Print() << "Diffusing drag components all together..." << std::endl; 

    for(int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab::Copy((*rhs[lev]),(*drag_in[lev]), 0, 0, 3, 0);

        MultiFab::Copy(*phi[lev],*drag_in[lev], 0, 0, 3, 1);
        phi[lev]->FillBoundary(geom[lev].periodicity());
        vel_matrix->setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);
    }

    MLMG solver(*vel_matrix);
    setSolverSettings(solver);

    // This ensures that ghost cells of sol are correctly filled when returned from the solver
    solver.setFinalFillBC(true);

    solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        phi[lev]->FillBoundary(geom[lev].periodicity());
        MultiFab::Copy(*drag_in[lev], *phi[lev], 0, 0, AMREX_SPACEDIM, 1);
    }

    if(verbose > 0)
        amrex::Print() << " Done diffusing all drag components" << std::endl;
}
