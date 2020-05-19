#include <DiffusionOp.H>

#include <MFIX_FLUID_Parms.H>

using namespace amrex;

//
// Implicit solve for scalar diffusion
//
void DiffusionOp::diffuse_temperature (      Vector< MultiFab* >  T_g_in,
                                       const Vector< MultiFab* > ep_g_in,
                                       const Vector< MultiFab* > ro_g_in,
                                       const Vector< MultiFab* >  h_g_in,
                                       const Vector< MultiFab* > cp_g_in,
                                       const Vector< MultiFab* >  k_g_in,
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
    //      b: eta

    if(verbose > 0)
        amrex::Print() << "Diffusing temperature ..." << std::endl; 

    // Set alpha and beta
    scal_matrix->setScalars(1.0, dt);

    Vector<BCRec> bcs_s; // This is just to satisfy the call to EB_interp...
    bcs_s.resize(3);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab ep_g_k_g(ep_g_in[lev]->boxArray(), ep_g_in[lev]->DistributionMap(),
            1, 1, MFInfo(), ep_g_in[lev]->Factory());
        // Initialize to 0
        ep_g_k_g.setVal(0.);

        MultiFab::Copy(ep_g_k_g, *ep_g_in[lev], 0, 0, 1, 1);
        MultiFab::Multiply(ep_g_k_g, *k_g_in[lev], 0, 0, 1, 1);

        EB_interp_CellCentroid_to_FaceCentroid (ep_g_k_g, GetArrOfPtrs(b[lev]), 0, 0, 1, geom[lev], bcs_s);

        // Turn "ep_g_in" into (rho * ep_g * cp_g)
        MultiFab::Multiply((*ep_g_in[lev]), (*ro_g_in[lev]), 0, 0, 1, 0);
        MultiFab::Multiply((*ep_g_in[lev]), (*cp_g_in[lev]), 0, 0, 1, 0);

        // This sets the coefficients
        scal_matrix->setACoeffs (lev, (*ep_g_in[lev]));
        scal_matrix->setBCoeffs (lev, GetArrOfConstPtrs(b[lev]),MLMG::Location::FaceCentroid);

        // Zero these out just to have a clean start because they have 3 components
        //      (due to re-use with velocity solve)
        phi[lev]->setVal(0.0);
        rhs[lev]->setVal(0.0);

        // Set the right hand side to equal rhs
        MultiFab::Copy((*rhs[lev]),(*T_g_in[lev]), 0, 0, 1, 0);

        // Multiply rhs by (rho * ep_g * cp_g) -- we are solving
        //
        //      rho ep_g c_p T_star = rho ep_g c_p T_old + dt ( - div (rho u ep_g h) + rho div (ep_g k_g (grad T)) )
        //
        MultiFab::Multiply((*rhs[lev]), (*ep_g_in[lev]), 0, 0, 1, 0);

        // Turn "ep_g_in" back into (ep_g)
        MultiFab::Divide((*ep_g_in[lev]), (*ro_g_in[lev]), 0, 0, 1, 0);
        MultiFab::Divide((*ep_g_in[lev]), (*cp_g_in[lev]), 0, 0, 1, 0);

        MultiFab::Copy(*phi[lev],*T_g_in[lev], 0, 0, 1, 1);
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
        MultiFab::Copy(*T_g_in[lev], *phi[lev], 0, 0, 1, 1);
        MultiFab::Copy(*h_g_in[lev], *phi[lev], 0, 0, 1, 1);
        MultiFab::Multiply(*h_g_in[lev], *cp_g_in[lev], 0, 0, 1, 1);
    }
}
