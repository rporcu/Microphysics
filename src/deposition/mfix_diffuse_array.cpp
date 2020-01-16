#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLTensorOp.H>
#include <AMReX_MLEBTensorOp.H>

#include <MFIX_BC_Parms.H>
//
// Implicit tensor solve
//
void
mfix::mfix_diffuse_array (const amrex::Vector< std::unique_ptr<MultiFab> > & mf_to_diffuse,
                          amrex::Real dcoeff)
{
   BL_PROFILE("mfix::mfix_diffuse_drag");


   //
   // First define the operator "ebtensorop"
   //
   //       (alpha * a - beta * (del dot b grad)) sol
   //
   LPInfo info;
   info.setMaxCoarseningLevel(diff_mg_max_coarsening_level);
   MLEBTensorOp ebtensorop(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));

   // It is essential that we set MaxOrder of the solver to 2
   // if we want to use the standard sol(i)-sol(i-1) approximation
   // for the gradient at Dirichlet boundaries.
   // The solver's default order is 3 and this uses three points for the
   // gradient at a Dirichlet boundary.
   ebtensorop.setMaxOrder(2);

   // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
   ebtensorop.setDomainBC(BC::diff_vel_lobc,
                          BC::diff_vel_hibc);

   // Solving (1.0 * a_coeff - dt * div (mu grad)) phi = rhs
   ebtensorop.setScalars(1.0, 1.0);

   // Compute the coefficients
   for (int lev = 0; lev < nlev; lev++)
   {

     MultiFab dcoeff_mf(mu_g[lev]->boxArray(), mu_g[lev]->DistributionMap(), 1,
                        mu_g[lev]->nGrow(), MFInfo(), *ebfactory[lev]);

     dcoeff_mf.setVal(dcoeff);

     average_cellcenter_to_face(GetArrOfPtrs(bcoeff[lev]), dcoeff_mf, geom[lev]);

     bcoeff[lev][0]->FillBoundary(geom[lev].periodicity());
     bcoeff[lev][1]->FillBoundary(geom[lev].periodicity());
     bcoeff[lev][2]->FillBoundary(geom[lev].periodicity());

     ebtensorop.setShearViscosity(lev, GetArrOfConstPtrs(bcoeff[lev]));
     ebtensorop.setEBShearViscosity(lev, dcoeff_mf);

     // This sets the spatially varying A coefficients
     MultiFab a_coeff(ro_g[lev]->boxArray(), ro_g[lev]->DistributionMap(), 1,
                      ro_g[lev]->nGrow(), MFInfo(), *ebfactory[lev]);

     a_coeff.setVal(1.0);

     ebtensorop.setACoeffs(lev, a_coeff);
   }

   amrex::Print() << "Diffusing all velocity components together " << std::endl;
   mfix_print_max_vel(0);

   MLMG solver(ebtensorop);

   // Set the verbosity
   solver.setVerbose(diff_mg_verbose);
   solver.setCGVerbose(diff_mg_cg_verbose);

   // Set the max number of iterations
   solver.setMaxIter(diff_mg_maxiter);
   solver.setCGMaxIter(diff_mg_cg_maxiter);

   if (diff_bottom_solver_type == "smoother")
   {
      solver.setBottomSolver(MLMG::BottomSolver::smoother);
   }
   else if (diff_bottom_solver_type == "hypre")
   {
      solver.setBottomSolver(MLMG::BottomSolver::hypre);
   }

   // By this point we must have filled the Dirichlet values of sol stored in the ghost cells
   for (int lev = 0; lev < nlev; lev++)
   {
       MultiFab::Copy((*diff_phi4[lev]), (*mf_to_diffuse[lev]), 0, 0, 4, diff_phi4[lev]->nGrow());

       EB_set_covered(*diff_phi4[lev], 0, diff_phi4[lev]->nComp(), diff_phi4[lev]->nGrow(), covered_val);
       diff_phi4[lev]->FillBoundary(geom[lev].periodicity());

       ebtensorop.setLevelBC(lev, GetVecOfConstPtrs(diff_phi4)[lev]);

      // Define RHS = drag force
       MultiFab::Copy((*diff_rhs4[lev]),(*mf_to_diffuse[lev]), 0, 0, 4, 0);
   }

   // This ensures that ghost cells of sol are correctly filled when returned from the solver
   solver.setFinalFillBC(true);

   //
   // Finally, solve the system
   //
   //  (1.0 - div dot nu grad) drag = RHS
   //
   solver.solve(GetVecOfPtrs(diff_phi4), GetVecOfConstPtrs(diff_rhs4), diff_mg_rtol, diff_mg_atol);

   for (int lev = 0; lev < nlev; lev++)
   {
       diff_phi4[lev]->FillBoundary(geom[lev].periodicity());
       MultiFab::Copy(*mf_to_diffuse[lev], *diff_phi4[lev], 0, 0, 4, 1);
   }

   amrex::Print() << "After diffusing all drag components " << std::endl;

}
