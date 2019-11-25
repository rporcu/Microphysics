#include <diffusion_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLTensorOp.H>
#include <AMReX_MLEBTensorOp.H>

//
// Implicit scalar solve
//
void
mfix::mfix_diffuse_scalar (const amrex::Vector< std::unique_ptr<MultiFab> > & mf_to_diffuse,
                           amrex::Real dcoeff)
{
   BL_PROFILE("mfix::mfix_diffuse_scalar");

   // The boundary conditions need only be set once -- we do this at level 0
   int bc_lo[3], bc_hi[3];

   // Whole domain
   Box domain(geom[0].Domain());

   mfix_set_epg_bcs(mf_to_diffuse);

   // Set BCs for Poisson's solver
   set_scal_diff_bc (bc_lo, bc_hi,
                     domain.loVect(), domain.hiVect(),
                     &nghost,
                     bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                     bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                     bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

   //
   // First define the operator "ebscalarop"
   //
   //       (alpha * a - beta * (del dot b grad)) sol
   //
   LPInfo                       info;
   info.setMaxCoarseningLevel(diff_mg_max_coarsening_level);
   MLEBABecLap ebscalarop(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));

   // It is essential that we set MaxOrder of the solver to 2
   // if we want to use the standard sol(i)-sol(i-1) approximation
   // for the gradient at Dirichlet boundaries.
   // The solver's default order is 3 and this uses three points for the
   // gradient at a Dirichlet boundary.
   ebscalarop.setMaxOrder(2);

   // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
   ebscalarop.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                            {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );

   // Solving (1.0 * a_coeff - dt * div (mu grad)) phi = rhs
   ebscalarop.setScalars(1.0, 1.0);

   // Compute the coefficients
   for (int lev = 0; lev < nlev; lev++)
   {

       MultiFab dcoeff_mf( mu_g[lev]->boxArray(), mu_g[lev]->DistributionMap(), 1, mu_g[lev]->nGrow(),
                         MFInfo(), *ebfactory[lev]);

       dcoeff_mf.setVal(dcoeff);

       average_cellcenter_to_face( GetArrOfPtrs(bcoeff[lev]), dcoeff_mf, geom[lev] );

       bcoeff[lev][0] -> FillBoundary(geom[lev].periodicity());
       bcoeff[lev][1] -> FillBoundary(geom[lev].periodicity());
       bcoeff[lev][2] -> FillBoundary(geom[lev].periodicity());

       ebscalarop.setBCoeffs (lev, GetArrOfConstPtrs(bcoeff[lev]));

       // This sets the spatially varying A coefficients
       MultiFab a_coeff( ro_g[lev]->boxArray(), ro_g[lev]->DistributionMap(), 1, ro_g[lev]->nGrow(),
                         MFInfo(), *ebfactory[lev]);

       a_coeff.setVal(1.0);

       ebscalarop.setACoeffs ( lev, a_coeff );

       Gpu::synchronize();
   }

   amrex::Print() << "Diffusing solids volume fraction " << std::endl;

   MLMG  solver(ebscalarop);

   // Set the verbosity
   solver.setVerbose   (diff_mg_verbose);
   solver.setCGVerbose (diff_mg_cg_verbose);

   // Set the max number of iterations
   solver.setMaxIter (diff_mg_maxiter);
   solver.setCGMaxIter (diff_mg_cg_maxiter);

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
       MultiFab::Copy((*diff_phi1[lev]),(*mf_to_diffuse[lev]), 0, 0, 1, diff_phi1[lev]->nGrow());

       EB_set_covered(*diff_phi1[lev], 0, diff_phi1[lev]->nComp(), diff_phi1[lev]->nGrow(), covered_val);
       diff_phi1[lev] -> FillBoundary (geom[lev].periodicity());

       ebscalarop.setLevelBC ( lev, GetVecOfConstPtrs(diff_phi1)[lev] );

      // Define RHS = eps
       MultiFab::Copy((*diff_rhs1[lev]),(*mf_to_diffuse[lev]), 0, 0, 1, 0);
   }

   // This ensures that ghost cells of sol are correctly filled when returned from the solver
   solver.setFinalFillBC(true);

   //
   // Finally, solve the system
   //
   //  (1.0 - div dot nu grad) eps = RHS
   //
   solver.solve ( GetVecOfPtrs(diff_phi1), GetVecOfConstPtrs(diff_rhs1), diff_mg_rtol, diff_mg_atol );

   for (int lev = 0; lev < nlev; lev++)
   {
       diff_phi1[lev]->FillBoundary (geom[lev].periodicity());
       MultiFab::Copy( *mf_to_diffuse[lev], *diff_phi1[lev], 0, 0, 1, 1);
   }

   amrex::Print() << "After diffusing volume fraction " << std::endl;

}
