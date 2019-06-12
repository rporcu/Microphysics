#include <mfix_diff_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>
#include <AMReX_EBMultiFabUtil.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLTensorOp.H>
#include <AMReX_MLEBTensorOp.H>

//
// Implicit tensor solve
//
void
mfix::mfix_diffuse_velocity_tensor (amrex::Real time, amrex::Real dt)
{
   BL_PROFILE("mfix::mfix_diffuse_velocity_tensor");

   // Swap ghost cells and apply BCs to velocity
   mfix_set_velocity_bcs (time, 0);

   // The boundary conditions need only be set once -- we do this at level 0
   int bc_lo[3], bc_hi[3];

   // Whole domain
   Box domain(geom[0].Domain());

   // Set BCs for Poisson's solver
   set_diff_bc (bc_lo, bc_hi,
                domain.loVect(), domain.hiVect(),
                &nghost,
                bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

   //
   // First define the operator "ebtensorop"
   //
   //       (alpha * a - beta * (del dot b grad)) sol
   //
   LPInfo                       info;
   MLEBTensorOp ebtensorop(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));

   // It is essential that we set MaxOrder of the solver to 2
   // if we want to use the standard sol(i)-sol(i-1) approximation
   // for the gradient at Dirichlet boundaries.
   // The solver's default order is 3 and this uses three points for the
   // gradient at a Dirichlet boundary.
   ebtensorop.setMaxOrder(2);

   // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
   ebtensorop.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                            {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );

   // Solving (1.0 * a_coeff - dt * div (mu grad)) phi = rhs
   ebtensorop.setScalars(1.0, dt);

   // Compute the coefficients
   for (int lev = 0; lev < nlev; lev++)
   {
       average_cellcenter_to_face( GetArrOfPtrs(bcoeff_cc[lev]), *mu_g[lev], geom[lev] );

       bcoeff_cc[lev][0] -> FillBoundary(geom[lev].periodicity());
       bcoeff_cc[lev][1] -> FillBoundary(geom[lev].periodicity());
       bcoeff_cc[lev][2] -> FillBoundary(geom[lev].periodicity());

       ebtensorop.setShearViscosity  (lev, GetArrOfConstPtrs(bcoeff_cc[lev]));
       ebtensorop.setEBShearViscosity(lev, (*mu_g[lev]));

       // This sets the spatially varying A coefficients
       MultiFab a_coeff( ro_g[lev]->boxArray(), ro_g[lev]->DistributionMap(), 1, ro_g[lev]->nGrow(),
                         MFInfo(), *ebfactory[lev]);

       MultiFab::Copy    ( a_coeff, *ro_g[lev], 0, 0, 1, ro_g[lev]->nGrow() );
       MultiFab::Multiply( a_coeff, *ep_g[lev], 0, 0, 1, ep_g[lev]->nGrow() );
 
       ebtensorop.setACoeffs ( lev, a_coeff );
   }

   amrex::Print() << "Diffusing all velocity components together " << std::endl;
   mfix_print_max_vel(0);

   MLMG  solver(ebtensorop);

   // Set the verbosity
   solver.setVerbose   (diff_mg_verbose);
   solver.setCGVerbose (diff_mg_cg_verbose);

   // Set the max number of iterations
   solver.setMaxIter (diff_mg_max_iter);
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
      diff_rhs[lev]->copy(*vel_g[lev],0,0,3,nghost,nghost);
      diff_phi[lev]->copy(*vel_g[lev],0,0,3,nghost,nghost);

      EB_set_covered(*diff_phi[lev], 0, diff_phi[lev]->nComp(), diff_phi[lev]->nGrow(), covered_val);
      diff_phi[lev] -> FillBoundary (geom[lev].periodicity());

      ebtensorop.setLevelBC ( lev, GetVecOfConstPtrs(diff_phi)[lev] );

      // Define RHS = (ro_g) * (ep_g) * (vel_g)
      for (int i = 0; i < 3; i++)
      {
         MultiFab::Multiply((*diff_rhs[lev]), (*ro_g[lev]), 0, i, 1, diff_rhs[lev]->nGrow());
         MultiFab::Multiply((*diff_rhs[lev]), (*ep_g[lev]), 0, i, 1, diff_rhs[lev]->nGrow());
      }
   }

   // This ensures that ghost cells of sol are correctly filled when returned from the solver
   solver.setFinalFillBC(true);

   //
   // Finally, solve the system
   //
   //  ((ro_g)*(ep_g) - div dot mu grad) u = RHS
   //
   solver.solve ( GetVecOfPtrs(diff_phi), GetVecOfConstPtrs(diff_rhs), diff_mg_rtol, diff_mg_atol );

   for (int lev = 0; lev < nlev; lev++)
   {
      diff_phi[lev]->FillBoundary (geom[lev].periodicity());
      vel_g[lev]->copy(*diff_phi[lev],0,0,3,nghost,nghost);
   }

   amrex::Print() << "After diffusing all velocity components " << std::endl;
   mfix_print_max_vel(0);

   // Swap ghost cells and apply BCs to velocity
   mfix_set_velocity_bcs (time, 0);
}
