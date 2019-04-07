#include <AMReX_ParmParse.H>

#include <mfix_diff_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

#include <AMReX_EBMultiFabUtil.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>

//
// Implicit diffusion
//
void
mfix::mfix_diffuse_velocity (amrex::Real time, amrex::Real dt)

{
   BL_PROFILE("mfix::mfix_diffuse_velocity");

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

   // Compute the coefficients
   mfix_compute_bcoeff_diff ();

   //
   // First define the matrix (operator).
   // Class MLABecLaplacian describes the following operator:
   //
   //       (alpha * a - beta * (del dot b grad)) sol
   //
   LPInfo                       info;
   MLEBABecLap matrix(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));
   array<MultiFab const*,AMREX_SPACEDIM>   b_tmp;

   // It is essential that we set MaxOrder of the solver to 2
   // if we want to use the standard sol(i)-sol(i-1) approximation
   // for the gradient at Dirichlet boundaries.
   // The solver's default order is 3 and this uses three points for the
   // gradient at a Dirichlet boundary.
   matrix.setMaxOrder(2);

   // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
   matrix.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                        {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );

   // This sets alpha = 1 and beta = dt
   matrix.setScalars ( 1.0, dt );

   // Copy the PPE coefficient into the proper data strutcure
   for (int lev = 0; lev < nlev; lev++)
   {
      // This sets the spatially varying A coefficients
      MultiFab a_coeff( ro_g[lev]->boxArray(), ro_g[lev]->DistributionMap(), 1, ro_g[lev]->nGrow(),
                        MFInfo(), *ebfactory[lev]);

      MultiFab::Copy    ( a_coeff, *ro_g[lev], 0, 0, 1, ro_g[lev]->nGrow() );
      MultiFab::Multiply( a_coeff, *ep_g[lev], 0, 0, 1, ep_g[lev]->nGrow() );

      matrix.setACoeffs ( lev, a_coeff );

      // This sets the spatially varying b coefficients
      matrix.setBCoeffs ( lev, GetArrOfConstPtrs(bcoeff_diff[lev]) );

      // This sets the coefficient on the wall and defines it as a homogeneous Dirichlet bc for the solve.
      matrix.setEBHomogDirichlet ( lev, (*mu_g[lev]) );

      // This tells the solver to use the higher order extrapolation to define d(phi)/dn at EB walls
      // This may not be robust in the presence of small cells so it is an option, not required
      //     (but does get Poiseuille flow right in the presence of walls at cell boundaries)
      if (eb_ho_dirichlet == 1)
         matrix.setEBHODirichlet ( );
   }

   // Loop over the velocity components
   for (int i = 0; i < 3; i++)
   {
      amrex::Print() << "Diffusing velocity component " << i << std::endl;

      MLMG  solver(matrix);

      // Set the verbosity
      solver.setVerbose   (diff_mg_verbose);
      solver.setCGVerbose (diff_mg_cg_verbose);

      // Set the max number of iterations
      solver.setMaxIter (diff_mg_max_iter);
      solver.setMaxFmgIter (diff_mg_max_fmg_iter);
      solver.setCGMaxIter (diff_mg_cg_maxiter);

      // By this point we must have filled the Dirichlet values of sol stored in the ghost cells
      for (int lev = 0; lev < nlev; lev++)
      {
         rhs_diff[lev]->copy(*vel_g[lev],i,0,1,nghost,nghost);
         phi_diff[lev]->copy(*vel_g[lev],i,0,1,nghost,nghost);

         EB_set_covered(*phi_diff[lev], 0, phi_diff[lev]->nComp(), phi_diff[lev]->nGrow(), covered_val);
         phi_diff[lev] -> FillBoundary (geom[lev].periodicity());

         matrix.setLevelBC ( lev, GetVecOfConstPtrs(phi_diff)[lev] );

         // Define RHS = (ro_g) * (ep_g) * (vel_g)
         MultiFab::Multiply((*rhs_diff[lev]), (*ro_g[lev]), 0, 0, 1, rhs_diff[lev]->nGrow());
         MultiFab::Multiply((*rhs_diff[lev]), (*ep_g[lev]), 0, 0, 1, rhs_diff[lev]->nGrow());
      }

      // This ensures that ghost cells of sol are correctly filled when returned from the solver
      solver.setFinalFillBC(true);

      //
      // Finally, solve the system
      //
      //  (1 - div dot mu grad) u = RHS
      //
      solver.solve ( GetVecOfPtrs(phi_diff), GetVecOfConstPtrs(rhs_diff), diff_mg_rtol, diff_mg_atol );

      for (int lev = 0; lev < nlev; lev++)
      {
         phi_diff[lev] -> FillBoundary (geom[lev].periodicity());
         vel_g[lev]->copy(*phi_diff[lev],0,i,1,nghost,nghost);
      }
   }

   // Swap ghost cells and apply BCs to velocity
   mfix_set_velocity_bcs (time, 0);
}

//
// Computes bcoeff = mu_g at the faces of the scalar cells
//
void
mfix::mfix_compute_bcoeff_diff ()
{
   BL_PROFILE("mfix::mfix_compute_bcoeff_diff");

   for (int lev = 0; lev < nlev; lev++)
   {
       average_cellcenter_to_face( GetArrOfPtrs(bcoeff_diff[lev]), *mu_g[lev], geom[lev] );

      bcoeff_diff[lev][0] -> FillBoundary(geom[lev].periodicity());
      bcoeff_diff[lev][1] -> FillBoundary(geom[lev].periodicity());
      bcoeff_diff[lev][2] -> FillBoundary(geom[lev].periodicity());
   }
}
