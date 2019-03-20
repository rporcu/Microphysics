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
   mfix_compute_bcoeff_diff();

   //
   // First define the matrix (operator).
   // Class MLABecLaplacian describes the following operator:
   //
   //       (alpha * a - beta * (del dot b grad)) sol
   //
   LPInfo                       info;
   MLEBABecLap matrix(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));
   Vector<const MultiFab*>      tmp;
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
      tmp = amrex::GetVecOfConstPtrs ( bcoeff_diff[lev] ) ;
      b_tmp[0] = tmp[0];
      b_tmp[1] = tmp[1];
      b_tmp[2] = tmp[2];

      // This sets the spatially varying A coefficients
      matrix.setACoeffs ( lev, (*rop_g[lev]) );

      // This sets the spatially varying b coefficients
      matrix.setBCoeffs ( lev, b_tmp );

      // This sets the coefficient on the wall and defines it as a homogeneous Dirichlet bc for the solve.
      matrix.setEBHomogDirichlet (lev);

      // This tells the solver to use the higher order extrapolation to define d(phi)/dn at EB walls
      // This may not be robust in the presence of small cells so it is an option, not required
      //     (but does get Poiseuille flow right in the presence of walls at cell boundaries)
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
      solver.setMaxIter (mg_max_iter);
      solver.setMaxFmgIter (mg_max_fmg_iter);
      solver.setCGMaxIter (mg_cg_maxiter);

      // By this point we must have filled the Dirichlet values of sol stored in the ghost cells
      for (int lev = 0; lev < nlev; lev++)
      {
         rhs_diff[lev]->copy(*vel_g[lev],i,0,1,nghost,nghost);
         phi_diff[lev]->copy(*vel_g[lev],i,0,1,nghost,nghost);

         EB_set_covered(*phi_diff[lev], 0, phi_diff[lev]->nComp(), phi_diff[lev]->nGrow(), covered_val);
         phi_diff[lev] -> FillBoundary (geom[lev].periodicity());

         matrix.setLevelBC ( lev, GetVecOfConstPtrs(phi_diff)[lev] );

         // Define RHS = (rop) * (vel_g)
         MultiFab::Multiply((*rhs_diff[lev]), (*rop_g[lev]), 0, 0, 1, rhs_diff[lev]->nGrow());
      }

      // This ensures that ghost cells of sol are correctly filled when returned from the solver
      solver.setFinalFillBC(true);

      //
      // Finally, solve the system
      //
      //  (1 - div dot mu grad) u = RHS
      //
      solver.solve ( GetVecOfPtrs(phi_diff), GetVecOfConstPtrs(rhs_diff), mg_rtol, mg_atol );

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

   // Directions
   int xdir = 1;
   int ydir = 2;
   int zdir = 3;

   for (int lev = 0; lev < nlev; lev++)
   {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*mu_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         // Tileboxes for staggered components
         Box ubx = mfi.tilebox (e_x);
         Box vbx = mfi.tilebox (e_y);
         Box wbx = mfi.tilebox (e_z);

         // X direction
         compute_bcoeff_diff (BL_TO_FORTRAN_BOX(ubx),
                              BL_TO_FORTRAN_ANYD((*(bcoeff_diff[lev][0]))[mfi]),
                              BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]), &xdir );

         // Y direction
         compute_bcoeff_diff (BL_TO_FORTRAN_BOX(vbx),
                              BL_TO_FORTRAN_ANYD((*(bcoeff_diff[lev][1]))[mfi]),
                              BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]), &ydir );

         // Z direction
         compute_bcoeff_diff (BL_TO_FORTRAN_BOX(wbx),
                              BL_TO_FORTRAN_ANYD((*(bcoeff_diff[lev][2]))[mfi]),
                              BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]), &zdir );
   
      }

      bcoeff_diff[lev][0] -> FillBoundary(geom[lev].periodicity());
      bcoeff_diff[lev][1] -> FillBoundary(geom[lev].periodicity());
      bcoeff_diff[lev][2] -> FillBoundary(geom[lev].periodicity());
   }
}
