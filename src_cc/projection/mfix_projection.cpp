#include <AMReX_ParmParse.H>

#include <mfix_proj_F.H>
#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLNodeLaplacian.H>

void
mfix_level::EvolveFluidProjection(int lev, int nstep, int steady_state, Real& dt,  Real& time, Real stop_time )
{
    BL_PROFILE_REGION_START("mfix::EvolveFluidProjection");
    BL_PROFILE("mfix::EvolveFluidProjection");

    amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";
    // Extrapolate boundary values for density and volume fraction
    // The subsequent call to mfix_set_scalar_bcs will only overwrite
    // rop_g and ep_g ghost values for PINF and POUT
    fill_mf_bc ( lev, *rop_g[lev] );
    fill_mf_bc ( lev, *ep_g[lev] );
    fill_mf_bc ( lev, *mu_g[lev] );

    // Fill ghost nodes and reimpose boundary conditions
    mfix_set_scalar_bcs (lev);
    mfix_set_velocity_bcs (lev,0);
    
    //
    // Start loop: if we are not seeking a steady state solution,
    // the loop will execute only once
    //
    int keep_looping = 1;
    int iter = 1;

    do
    {
        mfix_compute_dt(lev, time, stop_time, steady_state, dt, nodal_pressure);

        if (steady_state)
        {
           amrex::Print() << "\n   Iteration " << iter << " with dt = " << dt << "\n" << std::endl;
        } else {
           amrex::Print() << "\n   Step " << nstep+1 << ": from old_time " \
                          << time << " to new time " << time+dt
                          << " with dt = " << dt << "\n" << std::endl;
        }
 
        // Back up field
        // Backup field variable to old
	MultiFab::Copy (*ep_go[lev],  *ep_g[lev],  0, 0,  ep_g[lev]->nComp(),  ep_go[lev]->nGrow());
	MultiFab::Copy ( *p_go[lev],   *p_g[lev],  0, 0,   p_g[lev]->nComp(),   p_go[lev]->nGrow());
	MultiFab::Copy (*ro_go[lev],  *ro_g[lev],  0, 0,  ro_g[lev]->nComp(),  ro_go[lev]->nGrow());
	MultiFab::Copy (*rop_go[lev], *rop_g[lev], 0, 0, rop_g[lev]->nComp(), rop_go[lev]->nGrow());
	MultiFab::Copy (*vel_go[lev], *vel_g[lev], 0, 0, vel_g[lev]->nComp(), vel_go[lev]->nGrow());

        // User hooks
        for (MFIter mfi(*ep_g[lev], true); mfi.isValid(); ++mfi)
    	   mfix_usr2();

	//
	// Time integration step
	//
	
	// Calculate drag coefficient
	if (solve_dem)
	    mfix_calc_drag_fluid(lev);
	
        // Predictor step 
        mfix_apply_predictor ( lev, dt );

	// Print info about predictor step
        {
	    amrex::Print() << "\nAfter predictor step:\n";
    	    mfix_print_max_vel (lev);
    	    mfix_compute_diveu (lev);
	    amrex::Print() << "max(abs(diveu)) = " << diveu[lev] -> norm0 () << "\n";
        }

	// Calculate drag coefficient
	if (solve_dem)
	    mfix_calc_drag_fluid(lev);
	
	// Corrector step 
	mfix_apply_corrector ( lev, dt );

	// Print info about corrector step
        {
	    amrex::Print() << "\nAfter corrector step:\n";
    	    mfix_print_max_vel (lev);
    	    mfix_compute_diveu (lev);
	    amrex::Print() << "max(abs(diveu)) = " << diveu[lev] -> norm0 () << "\n";
        }
	    
	// 
        // Check whether to exit the loop or not
	// 
	if (steady_state) {
	    keep_looping = !steady_state_reached ( lev, dt );
	} else {
	    keep_looping = 0;
	}

	// Update interations count
	++iter;
    }
    while ( keep_looping );

    BL_PROFILE_REGION_STOP("mfix::EvolveFluidProjection");
}

void
mfix_level::mfix_project_velocity (int lev)
{
    // Project velocity field to make sure initial velocity is divergence-free
    Real dummy_dt = 1.0;
    mfix_apply_projection ( lev, dummy_dt );
}

void
mfix_level::mfix_initial_iterations (int lev, Real dt, Real stop_time, int steady_state)
{
    // Fill ghost cells
    mfix_set_scalar_bcs (lev);
    mfix_set_velocity_bcs (lev,0);

    // Copy vel_g into vel_go
    MultiFab::Copy (*vel_go[lev], *vel_g[lev],   0, 0, 3, vel_go[lev]->nGrow());

    Real time = 0.0;
    mfix_compute_dt(lev,time,stop_time,steady_state,dt,nodal_pressure);

    amrex::Print() << "Doing initial pressure iterations with dt = " << dt << std::endl;
    
    // Calculate drag coefficient
    if (solve_dem)
        mfix_calc_drag_fluid(lev);

    // Compute fluid acceleration (convection + diffusion)
    mfix_compute_velocity_slopes ( lev, vel_go );
    mfix_compute_fluid_acceleration ( lev, vel_go );

    for (int iter = 0; iter < 3; ++iter)
    {
       // First add the fluid acceleration
       MultiFab::Saxpy (*vel_g[lev], dt, *acc[lev], 0, 0, vel_g[lev]->nComp(), 0);

       // Add the forcing terms
       mfix_apply_forcing_terms ( lev, dt, vel_g );

       mfix_add_pressure_gradient ( lev, -dt);

       // Compute intermediate velocity
       if (solve_dem)
          mfix_compute_intermediate_velocity ( lev, dt );

       // Exchange halo nodes and apply BCs to velocity
       mfix_set_velocity_bcs (lev,0);

       // Project velocity field
       mfix_apply_projection ( lev, dt );

       // Recover pressure
       MultiFab::Add (*p_g[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());

       // Exchange halo nodes and apply BCs
       mfix_set_velocity_bcs (lev,0);

       {
           amrex::Print() << "After projection in initial_iterations \n";
           mfix_print_max_vel (lev);
           mfix_compute_diveu (lev);
           amrex::Print() << "max(abs(diveu)) = " << diveu[lev] -> norm0 () << "\n";
       }

       // Replace vel_g by the original values 
       MultiFab::Copy (*vel_g[lev], *vel_go[lev],   0, 0, 3, vel_g[lev]->nGrow());
    }
}

//
// Compute predictor:
//
//  1. Compute
//
//     vel_g = vel_go + dt * R
//
//  2. Add explicit forcing term ( AKA gravity, lagged pressure gradient,
//     and explicit part of particles momentum exchange )
//
//     vel_g = vel_g + dt * ( g + grad(p_g+p0)/ro_g + drag_u/rop_g )
//
//  3. Add implicit forcing term ( AKA implicit part of particles
//     momentum exchange )
//
//     vel_g = vel_g / ( 1 + dt * f_gds/rop_g )
//
//  4. Solve
//
//     div( ep_g * grad(phi) / ro_g ) = div( ep_g * vel_g )
//
//  5. Compute
//
//     vel_g = vel_g -  dt * grad(phi) / ro_g
//
//  6. Compute
//
//     p_g = p_g + phi
//
void
mfix_level::mfix_apply_predictor (int lev, amrex::Real dt)
{
    // Compute fluid acceleration (convection + diffusion) 
    mfix_compute_velocity_slopes    ( lev, vel_go);
    mfix_compute_fluid_acceleration ( lev, vel_go);
    
    // First add the fluid acceleration
    MultiFab::Saxpy (*vel_g[lev], dt, *acc[lev], 0, 0, 3, 0);

    // Add the forcing terms
    mfix_apply_forcing_terms ( lev, dt, vel_g);
    mfix_add_pressure_gradient ( lev, -dt);

    // Compute intermediate velocity if drag terms present
    if (solve_dem)
       mfix_compute_intermediate_velocity ( lev, dt );
 
    // Exchange halo nodes and apply BCs to velocity
    mfix_set_velocity_bcs (lev,0);
 
    // Project velocity field
    mfix_apply_projection ( lev, dt );
 
    // Recover pressure
    MultiFab::Add (*p_g[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
 
    // Exchange halo nodes and apply BCs
    mfix_set_velocity_bcs (lev,0);
}

//
// Compute corrector:
//
//  1. Compute
//
//     vel_g = vel_go + dt * ((R_u^* + R_u^n) / 2 - gradp*(1/rho))
//
//     where the starred variables are computed using "predictor-step" variables.
//
//  2. Add explicit forcing term ( AKA gravity, lagged pressure gradient,
//     and explicit part of particles momentum exchange )
//
//     vel_g = vel_g + dt * ( g + grad(p_g+p0)/rp_g + drag_u/rop_g )
//
//  3. Add implicit forcing term ( AKA implicit part of particles
//     momentum exchange )
//
//     vel_g = vel_g / ( 1 + dt * f_gds/rop_g )
//
//  4. Solve
//
//     div( ep_g * grad(phi) / ro_g ) = div( ep_g * vel_g )
//
//  5. Compute
//
//     vel_g = vel_g -  dt * grad(phi) / ro_g
//
//  6. Compute
//
//     p_g = p_g + phi
//
void
mfix_level::mfix_apply_corrector (int lev, amrex::Real dt)
{
    BL_PROFILE("mfix_level::mfix_apply_corrector");

    // Compute fluid acceleration (convection + diffusion)
    // using first predictor
    mfix_compute_velocity_slopes    (lev, vel_g);
    mfix_compute_fluid_acceleration (lev, vel_g);
        
    // Store u_go + dt * R_u^* / 2
    MultiFab::LinComb ( *vel_g[lev], 1.0, *vel_go[lev], 0, dt/2.0, *acc[lev], 0, 0, 3, 0 ); 
	
    // Compute fluid acceleration (convection + diffusion) 
    // using velocity at the beginning of time step
    acc[lev]->setVal(0);
    mfix_compute_velocity_slopes    (lev, vel_go);
    mfix_compute_fluid_acceleration (lev, vel_go);
    
    // Add dt/2 * R_u^n 
    MultiFab::Saxpy (*vel_g[lev], dt/2.0, *acc[lev], 0, 0, 3, 0);

    // Add forcing terms
    mfix_apply_forcing_terms ( lev, dt, vel_g);
 
    // Add pressure gradient
    mfix_add_pressure_gradient ( lev, -dt );
 
    // Compute intermediate velocity if drag terms present
    if (solve_dem)
       mfix_compute_intermediate_velocity ( lev, dt );
 
    // Fill ghost cells and reimpose boundary conditions
    mfix_set_velocity_bcs (lev,0);
 
    // Apply projection
    mfix_apply_projection ( lev, dt );
 
    // Recover pressure
    MultiFab::Add (*p_g[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
 
    // Exchange halo nodes and apply BCs
    mfix_set_velocity_bcs (lev,0);
}

//
// Perform the following operations:
//
//       vel_g = vel_g + coeff * ( grad(p_g) ) * (1/ro_g)
//
void
mfix_level::mfix_add_pressure_gradient (int lev, amrex::Real coeff)
{
    BL_PROFILE("mfix_level::mfix_add_pressure_gradient");
    
    if (nodal_pressure == 1)
    {
#ifdef _OPENMP
#pragma omp parallel 
#endif
       for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
       {
   	   // Cell-centered tilebox
	   Box bx = mfi.tilebox();

	   add_grad_pnd (
	       BL_TO_FORTRAN_BOX(bx),  
	       BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
	       BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
               (*p_g[lev])[mfi].dataPtr(),
	       BL_TO_FORTRAN_ANYD((*p0_g[lev])[mfi]),
	       geom[lev].CellSize(), &coeff);
       }
    } else {
#ifdef _OPENMP
#pragma omp parallel 
#endif
       for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
       {
   	   // Tilebox
	   Box bx = mfi.tilebox();

	   add_grad_pcc (
	       BL_TO_FORTRAN_BOX(bx),  
	       BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
	       BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
               (*p_g[lev])[mfi].dataPtr(),
               (*p0_g[lev])[mfi].dataPtr(),
	       geom[lev].CellSize(), &coeff);
       }
    }
}

void
mfix_level::mfix_add_grad_phi (int lev, amrex::Real coeff)
{
    BL_PROFILE("mfix_level::mfix_add_grad_phi");
    
    if (nodal_pressure == 1)
    {
#ifdef _OPENMP
#pragma omp parallel 
#endif
       for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
       {
   	   // Tilebox
   	   Box bx = mfi.tilebox();

	   add_grad_phind (
	       BL_TO_FORTRAN_BOX(bx),  
	       BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
	       BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
	       BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
	       geom[lev].CellSize(), &coeff);
       }
    } else {
#ifdef _OPENMP
#pragma omp parallel 
#endif
       for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
       {
   	   // Tilebox
   	   Box bx = mfi.tilebox();

	   add_grad_phicc (
	       BL_TO_FORTRAN_BOX(bx),  
	       BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
	       BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
	       (*phi[lev])[mfi].dataPtr(),
	       geom[lev].CellSize(), &coeff);
       }
    }
    
}

//
// Compute acc using vel_g
//
void
mfix_level::mfix_compute_fluid_acceleration ( int lev,
					      Vector< std::unique_ptr<MultiFab> >& vel) 
{
    BL_PROFILE("mfix_level::mfix_compute_fluid_acceleration");
    Box domain(geom[lev].Domain());
    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*vel[lev],true); mfi.isValid(); ++mfi)
    {
	// Tilebox
	Box bx = mfi.tilebox ();

	compute_fluid_acceleration (
	    BL_TO_FORTRAN_BOX(bx),  
	    BL_TO_FORTRAN_ANYD((*acc[lev])[mfi]),
	    (*xslopes[lev])[mfi].dataPtr(),
	    (*yslopes[lev])[mfi].dataPtr(),
	    (*zslopes[lev])[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
            (*mu_g[lev])[mfi].dataPtr(),
            (*lambda_g[lev])[mfi].dataPtr(),
            BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
            domain.loVect (), domain.hiVect (),
	    bc_ilo.dataPtr(), bc_ihi.dataPtr(),
	    bc_jlo.dataPtr(), bc_jhi.dataPtr(),
	    bc_klo.dataPtr(), bc_khi.dataPtr(),
	    geom[lev].CellSize(), &nghost);
    }
}

void
mfix_level::mfix_apply_forcing_terms (int lev, amrex::Real dt,
				      Vector< std::unique_ptr<MultiFab> >& vel) 

{
    BL_PROFILE("mfix_level::mfix_apply_forcing_terms");

    Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
    {
	// Tilebox
	Box bx = mfi.tilebox ();

	add_forcing ( BL_TO_FORTRAN_BOX(bx),  
		      BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*drag[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
                      (*rop_g[lev])[mfi].dataPtr(),
		      domain.loVect (), domain.hiVect (),
		      geom[lev].CellSize (), &dt);
    }
}

//
// Implicit solve for the intermediate velocity.
// Currently this means accounting for the implicit part of the fluid/particle
// momentum exchange
// 
void
mfix_level::mfix_compute_intermediate_velocity ( int lev, amrex::Real dt )

{
    BL_PROFILE("mfix_level::mfix_compute_intermediate_velocity");

    // Whole domain
    Box domain(geom[lev].Domain());
    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi) {
	
	// Tilebox
	Box bx = mfi.tilebox();

	compute_intermediate_velocity ( BL_TO_FORTRAN_BOX(bx),  
					BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*f_gds[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
					&dt );
    }
}

//
// Compute the slopes of each velocity component in the
// three directions.
// 
void
mfix_level::mfix_compute_velocity_slopes (int lev,
                                          Vector< std::unique_ptr<MultiFab> >& vel)

{
    BL_PROFILE("mfix_level::mfix_compute_velocity_slopes");

    Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
    {
	// Tilebox
	Box bx = mfi.tilebox ();

	compute_slopes ( BL_TO_FORTRAN_BOX(bx),
		         BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
			 (*xslopes[lev])[mfi].dataPtr (),
			 (*yslopes[lev])[mfi].dataPtr (),
		         BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
			 domain.loVect (), domain.hiVect (),
			 bc_ilo.dataPtr(), bc_ihi.dataPtr(),
			 bc_jlo.dataPtr(), bc_jhi.dataPtr(),
			 bc_klo.dataPtr(), bc_khi.dataPtr(),
			 &nghost);

    }

    // Fill halo cells
    xslopes[lev] -> FillBoundary(geom[lev].periodicity());
    yslopes[lev] -> FillBoundary(geom[lev].periodicity());
    zslopes[lev] -> FillBoundary(geom[lev].periodicity());

}

//
// Computes the following decomposition:
//
//    u + grad(phi)/ro_g = u*,     where div(eps*u) = 0
//
// where u* is a non-div-free velocity field, stored
// by components in u_g, v_g, and w_g. The resulting div-free
// velocity field, u, overwrites the value of u* in u_g, v_g, and w_g.
//
// phi is an auxiliary function related to the pressure p_g by the relation:
//
//     new p_g  = old p_g + phi
void 
mfix_level::mfix_apply_projection ( int lev, amrex::Real scaling_factor )
{
    BL_PROFILE("mfix_level::mfix_apply_projection");

    // Compute right hand side, AKA div(ep_g* u)
    mfix_compute_diveu (lev);
    diveu[lev] -> mult ( 1.0/scaling_factor, diveu[lev]->nGrow() );

    // Compute the PPE coefficients
    mfix_compute_bcoeff ( lev );

    // Set BCs for Poisson's solver
    int bc_lo[3], bc_hi[3];
    int singular;
    Box domain(geom[lev].Domain());
    
    set_ppe_bc (bc_lo, bc_hi,
		domain.loVect(), domain.hiVect(),
		&nghost, 
		bc_ilo.dataPtr(), bc_ihi.dataPtr(),
		bc_jlo.dataPtr(), bc_jhi.dataPtr(),
		bc_klo.dataPtr(), bc_khi.dataPtr(),
		&singular );

    // Initialize phi to zero (any non-zero bc's are stored in p0)
    phi[lev] -> setVal(0.);
 
    // Solve PPE
    solve_poisson_equation ( lev, bcoeff, phi, diveu, bc_lo, bc_hi );
 
    // Correct the velocity field
    mfix_add_grad_phi ( lev, -scaling_factor );
}

//
// Solve PPE:
//
//                  div( eps_g/rho * grad(phi) ) = div(eps_g*u) 
// 
void
mfix_level::solve_poisson_equation (  int lev,
				      Vector< Vector< std::unique_ptr<MultiFab> > >& b,
				      Vector< std::unique_ptr<MultiFab> >& phi,
				      Vector< std::unique_ptr<MultiFab> >& rhs,
				      int bc_lo[], int bc_hi[] )
{
    BL_PROFILE("mfix_level::solve_poisson_equation");
    
    if (nodal_pressure == 1)
    {

       // 
       // First define the matrix (operator).
       //
       //        (del dot b sigma grad)) phi
       //
       LPInfo                       info;
       MLNodeLaplacian              matrix(geom, grids, dmap, info);

       matrix.setGaussSeidel(true);
       matrix.setHarmonicAverage(false);
       matrix.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
   			    {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );

       matrix.setSigma(0, *(b[lev][0]));
       matrix.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::Sigma);

       // By this point we must have filled the Dirichlet values of phi stored in the ghost cells
       phi[lev]->setVal(0.);
       matrix.setLevelBC ( lev, GetVecOfConstPtrs(phi)[lev] );

       // 
       // Then setup the solver ----------------------
       //
       MLMG  solver(matrix);
	
       solver.setMaxIter (mg_max_iter);
       solver.setMaxFmgIter (mg_max_fmg_iter);
       solver.setVerbose (mg_verbose);
       solver.setCGVerbose (mg_cg_verbose);
       solver.setCGMaxIter (mg_cg_maxiter);

       // 
       // Finally, solve the system
       //
       solver.solve ( GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol );

       phi[lev] -> FillBoundary (geom[lev].periodicity());

    } else {

       // 
       // First define the matrix (operator).
       // Class MLABecLaplacian describes the following operator:
       //
       //       (alpha * a - beta * (del dot b grad)) phi
       //
       LPInfo                       info;
       MLABecLaplacian              matrix(geom, grids, dmap, info);
       Vector<const MultiFab*>      tmp;
       array<MultiFab const*,AMREX_SPACEDIM>   b_tmp;

       // Copy the PPE coefficient into the proper data strutcure
       tmp = amrex::GetVecOfConstPtrs ( b[lev] ) ;
       b_tmp[0] = tmp[0];
       b_tmp[1] = tmp[1];
       b_tmp[2] = tmp[2];

       // It is essential that we set MaxOrder of the solver to 2
       // if we want to use the standard phi(i)-phi(i-1) approximation
       // for the gradient at Dirichlet boundaries.
       // The solver's default order is 3 and this uses three points for the
       // gradient at a Dirichlet boundary.
       matrix.setMaxOrder(2);

       // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
       matrix.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
   			    {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );
    
       matrix.setScalars ( 0.0, -1.0 );    
       matrix.setBCoeffs ( lev, b_tmp );

       // By this point we must have filled the Dirichlet values of phi stored in the ghost cells
       phi[lev]->setVal(0.);
       matrix.setLevelBC ( lev, GetVecOfConstPtrs(phi)[lev] );

       // 
       // Then setup the solver ----------------------
       //
       MLMG  solver(matrix);
	
       solver.setMaxIter (mg_max_iter);
       solver.setMaxFmgIter (mg_max_fmg_iter);
       solver.setVerbose (mg_verbose);
       solver.setCGVerbose (mg_cg_verbose);
       solver.setCGMaxIter (mg_cg_maxiter);

       // This ensures that ghost cells of phi are correctly filled when returned from the solver
       solver.setFinalFillBC(true);

       // 
       // Finally, solve the system
       //
       solver.solve ( GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol );

       phi[lev] -> FillBoundary (geom[lev].periodicity());
    }
}

//
// Compute div(ep_g * u)
// 
void
mfix_level::mfix_compute_diveu (int lev)
{
    Box domain(geom[lev].Domain());

    if (nodal_pressure == 1)
    {

       // Create a temporary multifab to hold (ep_g * vel_g)
       MultiFab vec(vel_g[lev]->boxArray(), vel_g[lev]->DistributionMap(), vel_g[lev]->nComp(), vel_g[lev]->nGrow());

       // Fill it with (ep_g * vel_g)
       vec.copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),vel_g[lev]->nGrow(),vel_g[lev]->nGrow());
       MultiFab::Multiply(vec,(*ep_g[lev]),0,0,1,1);

#ifdef _OPENMP
#pragma omp parallel
#endif
      // Extrapolate Dirichlet values to ghost cells -- but do it differently in that 
      //  no-slip walls are treated exactly like slip walls -- this is only relevant
      //  when going into the projection
      for (MFIter mfi(vec, true); mfi.isValid(); ++mfi)
        {
          set_vec_bcs ( BL_TO_FORTRAN_ANYD(vec[mfi]),
                        bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                        bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                        bc_klo.dataPtr(), bc_khi.dataPtr(),
                        domain.loVect(), domain.hiVect(),
                        &nghost);
        }

        vec.FillBoundary (geom[lev].periodicity());
    
#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIter mfi(*diveu[lev],true); mfi.isValid(); ++mfi)
       {
           // Note: this box is nodal!
	   const Box& bx = mfi.tilebox();
	
	   compute_diveund ( BL_TO_FORTRAN_BOX(bx),
			     BL_TO_FORTRAN_ANYD((*diveu[lev])[mfi]),
			     BL_TO_FORTRAN_ANYD(vec[mfi]),
	                     geom[lev].CellSize());
       }

    } else {

    int extrap_dir_bcs = 1;
    mfix_set_velocity_bcs (lev,extrap_dir_bcs);
    vel_g[lev]->FillBoundary (geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIter mfi(*diveu[lev],true); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.tilebox();
	
	   compute_diveucc ( BL_TO_FORTRAN_BOX(bx),
			     BL_TO_FORTRAN_ANYD((*diveu[lev])[mfi]),
  			     (*ep_g[lev])[mfi].dataPtr (),
  			     BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
	                     geom[lev].CellSize());
       }
    }

    // Restore velocities to carry Dirichlet values on faces
    int extrap_dir_bcs = 0;
    mfix_set_velocity_bcs (lev,extrap_dir_bcs);
}

//
// Computes bcoeff = ep_g/ro_g at the faces of the scalar cells
// 
void
mfix_level::mfix_compute_bcoeff (int lev)
{
    BL_PROFILE("mfix_level::mfix_compute_bcoeff");

    // Directions
    int xdir = 1;
    int ydir = 2;
    int zdir = 3;
    
    if (nodal_pressure == 1)
    {
#ifdef _OPENMP
#pragma omp parallel 
#endif
       for (MFIter mfi(*ro_g[lev],true); mfi.isValid(); ++mfi)
       {
	   // Cell-centered tilebox
	   Box bx = mfi.tilebox();

	   // X direction
	   compute_bcoeff_nd (BL_TO_FORTRAN_BOX(bx),
		              BL_TO_FORTRAN_ANYD((*(bcoeff[lev][0]))[mfi]),
		              BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		              (*ep_g[lev])[mfi].dataPtr(), &xdir );

	   // Y direction
	   compute_bcoeff_nd (BL_TO_FORTRAN_BOX(bx),
		              BL_TO_FORTRAN_ANYD((*(bcoeff[lev][1]))[mfi]),
		              BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		              (*ep_g[lev])[mfi].dataPtr(), &ydir );

	   // Z direction
	   compute_bcoeff_nd (BL_TO_FORTRAN_BOX(bx),
		              BL_TO_FORTRAN_ANYD((*(bcoeff[lev][2]))[mfi]),
		              BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		              (*ep_g[lev])[mfi].dataPtr(), &zdir );
	
       }
    } else {
#ifdef _OPENMP
#pragma omp parallel 
#endif
       for (MFIter mfi(*ro_g[lev],true); mfi.isValid(); ++mfi)
       {
	   // Tileboxes for staggered components
	   Box ubx = mfi.tilebox (e_x);
	   Box vbx = mfi.tilebox (e_y);
	   Box wbx = mfi.tilebox (e_z);

	   // X direction
	   compute_bcoeff_cc (BL_TO_FORTRAN_BOX(ubx),
		              BL_TO_FORTRAN_ANYD((*(bcoeff[lev][0]))[mfi]),
  		              BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		              (*ep_g[lev])[mfi].dataPtr(), &xdir );

	   // Y direction
	   compute_bcoeff_cc (BL_TO_FORTRAN_BOX(vbx),
		              BL_TO_FORTRAN_ANYD((*(bcoeff[lev][1]))[mfi]),
		              BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		              (*ep_g[lev])[mfi].dataPtr(), &ydir );

	   // Z direction
	   compute_bcoeff_cc (BL_TO_FORTRAN_BOX(wbx),
		              BL_TO_FORTRAN_ANYD((*(bcoeff[lev][2]))[mfi]),
		              BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		              (*ep_g[lev])[mfi].dataPtr(), &zdir );
	
       }
    }

    bcoeff[lev][0] -> FillBoundary(geom[lev].periodicity());
    bcoeff[lev][1] -> FillBoundary(geom[lev].periodicity());
    bcoeff[lev][2] -> FillBoundary(geom[lev].periodicity());
}

//
// Check if steady state has been reached by verifying that
// 
//      max(abs( u^(n+1) - u^(n) )) < tol * dt
//      max(abs( v^(n+1) - v^(n) )) < tol * dt
//      max(abs( w^(n+1) - w^(n) )) < tol * dt
// 

int
mfix_level::steady_state_reached (int lev, Real dt)
{
    //
    // Count number of access 
    //
    static int naccess = 0;

    //
    // Make sure velocity is up to date
    // 
    mfix_set_velocity_bcs (lev,0);
    
    // 
    // Use temporaries to store the difference
    // between current and previous solution
    // 
    MultiFab temp_vel(vel_g[lev]->boxArray(), vel_g[lev]->DistributionMap(),3,0);
    MultiFab::LinComb (temp_vel, 1.0, *vel_g[lev], 0, -1.0, *vel_go[lev], 0, 0, 3, 0);

    MultiFab tmp;
    
    if (nodal_pressure)
    {
       const BoxArray & nd_grid = amrex::convert(grids[lev], IntVect{1,1,1});
       tmp.define( nd_grid, dmap[lev], 1, 0 );     
    }
    else
    {
	tmp.define( grids[lev], dmap[lev], 1, 0 );
    }

    MultiFab::LinComb (tmp, 1.0, *p_g[lev], 0, -1.0, *p_go[lev], 0, 0, 1, 0);
    
    Real delta_u = temp_vel.norm0 (0);
    Real delta_v = temp_vel.norm0 (1);
    Real delta_w = temp_vel.norm0 (2);
    Real delta_p = tmp.norm0 ();
    
    Real tol = steady_state_tol; 
    
    int condition1 = (delta_u < tol*dt) && (delta_v < tol*dt ) && (delta_w < tol*dt);

    //
    // Second stop condition
    //
    Real du_n1 = temp_vel.norm1 (0, geom[lev].periodicity());
    Real dv_n1 = temp_vel.norm1 (1, geom[lev].periodicity());
    Real dw_n1 = temp_vel.norm1 (2, geom[lev].periodicity());
    Real dp_n1 = tmp.norm1 (0, geom[lev].periodicity());
    Real uo_n1 = vel_go[lev] -> norm1 (0, geom[lev].periodicity());
    Real vo_n1 = vel_go[lev] -> norm1 (1, geom[lev].periodicity());
    Real wo_n1 = vel_go[lev] -> norm1 (2, geom[lev].periodicity());
    Real po_n1 = p_go[lev] -> norm1 (0, geom[lev].periodicity());
    Real tmp1, tmp2, tmp3, tmp4;

    Real local_tol = 1.0e-8;
    
    if ( uo_n1 < local_tol ) {
    	tmp1 = 0.0;
    } else {
    	tmp1 = du_n1 / uo_n1;
    };

    if ( vo_n1 < local_tol ) {
    	tmp2 = 0.0;
    } else {
    	tmp2 = dv_n1 / vo_n1;
    };
    
    if ( wo_n1 < local_tol ) {
    	tmp3 = 0.0;
    } else {
    	tmp3 = dw_n1 / wo_n1;
    };

    if ( po_n1 < local_tol ) {
    	tmp4 = 0.0;
    } else {
    	tmp4 = dp_n1 / po_n1;
    };

    int condition2 = (tmp1 < tol) && (tmp2 < tol) && (tmp3 < tol); // && (tmp4 < tol);

    //
    // Print out info on steady state checks
    //
    amrex::Print() << "\nSteady state check:\n";
    amrex::Print() << "||u-uo||/||uo|| , du/dt  = " << tmp1 <<" , "<< delta_u/dt << "\n";
    amrex::Print() << "||v-vo||/||vo|| , dv/dt  = " << tmp2 <<" , "<< delta_v/dt << "\n";
    amrex::Print() << "||w-wo||/||wo|| , dw/dt  = " << tmp3 <<" , "<< delta_w/dt << "\n";
    amrex::Print() << "||p-po||/||po|| , dp/dt  = " << tmp4 <<" , "<< delta_p/dt << "\n";

    // Count # access
    naccess++;

    // 
    //  Always return negative to first access. This way
    //  initial zero velocity field do not test for false positive
    //
    if ( naccess == 1 ) {
	return 0;
    } else {
	return condition1 || condition2;
    };
}

//
// Set the BCs for all the variables EXCEPT pressure or velocity.
//
void
mfix_level::mfix_set_scalar_bcs (int lev)
{
  BL_PROFILE("mfix_level::mfix_set_scalar_bcs()");

  Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(*ep_g[lev], true); mfi.isValid(); ++mfi)
    {
      set_scalar_bcs ( BL_TO_FORTRAN_ANYD((*ep_g[lev])[mfi]),
                      (*ro_g[lev])[mfi].dataPtr (),
                      (*rop_g[lev])[mfi].dataPtr (),
                      (*mu_g[lev])[mfi].dataPtr (),
                      (*lambda_g[lev])[mfi].dataPtr (),
                      bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                      bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                      bc_klo.dataPtr(), bc_khi.dataPtr(),
                      domain.loVect(), domain.hiVect(),
                      &nghost );
    }
        ep_g[lev] -> FillBoundary (geom[lev].periodicity());
        ro_g[lev] -> FillBoundary (geom[lev].periodicity());
       rop_g[lev] -> FillBoundary (geom[lev].periodicity());
        mu_g[lev] -> FillBoundary (geom[lev].periodicity());
    lambda_g[lev] -> FillBoundary (geom[lev].periodicity());
}

//
// Set the BCs for velocity only
//
void
mfix_level::mfix_set_velocity_bcs (int lev, int extrap_dir_bcs)
{
  BL_PROFILE("mfix_level::mfix_set_velocity_bcs()");

  vel_g[lev] -> FillBoundary (geom[lev].periodicity());

  Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(*vel_g[lev], true); mfi.isValid(); ++mfi)
    {
      set_velocity_bcs ( BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
                         bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                         bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                         bc_klo.dataPtr(), bc_khi.dataPtr(),
                         domain.loVect(), domain.hiVect(),
                         &nghost, &extrap_dir_bcs );
    }
}

//
// Fills ghost cell values of pressure appropriately for the BC type
//
void
mfix_level::mfix_extrap_pressure (int lev, std::unique_ptr<amrex::MultiFab>& p)
{
    BL_PROFILE("mfix_level::mfix_extrap_pressure()");
    if (nodal_pressure == 1) return;
 
    Box domain(geom[lev].Domain());
 
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for (MFIter mfi(*p, true); mfi.isValid(); ++mfi) {
 
        extrap_pressure_to_ghost_cells (
            BL_TO_FORTRAN_ANYD((*p)[mfi]),
            bc_ilo.dataPtr(), bc_ihi.dataPtr(),
            bc_jlo.dataPtr(), bc_jhi.dataPtr(),
            bc_klo.dataPtr(), bc_khi.dataPtr(),
            domain.loVect(), domain.hiVect(),
            &nghost);
    }
}

void
mfix_level::check_for_nans (int lev)
{
    bool ug_has_nans = vel_g[lev] -> contains_nan (0);
    bool vg_has_nans = vel_g[lev] -> contains_nan (1);
    bool wg_has_nans = vel_g[lev] -> contains_nan (2);
    bool pg_has_nans =   p_g[lev] -> contains_nan (0);
    bool ropg_has_nans = rop_g[lev] -> contains_nan (0);

    if (ug_has_nans)
	std::cout << "WARNING: u_g contains NaNs!!!";

    if (vg_has_nans)
	std::cout << "WARNING: v_g contains NaNs!!!";

    if (wg_has_nans)
	std::cout << "WARNING: w_g contains NaNs!!!";

    if (pg_has_nans)
	std::cout << "WARNING: p_g contains NaNs!!!";

    if (ropg_has_nans)
	std::cout << "WARNING: rop_g contains NaNs!!!";

}

//
// Print the maximum values of the velocity components
//
void
mfix_level::mfix_print_max_vel(int lev)
{
    amrex::Print() << "max(abs(u/v/w/p))  = " << 
        vel_g[lev] -> norm0 (0) << "  " <<
        vel_g[lev] -> norm0 (1) << "  " <<
        vel_g[lev] -> norm0 (2) << "  " <<
          p_g[lev] -> norm0 () << "  " << std::endl;
}
