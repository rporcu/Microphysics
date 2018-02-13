#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>


void
mfix_level::EvolveFluidProjection(int lev, int nstep, int steady_state, Real& dt, Real& prev_dt, Real time, Real stop_time )
{
    BL_PROFILE_REGION_START("mfix::EvolveFluidProjection");

    amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";
    
    // Extrapolate boundary values for density and volume fraction
    // The subsequent call to mfix_set_projection_bcs will only overwrite
    // rop_g and ep_g ghost values for PINF and POUT
    fill_mf_bc ( lev, *rop_g[lev] );
    fill_mf_bc ( lev, *ep_g[lev] );
    fill_mf_bc ( lev, *mu_g[lev] );
    
    // Fill ghost nodes and reimpose boundary conditions
    mfix_set_projection_bcs (lev);

    //
    // Start loop: if we are not seeking a steady state solution,
    // the loop will execute only one time
    //
    int keep_looping = 1;
    int iter = 1;
    
    do
    {
	// Here we should check the CFL condition
	// Compute dt for this time step
	Real umax  = u_g[lev] -> norm0 ();
	Real vmax  = v_g[lev] -> norm0 ();
	Real wmax  = w_g[lev] -> norm0 ();
	Real romin = rop_g[lev] -> min (0);
	Real mumax = mu_g[lev] -> max (0);
	compute_new_dt ( &umax, &vmax, &wmax, &romin, &mumax,
			 geom[lev].CellSize(), &cfl, &time, &stop_time, &dt );
	prev_dt = dt ;

	if (steady_state)
	{
	   amrex::Print() << "\n   Iteration " << iter << " with dt = " << dt << "\n" << std::endl;
	} else {
	   amrex::Print() << "\n   Step " << nstep+1 <<": from old_time " \
	   	          << time << " to new_time " << time+dt 
                          << " with dt = " << dt << "\n" << std::endl;
	}

	// Back up field
	// Backup field variable to old
	int nghost = ep_go[lev] -> nGrow();
	MultiFab::Copy (*ep_go[lev],  *ep_g[lev],  0, 0, 1, nghost);
	MultiFab::Copy ( *p_go[lev],   *p_g[lev],  0, 0, 1, nghost);
	MultiFab::Copy (*ro_go[lev],  *ro_g[lev],  0, 0, 1, nghost);
	MultiFab::Copy (*rop_go[lev], *rop_g[lev], 0, 0, 1, nghost);
	MultiFab::Copy (*u_go[lev],   *u_g[lev],   0, 0, 1, nghost);
	MultiFab::Copy (*v_go[lev],   *v_g[lev],   0, 0, 1, nghost);
	MultiFab::Copy (*w_go[lev],   *w_g[lev],   0, 0, 1, nghost);
  
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
	amrex::Print() << "\nAfter predictor step:\n";
	
	mfix_print_max_vel (lev);

	mfix_compute_diveu (lev);
	
	amrex::Print() << "max(abs(diveu)) = " << diveu[lev] -> norm0 () << "\n";

	// Calculate drag coefficient
	if (solve_dem)
	    mfix_calc_drag_fluid(lev);

	// Corrector step 
	mfix_apply_corrector ( lev, dt );

	// Print info about correction step
	amrex::Print() << "\nAfter corrector step:\n";
	
	mfix_print_max_vel (lev);

	mfix_compute_diveu (lev);
	
	amrex::Print() << "max(abs(diveu)) = " << diveu[lev] -> norm0 () << "\n";
	    
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
	
	// Just for debugging
	// keep_looping = 0;
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
mfix_level::mfix_initial_iterations (int lev, Real stop_time)
{
    // Copy u_g into u_go
    MultiFab::Copy (*u_go[lev],   *u_g[lev],   0, 0, 1, u_go[lev]->nGrow());
    MultiFab::Copy (*v_go[lev],   *v_g[lev],   0, 0, 1, v_go[lev]->nGrow());
    MultiFab::Copy (*w_go[lev],   *w_g[lev],   0, 0, 1, w_go[lev]->nGrow());
  
    //  Here we should check the CFL condition
    // Compute dt for this time step
    Real umax  = u_g[lev] -> norm0 ();
    Real vmax  = v_g[lev] -> norm0 ();
    Real wmax  = w_g[lev] -> norm0 ();
    Real romin = rop_g[lev] -> min (0);
    Real mumax = mu_g[lev] -> max (0);
    
    Real time = 0.0;
    Real dt   = 1.e20;
    compute_new_dt ( &umax, &vmax, &wmax, &romin, &mumax,
		     geom[lev].CellSize(), &cfl, &time, &stop_time, &dt );

    // Calculate drag coefficient
    if (solve_dem)
        mfix_calc_drag_fluid(lev);

    // Compute fluid acceleration (convection + diffusion)
    mfix_compute_velocity_slopes ( lev, u_go, v_go, w_go );
    mfix_compute_fluid_acceleration ( lev, u_go, v_go, w_go );

    for (int iter = 0; iter < 3; ++iter)
    {
       // First add the fluid acceleration
       MultiFab::Saxpy (*u_g[lev], dt, *uacc[lev], 0, 0, 1, 0);
       MultiFab::Saxpy (*v_g[lev], dt, *vacc[lev], 0, 0, 1, 0);
       MultiFab::Saxpy (*w_g[lev], dt, *wacc[lev], 0, 0, 1, 0);

       // Add the forcing terms
       mfix_apply_forcing_terms ( lev, dt, u_g, v_g, w_g );
       mfix_add_pressure_gradient ( lev, -dt);

       // Compute intermediate velocity
       mfix_compute_intermediate_velocity ( lev, dt );

       // Exchange halo nodes and apply BCs to velocity
       mfix_set_projection_bcs (lev);

       // Project velocity field
       mfix_apply_projection ( lev, dt );

       // Recover pressure
       MultiFab::Add (*p_g[lev], *phi[lev], 0, 0, 1, 1);

       // Exchange halo nodes and apply BCs
       mfix_set_projection_bcs (lev);

       mfix_print_max_vel (lev);
       mfix_compute_diveu (lev);
       amrex::Print() << "max(abs(diveu)) = " << diveu[lev] -> norm0 () << "\n";

       // Replace u_g by the original values 
       MultiFab::Copy (*u_g[lev],   *u_go[lev],   0, 0, 1, u_g[lev]->nGrow());
       MultiFab::Copy (*v_g[lev],   *v_go[lev],   0, 0, 1, v_g[lev]->nGrow());
       MultiFab::Copy (*w_g[lev],   *w_go[lev],   0, 0, 1, w_g[lev]->nGrow());
    }
}

//
// Compute predictor.
//
// This routine solves:
//
//      du/dt  + grad(p)/rho  = RHS
//
// by using a second order discretization in space and
// a first order discretization in time + 
// non-incremental projection
//
//  1. Compute
// 
//     u_g = u_go + dt * R_u 
//     v_g = v_go + dt * R_v 
//     w_g = w_go + dt * R_w  
//
//  2. Solve
//
//     div( grad(phi) / rho ) = du_g/dx + dv_g/dy + dw_g/dz
//
//  3. Compute
//
//     u_g = u_g -  (dphi/dx) / rho 
//     v_g = v_g -  (dphi/dy) / rho 
//     w_g = w_g -  (dphi/dz) / rho
//
//  4. Compute
//
//     p_g = phi / dt
//
//
//  This is the predictor step of the Heun's integration
//  scheme, AKA Predictor-Corrector Method (PCM).
//  This step is first order in time and second order in space
// 
void
mfix_level::mfix_apply_predictor (int lev, amrex::Real dt)
{
    // Compute fluid acceleration (convection + diffusion) 
    mfix_compute_velocity_slopes ( lev, u_go, v_go, w_go );
    mfix_compute_fluid_acceleration ( lev, u_go, v_go, w_go );
    
    // First add the fluid acceleration
    MultiFab::Saxpy (*u_g[lev], dt, *uacc[lev], 0, 0, 1, 0);
    MultiFab::Saxpy (*v_g[lev], dt, *vacc[lev], 0, 0, 1, 0);
    MultiFab::Saxpy (*w_g[lev], dt, *wacc[lev], 0, 0, 1, 0);

    // Add the forcing terms
    mfix_apply_forcing_terms ( lev, dt, u_g, v_g, w_g );
    mfix_add_pressure_gradient ( lev, -dt);
    
    // Compute intermediate velocity
    mfix_compute_intermediate_velocity ( lev, dt );
    
    // Exchange halo nodes and apply BCs to velocity
    mfix_set_projection_bcs (lev);
 
    // Project velocity field
    mfix_apply_projection ( lev, dt );

    // Recover pressure
    MultiFab::Add (*p_g[lev], *phi[lev], 0, 0, 1, 1);
    
    // Exchange halo nodes and apply BCs
    mfix_set_projection_bcs (lev); 
}


//
// Compute corrector:
//
//  1. Compute
//
//     u_g = u_go + dt * (R_u^* + R_u^n - (dp*/dx)*(1/rho)) / 2
//     v_g = v_go + dt * (R_v^* + R_v^n - (dp*/dy)*(1/rho)) / 2
//     w_g = w_go + dt * (R_w^* + R_w^n - (dp*/dz)*(1/rho)) / 2
//
//     where the starred variables are the "predictor-step" variables. 
//     
//  2. Solve
//
//     div( grad(phi) / rho ) = du_g/dx + dv_g/dy + dw_g/dz
//
//  3. Compute
//
//     u_g = u_g - (dphi/dx) / rho 
//     v_g = v_g - (dphi/dy) / rho 
//     w_g = w_g - (dphi/dz) / rho
//
//  4. Compute
//
//     p_g = 2 * phi / dt
//
//
//  This is the correction step of the Heun's integration
//  scheme, AKA Predictor-Corrector Method (PCM).
//  This step is second order in time and space.
// 
void
mfix_level::mfix_apply_corrector (int lev, amrex::Real dt)
{
    BL_PROFILE("mfix_level::mfix_compute_second_predictor");

    // Compute fluid acceleration (convection + diffusion)
    // using first predictor
    mfix_compute_velocity_slopes ( lev, u_g, v_g, w_g );
    mfix_compute_fluid_acceleration ( lev, u_g, v_g, w_g );
        
    // Store u_go + dt * R_u^* / 2
    MultiFab::LinComb ( *u_g[lev], 1.0, *u_go[lev], 0, dt/2.0, *uacc[lev], 0, 0, 1, 0 ); 
    MultiFab::LinComb ( *v_g[lev], 1.0, *v_go[lev], 0, dt/2.0, *vacc[lev], 0, 0, 1, 0 );
    MultiFab::LinComb ( *w_g[lev], 1.0, *w_go[lev], 0, dt/2.0, *wacc[lev], 0, 0, 1, 0 );
	
    // Compute fluid acceleration (convection + diffusion) 
    // using velocity at the beginning of time step
    uacc[lev]->setVal(0.0);
    vacc[lev]->setVal(0.0);
    wacc[lev]->setVal(0.0);
    mfix_compute_velocity_slopes ( lev, u_go, v_go, w_go );
    mfix_compute_fluid_acceleration ( lev, u_go, v_go, w_go );
    
    // Add dt/2 * R_u^n 
    MultiFab::Saxpy (*u_g[lev], dt/2.0, *uacc[lev], 0, 0, 1, 0);
    MultiFab::Saxpy (*v_g[lev], dt/2.0, *vacc[lev], 0, 0, 1, 0);
    MultiFab::Saxpy (*w_g[lev], dt/2.0, *wacc[lev], 0, 0, 1, 0);

    // Add forcing terms
    mfix_apply_forcing_terms ( lev, dt, u_g, v_g, w_g );
    
    // Add pressure gradient
    mfix_add_pressure_gradient ( lev, -dt );

    // Compute intermediate velocity
    mfix_compute_intermediate_velocity ( lev, dt );
    
    // Fill ghost cells and reimpose boundary conditions
    mfix_set_projection_bcs (lev);

    // Apply projection
    mfix_apply_projection ( lev, dt );

    // Recover pressure
    MultiFab::Add (*p_g[lev], *phi[lev], 0, 0, 1, 1);
    
    // Exchange halo nodes and apply BCs
    mfix_set_projection_bcs (lev);
}


//
// Perform the following operations:
//
//       u_g = u_g + coeff * ( dp_g/dx ) * (1/ro_g)
//       v_g = v_g + coeff * ( dp_g/dy ) * (1/ro_g)
//       w_g = w_g + coeff * ( dp_g/dz ) * (1/ro_g)
//
void
mfix_level::mfix_add_pressure_gradient (int lev, amrex::Real coeff)
{
    BL_PROFILE("mfix_level::mfix_add_pressure_gradient");

    int xdir = 1;
    int ydir = 2;
    int zdir = 3;
    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi)
    {
	// Boxes for staggered components
	Box ubx = mfi.tilebox (e_x);
	Box vbx = mfi.tilebox (e_y);
	Box wbx = mfi.tilebox (e_z);

	add_grad_p (
	    BL_TO_FORTRAN_BOX(ubx),  
	    BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
	    (*p_g[lev])[mfi].dataPtr(),
	    (*p0_g[lev])[mfi].dataPtr(),
	    geom[lev].CellSize(), &coeff, &xdir );
	
	add_grad_p (
	    BL_TO_FORTRAN_BOX(vbx),  
	    BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
	    (*p_g[lev])[mfi].dataPtr(),
	    (*p0_g[lev])[mfi].dataPtr(),
	    geom[lev].CellSize(), &coeff, &ydir );
	

	add_grad_p (
	    BL_TO_FORTRAN_BOX(wbx),  
	    BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
	    (*p_g[lev])[mfi].dataPtr(),
	    (*p0_g[lev])[mfi].dataPtr(),
	    geom[lev].CellSize(), &coeff, &zdir );
    }
}

//
// Add gradient of phi to velocity field
//
void
mfix_level::mfix_add_grad_phi (int lev, amrex::Real coeff)
{
    BL_PROFILE("mfix_level::mfix_add_grad_phi");

    int xdir = 1;
    int ydir = 2;
    int zdir = 3;

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for (MFIter mfi(*phi[lev],true); mfi.isValid(); ++mfi)
    {
	// Boxes for staggered components
	Box ubx = mfi.tilebox (e_x);
	Box vbx = mfi.tilebox (e_y);
	Box wbx = mfi.tilebox (e_z);
	
	add_grad_phi (
	    BL_TO_FORTRAN_BOX(ubx),
	    BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
	    (*ro_g[lev])[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
	    geom[lev].CellSize(), &coeff, &xdir );

	add_grad_phi (
	    BL_TO_FORTRAN_BOX(vbx),
	    BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
	    (*ro_g[lev])[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
	    geom[lev].CellSize(), &coeff, &ydir );

	add_grad_phi (
	    BL_TO_FORTRAN_BOX(wbx),
	    BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
	    (*ro_g[lev])[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
	    geom[lev].CellSize(), &coeff, &zdir );
    }
    
}


//
// Compute uacc, vacc, and wacc by using u_g, v_g, and w_g
//
void
mfix_level::mfix_compute_fluid_acceleration ( int lev,
					      Vector< std::unique_ptr<MultiFab> >& u, 
					      Vector< std::unique_ptr<MultiFab> >& v,
					      Vector< std::unique_ptr<MultiFab> >& w )
{
    BL_PROFILE("mfix_level::mfix_compute_fluid_acceleration");

    int xdir = 1;
    int ydir = 2;
    int zdir = 3;
    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi)
    {
	// Boxes for staggered components
	Box ubx = mfi.tilebox (e_x);
	Box vbx = mfi.tilebox (e_y);
	Box wbx = mfi.tilebox (e_z);

	// x direction
	compute_fluid_acceleration (
	    BL_TO_FORTRAN_BOX(ubx),  
	    BL_TO_FORTRAN_ANYD((*uacc[lev])[mfi]),
	    (*slopes_u[lev])[mfi].dataPtr (),
	    BL_TO_FORTRAN_ANYD((*u[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*rop_g[lev])[mfi].dataPtr (),
	    geom[lev].CellSize (), &xdir );

	// y direction
	compute_fluid_acceleration (
	    BL_TO_FORTRAN_BOX(vbx),  
	    BL_TO_FORTRAN_ANYD((*vacc[lev])[mfi]),
	    (*slopes_v[lev])[mfi].dataPtr (),
	    BL_TO_FORTRAN_ANYD((*u[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*rop_g[lev])[mfi].dataPtr (),
	    geom[lev].CellSize (), &ydir );

	// z direction
	compute_fluid_acceleration (
	    BL_TO_FORTRAN_BOX(wbx),  
	    BL_TO_FORTRAN_ANYD((*wacc[lev])[mfi]),
	    (*slopes_w[lev])[mfi].dataPtr (),
	    BL_TO_FORTRAN_ANYD((*u[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*rop_g[lev])[mfi].dataPtr (),
	    geom[lev].CellSize (), &zdir );
    }
}


//
// Add explicit forcing terms. These include the body forces and
// the explicit part of the particle/fluid momentum exchange 
// 
void
mfix_level::mfix_apply_forcing_terms (int lev, amrex::Real dt,
				      Vector< std::unique_ptr<MultiFab> >& u, 
				      Vector< std::unique_ptr<MultiFab> >& v,
				      Vector< std::unique_ptr<MultiFab> >& w )

{
    BL_PROFILE("mfix_level::mfix_apply_forcing_terms");

    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi) {
	
	// Directions
	int xdir = 1;
	int ydir = 2;
	int zdir = 3;
	
	// Whole domain
	Box domain(geom[lev].Domain());
	
	// Boxes for staggered components
	Box ubx = mfi.tilebox (e_x);
	Box vbx = mfi.tilebox (e_y);
	Box wbx = mfi.tilebox (e_z);
	

	add_forcing ( BL_TO_FORTRAN_BOX(ubx),  
		      BL_TO_FORTRAN_ANYD((*u[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*drag_u[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		      (*rop_g[lev])[mfi].dataPtr(),
		      domain.loVect (), domain.hiVect (),
		      geom[lev].CellSize (), &dt, &xdir );

	add_forcing ( BL_TO_FORTRAN_BOX(vbx),  
		      BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*drag_v[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		      (*rop_g[lev])[mfi].dataPtr(),
		      domain.loVect (), domain.hiVect (),
		      geom[lev].CellSize (), &dt, &ydir );

	add_forcing ( BL_TO_FORTRAN_BOX(wbx),  
		      BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*drag_w[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		      (*rop_g[lev])[mfi].dataPtr(),
		      domain.loVect (), domain.hiVect (),
		      geom[lev].CellSize (), &dt, &zdir );
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
    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi) {
	
	// Directions
	int xdir = 1;
	int ydir = 2;
	int zdir = 3;
	
	// Whole domain
	Box domain(geom[lev].Domain());
	
	// Boxes for staggered components
	Box ubx = mfi.tilebox (e_x);
	Box vbx = mfi.tilebox (e_y);
	Box wbx = mfi.tilebox (e_z);
	
	compute_intermediate_velocity ( BL_TO_FORTRAN_BOX(ubx),  
					BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*f_gds_u[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
					&xdir, &dt );

	compute_intermediate_velocity ( BL_TO_FORTRAN_BOX(vbx),  
					BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*f_gds_v[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
					&ydir, &dt );
	
	compute_intermediate_velocity ( BL_TO_FORTRAN_BOX(wbx),  
					BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*f_gds_w[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
					&zdir, &dt );
    }

}




//
// Compute the slopes of each velocity component in the
// three directions.
// 
void
mfix_level::mfix_compute_velocity_slopes (int lev,
                                          Vector< std::unique_ptr<MultiFab> >& u,
                                          Vector< std::unique_ptr<MultiFab> >& v,
                                          Vector< std::unique_ptr<MultiFab> >& w )

{
    BL_PROFILE("mfix_level::mfix_compute_velocity_slopes");

#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi)
    {
	// Boxes for staggered components
	Box domain(geom[lev].Domain());
	Box ubx = mfi.tilebox (e_x);
	Box vbx = mfi.tilebox (e_y);
	Box wbx = mfi.tilebox (e_z);

	compute_u_slopes ( BL_TO_FORTRAN_BOX(ubx),
			   BL_TO_FORTRAN_ANYD((*u[lev])[mfi]),
			   (*slopes_u[lev])[mfi].dataPtr (),
			   domain.loVect (), domain.hiVect (),
			   bc_ilo.dataPtr(), bc_ihi.dataPtr() );

	compute_v_slopes ( BL_TO_FORTRAN_BOX(vbx),
			   BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
			   (*slopes_v[lev])[mfi].dataPtr (),
			   domain.loVect (), domain.hiVect (),
			   bc_jlo.dataPtr(), bc_jhi.dataPtr() );

	compute_w_slopes ( BL_TO_FORTRAN_BOX(wbx),
			   BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
			   (*slopes_w[lev])[mfi].dataPtr (),
			   domain.loVect (), domain.hiVect (),
			   bc_klo.dataPtr(), bc_khi.dataPtr() );
    }

    // Fill halo cells
    slopes_u[lev] -> FillBoundary(geom[lev].periodicity());
    slopes_v[lev] -> FillBoundary(geom[lev].periodicity());
    slopes_w[lev] -> FillBoundary(geom[lev].periodicity());
}


//
// Computes the following decomposition:
// 
//    u + grad(phi)/rho = u*,     where div(eps*u) = 0
//
// where u* is a non-div-free velocity field, stored
// by components in u_g, v_g, and w_g. The resulting div-free
// velocity field, u, overwrites the value of u* in u_g, v_g, and w_g.
//
// phi is an auxiliary function related to the pressure p_g by the relation:
// 
//     phi = scaling_factor * p_g
//
// p_g is not required to have any value set, except at the Dirichlet's boundary.
// 
void 
mfix_level::mfix_apply_projection ( int lev, amrex::Real scaling_factor )
{
    BL_PROFILE("mfix_level::mfix_apply_projection");

    // Compute right hand side, AKA div(ep_g* u)
    mfix_compute_diveu (lev);
    diveu[lev] -> mult ( 1.0/scaling_factor, 1 );
    
    // Compute the PPE coefficients
    mfix_compute_bcoeff ( lev );

    // Set BCs for Poisson's solver
    int bc_lo[3], bc_hi[3];
    int singular;
    Box domain(geom[lev].Domain());
    
    set_ppe_bc (bc_lo, bc_hi,
		domain.loVect(), domain.hiVect(),
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
    
    if (singular) {
	Real phi_mean = ( phi[lev] -> sum () ) / domain.numPts () ;
	phi[lev] -> plus ( -phi_mean, 0 ); // pg_mean is 0 for non-singular case
    }
}

//
// Solve PPE:
//
//                  div( b * grad(phi) ) = div(u) 
// 
void
mfix_level::solve_poisson_equation (  int lev,
				      Vector< Vector< std::unique_ptr<MultiFab> > >& b,
				      Vector< std::unique_ptr<MultiFab> >& phi,
				      Vector< std::unique_ptr<MultiFab> >& rhs,
				      int bc_lo[], int bc_hi[] )
{
    BL_PROFILE("mfix_level::solve_poisson_equation");
    
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
			 {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]});

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
    // solver.setFinalFillBC(true);

    // 
    // Finally, solve the system
    //
    solver.solve ( GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol );

    phi[lev] -> FillBoundary (geom[lev].periodicity());

}



//
// Compute div(ep_g * {u,v,w})
// 
void
mfix_level::mfix_compute_diveu (int lev)
{

    u_g[lev] -> FillBoundary (geom[lev].periodicity());
    v_g[lev] -> FillBoundary (geom[lev].periodicity());
    w_g[lev] -> FillBoundary (geom[lev].periodicity());
    mfix_set_projection_bcs (lev);

    fill_mf_bc (lev,*ep_g[lev]);


    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*diveu[lev],true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();
	
	compute_diveu ( BL_TO_FORTRAN_BOX(bx),
			BL_TO_FORTRAN_ANYD((*diveu[lev])[mfi]),
			(*ep_g[lev])[mfi].dataPtr (),
			BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
			BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
			BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
			geom[lev].CellSize() );		     
	
    }

    // BCs 
    fill_mf_bc (lev,*diveu[lev]);

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
    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi)
    {
	// Boxes for staggered components
	Box ubx = mfi.tilebox (e_x);
	Box vbx = mfi.tilebox (e_y);
	Box wbx = mfi.tilebox (e_z);

	// X direction
	compute_bcoeff (BL_TO_FORTRAN_BOX(ubx),
		       BL_TO_FORTRAN_ANYD((*(bcoeff[lev][0]))[mfi]),
		       BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		       (*ep_g[lev])[mfi].dataPtr(), &xdir );

	// Y direction
	compute_bcoeff (BL_TO_FORTRAN_BOX(vbx),
		       BL_TO_FORTRAN_ANYD((*(bcoeff[lev][1]))[mfi]),
		       BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		       (*ep_g[lev])[mfi].dataPtr(), &ydir );

	// Z direction
	compute_bcoeff (BL_TO_FORTRAN_BOX(wbx),
		       BL_TO_FORTRAN_ANYD((*(bcoeff[lev][2]))[mfi]),
		       BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		       (*ep_g[lev])[mfi].dataPtr(), &zdir );
	
	
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
    // Make sure all ghost nodes are up to date
    // 
    mfix_set_projection_bcs (lev);
    
    // 
    // Use temporaries to store the difference
    // between current and previous solution
    // 
    MultiFab::LinComb (*u_gt[lev], 1.0, *u_g[lev], 0, -1.0, *u_go[lev], 0, 0, 1, 0);
    MultiFab::LinComb (*v_gt[lev], 1.0, *v_g[lev], 0, -1.0, *v_go[lev], 0, 0, 1, 0);
    MultiFab::LinComb (*w_gt[lev], 1.0, *w_g[lev], 0, -1.0, *w_go[lev], 0, 0, 1, 0);

    MultiFab tmp( grids[lev], dmap[lev], 1, 0 );
    MultiFab::LinComb (tmp, 1.0, *p_g[lev], 0, -1.0, *p_go[lev], 0, 0, 1, 0);
    
    Real delta_u = u_gt[lev] -> norm0 ();
    Real delta_v = v_gt[lev] -> norm0 ();
    Real delta_w = w_gt[lev] -> norm0 ();
    Real delta_p = tmp.norm0 ();
    
    Real tol = steady_state_tol; 

    int condition1 = (delta_u < tol*dt) && (delta_v < tol*dt ) && (delta_w < tol*dt);


    //
    // Second stop condition
    //
    Real du_n1 = u_gt[lev] -> norm1 (0, geom[lev].periodicity());
    Real dv_n1 = v_gt[lev] -> norm1 (0, geom[lev].periodicity());
    Real dw_n1 = w_gt[lev] -> norm1 (0, geom[lev].periodicity());
    Real dp_n1 = tmp.norm1 (0, geom[lev].periodicity());
    Real uo_n1 = u_go[lev] -> norm1 (0, geom[lev].periodicity());
    Real vo_n1 = v_go[lev] -> norm1 (0, geom[lev].periodicity());
    Real wo_n1 = w_go[lev] -> norm1 (0, geom[lev].periodicity());
    Real po_n1 = p_go[lev] -> norm1 (0, geom[lev].periodicity());
    Real tmp1, tmp2, tmp3, tmp4;
    
    if ( uo_n1 < 1.0e-8 ) {
    	tmp1 = 0.0;
    } else {
    	tmp1 = du_n1 / uo_n1;
    };

    if ( vo_n1 < 1.0e-8 ) {
    	tmp2 = 0.0;
    } else {
    	tmp2 = dv_n1 / vo_n1;
    };
    
    if ( wo_n1 < 1.0e-8 ) {
    	tmp3 = 0.0;
    } else {
    	tmp3 = dw_n1 / wo_n1;
    };

    if ( po_n1 < 1.0e-8 ) {
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
// Set the BCs for all the variables EXCEPT pressure
// Call set_pressure_bcs to impose BCs on pressure
// field
// 
void
mfix_level::mfix_set_projection_bcs (int lev)
{
  BL_PROFILE("mfix_level::mfix_set_projection_bcs()");

  u_g[lev] -> FillBoundary (geom[lev].periodicity());
  v_g[lev] -> FillBoundary (geom[lev].periodicity());
  w_g[lev] -> FillBoundary (geom[lev].periodicity());

  
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(*p_g[lev], true); mfi.isValid(); ++mfi)
    {
      Box domain(geom[lev].Domain());
      const Box& sbx = (*p_g[lev])[mfi].box();
      Box ubx((*u_g[lev])[mfi].box());
      Box vbx((*v_g[lev])[mfi].box());
      Box wbx((*w_g[lev])[mfi].box());


      set_projection_bcs ( BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*ep_g[lev])[mfi]),
			   (*ro_g[lev])[mfi].dataPtr (),
			   (*rop_g[lev])[mfi].dataPtr (),
			   (*mu_g[lev])[mfi].dataPtr (),
			   (*lambda_g[lev])[mfi].dataPtr (),
			   bc_ilo.dataPtr(), bc_ihi.dataPtr(),
			   bc_jlo.dataPtr(), bc_jhi.dataPtr(),
			   bc_klo.dataPtr(), bc_khi.dataPtr(),
			   domain.loVect(), domain.hiVect());
    }			   

}

//
// Fills ghost cell values of pressure appropriately for the BC type
//
void
mfix_level::mfix_extrap_pressure (int lev, std::unique_ptr<amrex::MultiFab>& p)
{
    BL_PROFILE("mfix_level::mfix_extrap_pressure()");

    Box domain(geom[lev].Domain());

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for (MFIter mfi(*p, true); mfi.isValid(); ++mfi) {

	const Box& sbx = (*p)[mfi].box();

	extrap_pressure_to_ghost_cells (
	    BL_TO_FORTRAN_ANYD((*p)[mfi]),
	    bc_ilo.dataPtr(), bc_ihi.dataPtr(),
	    bc_jlo.dataPtr(), bc_jhi.dataPtr(),
	    bc_klo.dataPtr(), bc_khi.dataPtr(),
	    domain.loVect(), domain.hiVect());
    }
}


void
mfix_level::check_for_nans (int lev)
{
    bool ug_has_nans = u_g[lev] -> contains_nan ();
    bool vg_has_nans = v_g[lev] -> contains_nan ();
    bool wg_has_nans = w_g[lev] -> contains_nan ();
    bool pg_has_nans = p_g[lev] -> contains_nan ();
    bool ropg_has_nans = rop_g[lev] -> contains_nan ();

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
    amrex::Print() << "max(abs(u/v/w/p))  = " << u_g[lev] -> norm0 () << "  " <<
	v_g[lev] -> norm0 () << "  " <<
	w_g[lev] -> norm0 () << "  " <<
	p_g[lev] -> norm0 () << "  " << std::endl;
}
