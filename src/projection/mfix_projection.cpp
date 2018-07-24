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
mfix_level::EvolveFluidProjection(int lev, int nstep, int steady_state, Real& dt, Real time, Real stop_time )
{
    BL_PROFILE_REGION_START("mfix::EvolveFluidProjection");
    BL_PROFILE("mfix::EvolveFluidProjection");
    
    amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";
    
    // Extrapolate boundary values for density and volume fraction
    // The subsequent call to mfix_set_projection_bcs will only overwrite
    // ep_g ghost values for PINF and POUT
    fill_mf_bc ( lev, *ep_g[lev] );
    fill_mf_bc ( lev, *mu_g[lev] );
    
    // Fill ghost nodes and reimpose boundary conditions
    mfix_set_projection_bcs (lev);
    mfix_set_velocity_bcs (lev);

    //
    // Start loop: if we are not seeking a steady state solution,
    // the loop will execute only once
    //
    int keep_looping = 1;
    int iter = 1;
    
    do
    {
        mfix_compute_dt(lev,time,stop_time,steady_state,dt);

	if (steady_state)
	{
	   amrex::Print() << "\n   Iteration " << iter << " with dt = " << dt << "\n" << std::endl;
	} else {
	   amrex::Print() << "\n   Step " << nstep+1 <<": from old time " \
	   	          << time << " to new time " << time+dt 
                          << " with dt = " << dt << "\n" << std::endl;
	}

	// Back up field
	// Backup field variable to old
	MultiFab::Copy(*ep_go[lev],  *ep_g[lev],  0, 0, 1,  ep_go[lev]->nGrow() );
	MultiFab::Copy( *p_go[lev],   *p_g[lev],  0, 0, 1,   p_go[lev]->nGrow() );
	MultiFab::Copy(*ro_go[lev],  *ro_g[lev],  0, 0, 1,  ro_go[lev]->nGrow() );
	MultiFab::Copy(*u_go[lev],   *u_g[lev],   0, 0, 1,   u_go[lev]->nGrow() );
	MultiFab::Copy(*v_go[lev],   *v_g[lev],   0, 0, 1,   v_go[lev]->nGrow() );
	MultiFab::Copy(*w_go[lev],   *w_g[lev],   0, 0, 1,   w_go[lev]->nGrow() );

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


	// Calculate drag coefficient
	if (solve_dem)
	    mfix_calc_drag_fluid(lev);

	// Corrector step 
	mfix_apply_corrector ( lev, dt );

	// Print info about correction step
	amrex::Print() << "\nAfter corrector step:\n";
	
	mfix_print_max_vel (lev);

	    
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

//
// Project velocity field to make sure initial velocity is divergence-free
// 
void
mfix_level::mfix_project_velocity (int lev)
{
    Real dummy_dt = 1.0;

    mac_projection -> apply_projection (u_g, v_g, w_g, ep_g, ro_g, dummy_dt );

    amrex::Print() << "After initial projection:\n";

    mfix_print_max_vel (lev);
	   
}


//
// Perform 3 iterations to compute an initial guess for p at time t=1/2
// Each iteration is an Euler time integration step, except that velocity
// is never updated. Basically, what we are doing is:
//
//  1)
//    u0, p0, rop0, R0 (RHS in momentum eq.) are initialized with
//    initial conditions and never modified for the whole iteration process,
//
//  2)
//     perform  the following iteration (only pressure get updated):
// 
//     do iter = 1, 3
//        u* = u0 + dt * R0 - grad(p^{k-1})/rop0
//        solve grad(eps_g/ro_g* grad(phi)) = div(u*)/dt
//        p^{k} = p^{k-1} + phi
//     end do
//
//  3) p at time t=1/2 is set to p^3 
//         
// 
void
mfix_level::mfix_initial_iterations (int lev, Real dt, Real stop_time, int steady_state)
{
    // Copy u_g into u_go
    MultiFab::Copy (*u_go[lev],   *u_g[lev],   0, 0, 1, u_go[lev]->nGrow());
    MultiFab::Copy (*v_go[lev],   *v_g[lev],   0, 0, 1, v_go[lev]->nGrow());
    MultiFab::Copy (*w_go[lev],   *w_g[lev],   0, 0, 1, w_go[lev]->nGrow());
    MultiFab::Copy (*ro_go[lev],  *ro_g[lev],   0, 0, 1, ro_go[lev]->nGrow());
    MultiFab::Copy (*ep_go[lev],  *ep_g[lev],   0, 0, 1, ep_go[lev]->nGrow());

    Real time = 0.0;
    mfix_compute_dt(lev,time,stop_time,steady_state,dt);

    amrex::Print() << "Doing initial pressure iterations with dt = " << dt << std::endl;

    // Calculate drag coefficient
    if (solve_dem)
        mfix_calc_drag_fluid(lev);

    // Compute fluid acceleration (convection + diffusion)
    mfix_compute_velocity_slopes ( lev, u_go, v_go, w_go );
    mfix_compute_fluid_acceleration ( lev, u_go, v_go, w_go, ro_go, ep_go );

    for (int iter = 0; iter < 3; ++iter)
    {
       // First add the fluid acceleration
       MultiFab::Saxpy (*u_g[lev], dt, *uacc[lev], 0, 0, 1, 0);
       MultiFab::Saxpy (*v_g[lev], dt, *vacc[lev], 0, 0, 1, 0);
       MultiFab::Saxpy (*w_g[lev], dt, *wacc[lev], 0, 0, 1, 0);

       // Add the forcing terms
       mfix_apply_forcing_terms ( lev, dt, u_g, v_g, w_g, ro_go, ep_go );

       mfix_add_pressure_gradient ( lev, -dt);

       // Compute intermediate velocity
       mfix_compute_intermediate_velocity ( lev, dt, ro_go, ep_go );

       // Projection
       mac_projection -> apply_projection (u_g, v_g, w_g, ep_g, ro_g, dt );

       // Recover pressure
       const MultiFab& phi_mac = mac_projection -> get_phi(lev);
       MultiFab::Add (*p_g[lev], phi_mac, 0, 0, 1, p_g[lev] -> nGrow() );

       mfix_print_max_vel (lev);

       // Replace u_g by the original values 
       MultiFab::Copy (*u_g[lev],   *u_go[lev],   0, 0, 1, u_g[lev]->nGrow());
       MultiFab::Copy (*v_g[lev],   *v_go[lev],   0, 0, 1, v_g[lev]->nGrow());
       MultiFab::Copy (*w_g[lev],   *w_go[lev],   0, 0, 1, w_g[lev]->nGrow());
    }

}

//
// Compute predictor:
//
//  1. Compute
// 
//     u_g = u_go + dt * R_u 
//     v_g = v_go + dt * R_v 
//     w_g = w_go + dt * R_w  
//
//
//  2. Add explicit forcing term ( AKA gravity, lagged pressure gradient,
//     and explicit part of particles momentum exchange )
// 
//     u_g = u_g + dt * ( g_x + d(p_g+p0)/dx/ro_g + drag_u/( ro_go + ep_go ) ) 
//     v_g = v_g + dt * ( g_y + d(p_g+p0)/dy/ro_g + drag_v/( ro_go + ep_go ) ) 
//     w_g = w_g + dt * ( g_z + d(p_g+p0)/dz/ro_g + drag_w/( ro_go + ep_go ) ) 
//
//  3. Add implicit forcing term ( AKA implicit part of particles
//     momentum exchange )
// 
//     u_g = u_g / ( 1 + dt * f_gds_u/(ro_go + ep_go) )
//     v_g = v_g / ( 1 + dt * f_gds_v/(ro_go + ep_go) )
//     w_g = w_g / ( 1 + dt * f_gds_w/(ro_go + ep_go) )
//  
//  4. Solve
//
//     div( ep_g * grad(phi) / ro_g ) = div( ep_g * {u_g,v_g,w_g} ) 
//
//  5. Compute
//
//     u_g = u_g -  dt * (dphi/dx) / ro_g 
//     v_g = v_g -  dt * (dphi/dy) / ro_g 
//     w_g = w_g -  dt * (dphi/dz) / ro_g
//
//  6. Compute
//
//     p_g = p_g + phi
//
void
mfix_level::mfix_apply_predictor (int lev, amrex::Real dt)
{
    // Compute fluid acceleration (convection + diffusion) 
    mfix_compute_velocity_slopes ( lev, u_go, v_go, w_go );
    mfix_compute_fluid_acceleration ( lev, u_go, v_go, w_go, ro_go, ep_go );
    
    // First add the fluid acceleration
    MultiFab::Saxpy (*u_g[lev], dt, *uacc[lev], 0, 0, 1, 0);
    MultiFab::Saxpy (*v_g[lev], dt, *vacc[lev], 0, 0, 1, 0);
    MultiFab::Saxpy (*w_g[lev], dt, *wacc[lev], 0, 0, 1, 0);

    // Add the forcing terms
    mfix_apply_forcing_terms ( lev, dt, u_g, v_g, w_g, ro_go, ep_go );
    mfix_add_pressure_gradient ( lev, -dt);
    
    // Compute intermediate velocity
    mfix_compute_intermediate_velocity ( lev, dt, ro_go, ep_go );

    // Projection
    mac_projection -> apply_projection (u_g, v_g, w_g, ep_g, ro_g, dt );

    // Recover pressure
    const MultiFab& phi_mac = mac_projection -> get_phi(lev);
    MultiFab::Add (*p_g[lev], phi_mac, 0, 0, 1, p_g[lev] -> nGrow() );
}


//
// Compute corrector:
//
//  1. Compute
// 
//     u_g = u_go + dt * ( R_u^n + R_u^* ) / 2
//     v_g = v_go + dt * ( R_v^n + R_v^* ) / 2
//     w_g = w_go + dt * ( R_w^n + R_w^* ) / 2
//
//     where the starred variables are computed using "predictor-step" variables. 
//
//  2. Add explicit forcing term ( AKA gravity, lagged pressure gradient,
//     and explicit part of particles momentum exchange )
// 
//     u_g = u_g + dt * ( g_x + d(p_g+p0)/dx/ro_g + drag_u/(ro_g + ep_g) ) 
//     v_g = v_g + dt * ( g_y + d(p_g+p0)/dy/ro_g + drag_v/(ro_g + ep_g) ) 
//     w_g = w_g + dt * ( g_z + d(p_g+p0)/dz/ro_g + drag_w/(ro_g + ep_g) ) 
//
//  3. Add implicit forcing term ( AKA implicit part of particles
//     momentum exchange )
// 
//     u_g = u_g / ( 1 + dt * f_gds_u/(ro_g + ep_g) )
//     v_g = v_g / ( 1 + dt * f_gds_v/(ro_g + ep_g) )
//     w_g = w_g / ( 1 + dt * f_gds_w/(ro_g + ep_g) )
//  
//  4. Solve
//
//     div( ep_g * grad(phi) / ro_g ) = div( ep_g * {u_g,v_g,w_g} ) 
//
//  5. Compute
//
//     u_g = u_g -  dt * (dphi/dx) / ro_g 
//     v_g = v_g -  dt * (dphi/dy) / ro_g 
//     w_g = w_g -  dt * (dphi/dz) / ro_g
//
//  6. Compute
//
//     p_g = p_g + phi
//
void
mfix_level::mfix_apply_corrector (int lev, amrex::Real dt)
{
    BL_PROFILE("mfix_level::mfix_compute_second_predictor");

    // Compute fluid acceleration (convection + diffusion)
    // using first predictor
    mfix_compute_velocity_slopes ( lev, u_g, v_g, w_g );
    mfix_compute_fluid_acceleration ( lev, u_g, v_g, w_g, ro_g, ep_g );
        
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
    mfix_compute_fluid_acceleration ( lev, u_go, v_go, w_go, ro_go, ep_go );
    
    // Add dt/2 * R_u^n 
    MultiFab::Saxpy (*u_g[lev], dt/2.0, *uacc[lev], 0, 0, 1, 0);
    MultiFab::Saxpy (*v_g[lev], dt/2.0, *vacc[lev], 0, 0, 1, 0);
    MultiFab::Saxpy (*w_g[lev], dt/2.0, *wacc[lev], 0, 0, 1, 0);

    // Add forcing terms
    mfix_apply_forcing_terms ( lev, dt, u_g, v_g, w_g, ro_g, ep_g );
    
    // Add pressure gradient
    mfix_add_pressure_gradient ( lev, -dt );

    // Compute intermediate velocity
    mfix_compute_intermediate_velocity ( lev, dt, ro_go, ep_go );

    // Projection
    mac_projection -> apply_projection (u_g, v_g, w_g, ep_g, ro_g, dt );

    // Recover pressure
    const MultiFab& phi_mac = mac_projection -> get_phi(lev);
    MultiFab::Add (*p_g[lev], phi_mac, 0, 0, 1, p_g[lev] -> nGrow() );
 
}


//
// Perform the following operations:
//
//       u_g = u_g + coeff * ( dp_g/dx + dp0/dx ) * (1/ro_g)
//       v_g = v_g + coeff * ( dp_g/dy + dp0/dy ) * (1/ro_g)
//       w_g = w_g + coeff * ( dp_g/dz + dp0/dz ) * (1/ro_g)
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
// Compute uacc, vacc, and wacc by using u_g, v_g, and w_g
//
void
mfix_level::mfix_compute_fluid_acceleration ( int lev,
					      Vector< std::unique_ptr<MultiFab> >& u, 
					      Vector< std::unique_ptr<MultiFab> >& v,
					      Vector< std::unique_ptr<MultiFab> >& w,
					      Vector< std::unique_ptr<MultiFab> >& ro,
					      Vector< std::unique_ptr<MultiFab> >& ep )
{
    BL_PROFILE("mfix_level::mfix_compute_fluid_acceleration");

    int xdir = 1;
    int ydir = 2;
    int zdir = 3;

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*ep_g[lev],true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();
	const Box& sbx = (*ep_g[lev])[mfi].box();

	Box ubx((*u_g[lev])[mfi].box());
	Box vbx((*v_g[lev])[mfi].box());
	Box wbx((*w_g[lev])[mfi].box());

	calc_trd_g(sbx.loVect(), sbx.hiVect(),
		   ubx.loVect(), ubx.hiVect(),
		   vbx.loVect(), vbx.hiVect(),
		   wbx.loVect(), wbx.hiVect(),
		   bx.loVect(),  bx.hiVect(),
		   (*trD_g[lev])[mfi].dataPtr(),
		   (*u_g[lev])[mfi].dataPtr(),
		   (*v_g[lev])[mfi].dataPtr(),
		   (*w_g[lev])[mfi].dataPtr(),
		   &dx, &dy, &dz);
    }
    fill_mf_bc(lev,*trD_g[lev]);
    trD_g[lev] -> FillBoundary(geom[lev].periodicity());
    
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
            (*ro[lev])[mfi].dataPtr (),
	    (*ep[lev])[mfi].dataPtr (),
	    (*trD_g[lev])[mfi].dataPtr(),
	    ((*lambda_g[lev])[mfi]).dataPtr(),
	    geom[lev].CellSize (), &xdir, &upwind_ep );


	// y direction
	compute_fluid_acceleration (
	    BL_TO_FORTRAN_BOX(vbx),  
	    BL_TO_FORTRAN_ANYD((*vacc[lev])[mfi]),
	    (*slopes_v[lev])[mfi].dataPtr (),
	    BL_TO_FORTRAN_ANYD((*u[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*ro[lev])[mfi].dataPtr (),
	    (*ep[lev])[mfi].dataPtr (),
	    (*trD_g[lev])[mfi].dataPtr(),
	    ((*lambda_g[lev])[mfi]).dataPtr(),
	    geom[lev].CellSize (), &ydir, &upwind_ep );

	// z direction
	compute_fluid_acceleration (
	    BL_TO_FORTRAN_BOX(wbx),  
	    BL_TO_FORTRAN_ANYD((*wacc[lev])[mfi]),
	    (*slopes_w[lev])[mfi].dataPtr (),
	    BL_TO_FORTRAN_ANYD((*u[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
	    (*ro[lev])[mfi].dataPtr (),
	    (*ep[lev])[mfi].dataPtr (),
	    (*trD_g[lev])[mfi].dataPtr(),
	    ((*lambda_g[lev])[mfi]).dataPtr(),
	    geom[lev].CellSize (), &zdir, &upwind_ep );
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
				      Vector< std::unique_ptr<MultiFab> >& w,
				      Vector< std::unique_ptr<MultiFab> >& ro,
				      Vector< std::unique_ptr<MultiFab> >& ep)

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
		      BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
		      (*ep_g[lev])[mfi].dataPtr(),
		      domain.loVect (), domain.hiVect (),
		      geom[lev].CellSize (), &dt, &xdir,
		      &upwind_ep );

	add_forcing ( BL_TO_FORTRAN_BOX(vbx),  
		      BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*drag_v[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
		      (*ep_g[lev])[mfi].dataPtr(),
		      domain.loVect (), domain.hiVect (),
		      geom[lev].CellSize (), &dt, &ydir,
		      &upwind_ep );

	add_forcing ( BL_TO_FORTRAN_BOX(wbx),  
		      BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*drag_w[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
		      (*ep[lev])[mfi].dataPtr(),
		      domain.loVect (), domain.hiVect (),
		      geom[lev].CellSize (), &dt, &zdir,
		      &upwind_ep );
    }

}

//
// Implicit solve for the intermediate velocity.
// Currently this means accounting for the implicit part of the fluid/particle
// momentum exchange
// 
void
mfix_level::mfix_compute_intermediate_velocity ( int lev,
						 amrex::Real dt, 
						 Vector< std::unique_ptr<MultiFab> >& ro,
						 Vector< std::unique_ptr<MultiFab> >& ep )
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
					BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
					(*ep[lev])[mfi].dataPtr (),
					&xdir, &dt, &upwind_ep );

	compute_intermediate_velocity ( BL_TO_FORTRAN_BOX(vbx),  
					BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*f_gds_v[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
					(*ep[lev])[mfi].dataPtr (),
					&ydir, &dt, &upwind_ep );
	
	compute_intermediate_velocity ( BL_TO_FORTRAN_BOX(wbx),  
					BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*f_gds_w[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
					(*ep[lev])[mfi].dataPtr (),
					&zdir, &dt, &upwind_ep );
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
			   &nghost,
			   bc_ilo.dataPtr(), bc_ihi.dataPtr() );

	compute_v_slopes ( BL_TO_FORTRAN_BOX(vbx),
			   BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
			   (*slopes_v[lev])[mfi].dataPtr (),
			   domain.loVect (), domain.hiVect (),
			   &nghost,
			   bc_jlo.dataPtr(), bc_jhi.dataPtr() );

	compute_w_slopes ( BL_TO_FORTRAN_BOX(wbx),
			   BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
			   (*slopes_w[lev])[mfi].dataPtr (),
			   domain.loVect (), domain.hiVect (),
			   &nghost,
			   bc_klo.dataPtr(), bc_khi.dataPtr() );
    }

    // Fill halo cells
    slopes_u[lev] -> FillBoundary(geom[lev].periodicity());
    slopes_v[lev] -> FillBoundary(geom[lev].periodicity());
    slopes_w[lev] -> FillBoundary(geom[lev].periodicity());
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
    mfix_set_velocity_bcs (lev);
    
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
// Set the BCs for all the variables EXCEPT pressure or velocity.
// 
void
mfix_level::mfix_set_projection_bcs (int lev)
{
  BL_PROFILE("mfix_level::mfix_set_projection_bcs()");
  
  Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(*p_g[lev], true); mfi.isValid(); ++mfi)
    {
      set_projection_bcs ( BL_TO_FORTRAN_ANYD((*ep_g[lev])[mfi]),
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
}

//
// Fills ghost cell values of pressure appropriately for the BC type
//
void
mfix_level::mfix_extrap_pressure (int lev, std::unique_ptr<amrex::MultiFab>& p)
{
    BL_PROFILE("mfix_level::mfix_extrap_pressure()");

    Box domain(geom[lev].Domain());

    for (MFIter mfi(*p); mfi.isValid(); ++mfi) {

	extrap_pressure_to_ghost_cells (
	    BL_TO_FORTRAN_ANYD((*p)[mfi]),
	    bc_ilo.dataPtr(), bc_ihi.dataPtr(),
	    bc_jlo.dataPtr(), bc_jhi.dataPtr(),
	    bc_klo.dataPtr(), bc_khi.dataPtr(),
	    domain.loVect(), domain.hiVect(),
	    &nghost );
    }
}


void
mfix_level::check_for_nans (int lev)
{
    bool ug_has_nans = u_g[lev] -> contains_nan ();
    bool vg_has_nans = v_g[lev] -> contains_nan ();
    bool wg_has_nans = w_g[lev] -> contains_nan ();
    bool pg_has_nans = p_g[lev] -> contains_nan ();

    if (ug_has_nans)
	std::cout << "WARNING: u_g contains NaNs!!!";

    if (vg_has_nans)
	std::cout << "WARNING: v_g contains NaNs!!!";

    if (wg_has_nans)
	std::cout << "WARNING: w_g contains NaNs!!!";

    if (pg_has_nans)
	std::cout << "WARNING: p_g contains NaNs!!!";

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
