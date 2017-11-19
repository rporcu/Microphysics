#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

// For multigrid
#include <AMReX_FMultiGrid.H>
#include <AMReX_stencil_types.H>


void
mfix_level::EvolveFluidProjection(int lev, int nstep, int set_normg,
				  Real& dt, Real& prev_dt, Real time, Real normg)
{

    amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";
      
    // Just for testing purposes, call a subroutine to initialize the
    // velocity field. The subroutine will not do anything except the first
    // time it is called, i.e. it will work only as initialization.
    // This fails if the code is restarted since it will initialize the velocity
    // field to the initial value.
    init_tests_projection (lev);

    check_for_nans (lev);
    
    // Extrapolate boundary values for density and volume fraction
    // The subsequent call to mfix_set_bc1 will only overwrite
    // rop_g and ep_g ghost values for PINF and POUT
    fill_mf_bc ( lev, *rop_g[lev] );
    fill_mf_bc ( lev, *ep_g[lev] );
    
    // Fill ghost nodes and reimpose boundary conditions
    u_g[lev] -> FillBoundary (geom[lev].periodicity());
    v_g[lev] -> FillBoundary (geom[lev].periodicity());
    w_g[lev] -> FillBoundary (geom[lev].periodicity());
    mfix_set_bc1(lev);

    fill_mf_bc ( lev, *p_g[lev] );
    
    // Calculate transport coefficients
    int calc_flag = 2;
    mfix_calc_coeffs (lev,calc_flag);
  
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
  

 
//     // User hooks
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//     for (MFIter mfi(*ep_g[lev], true); mfi.isValid(); ++mfi)
// 	mfix_usr2();

//     // Calculate drag coefficient
//     if (solve_dem)
// 	mfix_calc_drag_fluid(lev);
  

    // Here we should check the CFL condition
    // Compute dt for this time step
    Real umax  = u_g[lev] -> norm0 ();
    Real vmax  = v_g[lev] -> norm0 ();
    Real wmax  = w_g[lev] -> norm0 ();
    Real romin = rop_g[lev] -> min (0);
    Real mumax = mu_g[lev] -> max (0);
    
    compute_new_dt ( &umax, &vmax, &wmax, &romin, &mumax,
		     geom[lev].CellSize(), &dt );
    prev_dt = dt ;


    std::cout << "\nTentative velocity computation at  time = " << time
	      << " ( dt = "<< dt << " )\n";
    std::cout << "At beginning of time step :\n";
    std::cout << "max(abs(u))  = " << u_g[lev] -> norm0 () << "\n";
    std::cout << "max(abs(v))  = " << v_g[lev] -> norm0 () << "\n";	
    std::cout << "max(abs(w))  = " << w_g[lev] -> norm0 () << "\n";

    //
    // Time integration step
    //
    
    // Step 1: compute u* (predictor step) and store it in u_g,
    // v_g, and w_g
    mfix_compute_velocity_slopes ( lev ); 
    mfix_compute_fluid_acceleration ( lev ); 
    mfix_apply_pcm_prediction ( lev, dt );

    std::cout << "\n Fluid accelerations :\n";
    std::cout << "max(abs(ru))  = " << uacc[lev] -> norm0 () << "\n";
    std::cout << "max(abs(rv))  = " << vacc[lev] -> norm0 () << "\n";	
    std::cout << "max(abs(rw))  = " << wacc[lev] -> norm0 () << "\n\n";

  
    // Exchange halo nodes and apply BCs
    u_g[lev] -> FillBoundary (geom[lev].periodicity());
    v_g[lev] -> FillBoundary (geom[lev].periodicity());
    w_g[lev] -> FillBoundary (geom[lev].periodicity());
    mfix_set_bc1 ( lev );
    
    // Step 2: compute u** (corrector step) and store it in
    // u_g, v_g, and w_g.
    mfix_compute_velocity_slopes ( lev );
    mfix_compute_fluid_acceleration ( lev );     
    mfix_apply_pcm_correction ( lev, dt );


    // Add forcing terms ( gravity and/or momentum
    // exchange with particles )
//    mfix_apply_forcing_terms ( lev, dt );

    // Fill ghost cells and reimpose boundary conditions
    u_g[lev] -> FillBoundary (geom[lev].periodicity());
    v_g[lev] -> FillBoundary (geom[lev].periodicity());
    w_g[lev] -> FillBoundary (geom[lev].periodicity());
    mfix_set_bc1(lev);

    //
    std::cout << "\nBefore projection step :\n";
    std::cout << "max(abs(u))  = " << u_g[lev] -> norm0 () << "\n";
    std::cout << "max(abs(v))  = " << v_g[lev] -> norm0 () << "\n";	
    std::cout << "max(abs(w))  = " << w_g[lev] -> norm0 () << "\n";
    std::cout << "max(abs(p))  = " << p_g[lev] -> norm0 () << "\n\n";

    check_for_nans (lev);

    // 
    //  Projection Step
    // 
    mfix_apply_projection ( lev, dt ); 
    
    //
    // std::cout << "\nAfter projection step :\n";
    std::cout << "Final max(abs(u))  = " << u_g[lev] -> norm0 () << "\n";
    std::cout << "Final max(abs(v))  = " << v_g[lev] -> norm0 () << "\n";	
    std::cout << "Final max(abs(w))  = " << w_g[lev] -> norm0 () << "\n";
    std::cout << "Final max(abs(p))  = " << p_g[lev] -> norm0 () << "\n\n";
    // std::cout << "\nAfter projection step : ";
    // std::cout << "max(abs(u)), max(abs(v)),max(abs(w)) =  ";
    // std::cout << u_g[lev] -> norm0 () << "   ";
    // std::cout << v_g[lev] -> norm0 () << "   ";	
    // std::cout << w_g[lev] -> norm0 () << "\n \n";
   
    
    check_for_nans (lev);
    
    // // Calculate transport coefficients
    // mfix_physical_prop(lev,0);

    // // Update fluid density
    // mfix_physical_prop(lev,0);



    // Compute the divergence of the velocity field, div(u).
    // to check if div(u) = 0 is satisfied
    u_g[lev] -> FillBoundary (geom[lev].periodicity());
    v_g[lev] -> FillBoundary (geom[lev].periodicity());
    w_g[lev] -> FillBoundary (geom[lev].periodicity());
    mfix_set_bc1(lev);
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*trD_g[lev],true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	compute_divu ( BL_TO_FORTRAN_BOX(bx),
		       BL_TO_FORTRAN_ANYD((*trD_g[lev])[mfi]),
		       BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
		       BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
		       BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
		       geom[lev].CellSize() );		     

    }

    // BCs 
    fill_mf_bc(lev,*trD_g[lev]);

    amrex::Print() << "Max(abs(divu)) = "<< trD_g[lev] -> norm0 () << "\n";
}

//
//
// 
void
mfix_level::init_tests_projection (int lev) 
{

    static bool first_access = true;

    if (!first_access) return;

    // Read in direction
    ParmParse pp;
    int vp = 1;			//  1=xy 2=xz 3=yz
    pp.query ("vortices_plane", vp );

    std::cout << "2D vortices initialized on plane " << vp << std::endl;
    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*p_go[lev],true); mfi.isValid(); ++mfi)
    {
	Box ubx = amrex::convert(mfi.tilebox(),e_x);
	Box vbx = amrex::convert(mfi.tilebox(),e_y);
	Box wbx = amrex::convert(mfi.tilebox(),e_z);

	init_periodic_vortices ( BL_TO_FORTRAN_BOX(ubx),  
				 BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
				 BL_TO_FORTRAN_BOX(vbx),  
				 BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
				 BL_TO_FORTRAN_BOX(wbx),  
				 BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
				 geom[lev].CellSize (),
				 (geom[lev].Domain ()).loVect (),
				 &vp ); 
	
    }

    std::cout << "Entering the init routine\n";


    u_g[lev] -> FillBoundary (geom[lev].periodicity());
    v_g[lev] -> FillBoundary (geom[lev].periodicity());
    w_g[lev] -> FillBoundary (geom[lev].periodicity());

    //mfix_set_bc1 (lev);

    first_access = false;
	
}



//
// Computes:
//
//           u_g = u_g + dt * R_u - dt * (dp/dx) / ro_g
//           v_g = v_g + dt * R_v - dt * (dp/dy) / ro_g
//           w_g = w_g + dt * R_w - dt * (dp/dz) / ro_g 
//
//  This is the prediction step of the Heun's integration
//  scheme, AKA Predictor-Corrector Method (PCM)
// 
void
mfix_level::mfix_apply_pcm_prediction (int lev, amrex::Real dt)
{
    BL_PROFILE("mfix_level::mfix_apply_pcm_prediction");


    // First add the fluid acceleration
    MultiFab::Saxpy (*u_g[lev], dt, *uacc[lev], 0, 0, 1, 0);
    MultiFab::Saxpy (*v_g[lev], dt, *vacc[lev], 0, 0, 1, 0);
    MultiFab::Saxpy (*w_g[lev], dt, *wacc[lev], 0, 0, 1, 0);

    // The add the pressure gradient
    mfix_add_pressure_gradient ( lev, -dt ); 

}


//
// Computes:
//
//   u_g = 0.5 * ( u_g + u_go + dt * uacc) 
//   v_g = 0.5 * ( v_g + v_go + dt * vacc) 
//   w_g = 0.5 * ( w_g + w_go + dt * wacc) 
//
//  This is the correction step of the Heun's integration
//  scheme, AKA Predictor-Corrector Method (PCM).
//  
void
mfix_level::mfix_apply_pcm_correction (int lev, amrex::Real dt)
{
    BL_PROFILE("mfix_level::mfix_apply_pcm_correction");

    // First add the fluid acceleration
    MultiFab::Saxpy (*u_g[lev], dt, *uacc[lev], 0, 0, 1, 0);
    MultiFab::Saxpy (*v_g[lev], dt, *vacc[lev], 0, 0, 1, 0);
    MultiFab::Saxpy (*w_g[lev], dt, *wacc[lev], 0, 0, 1, 0);

    // Then add the old velocity
    MultiFab::Add (*u_g[lev], *u_go[lev], 0, 0, 1, 0);
    MultiFab::Add (*v_g[lev], *v_go[lev], 0, 0, 1, 0);
    MultiFab::Add (*w_g[lev], *w_go[lev], 0, 0, 1, 0);
    
    // Multiply result by 0.5
    u_g[lev] -> mult ( 0.5, 0 );
    v_g[lev] -> mult ( 0.5, 0 );
    w_g[lev] -> mult ( 0.5, 0 );
}



//
// Perform the following operations:
//
//       u_g = u_g + coeff * ( dp_g/dx ) * (1/ro_g)
//       v_g = v_g + coeff * ( dp_g/dy ) * (1/ro_g)
//       w_g = w_g + coeff * ( dp_g/dz ) * (1/ro_g)
//
// 1/ro_g is stored in oro_g[lev][<0,1,2>] 
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
	Box ubx = amrex::convert ( mfi.tilebox(), e_x );
	Box vbx = amrex::convert ( mfi.tilebox(), e_y );
	Box wbx = amrex::convert ( mfi.tilebox(), e_z );

	add_gradient (
	    BL_TO_FORTRAN_BOX(ubx),  
	    BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
	    (*(oro_g[lev][0]))[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*p_g[lev])[mfi]),
	    geom[lev].CellSize(), &coeff, &xdir );
	
	add_gradient (
	    BL_TO_FORTRAN_BOX(vbx),  
	    BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
	    (*(oro_g[lev][1]))[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*p_g[lev])[mfi]),
	    geom[lev].CellSize(), &coeff, &ydir );
	
	add_gradient (
	    BL_TO_FORTRAN_BOX(wbx),  
	    BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
	    (*(oro_g[lev][2]))[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*p_g[lev])[mfi]),
	    geom[lev].CellSize(), &coeff, &zdir );

	
    }
}



//
// Compute uacc, vacc, and wacc by using u_g, v_g, and w_g
//
void
mfix_level::mfix_compute_fluid_acceleration (int lev)
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
	Box ubx = amrex::convert ( mfi.tilebox(), e_x );
	Box vbx = amrex::convert ( mfi.tilebox(), e_y );
	Box wbx = amrex::convert ( mfi.tilebox(), e_z );

	// x direction
	compute_fluid_acceleration (
	    BL_TO_FORTRAN_BOX(ubx),  
	    BL_TO_FORTRAN_ANYD((*uacc[lev])[mfi]),
	    (*slopes_u[lev])[mfi].dataPtr (),
	    BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*rop_g[lev])[mfi].dataPtr (),
	    geom[lev].CellSize (), &xdir );

	// y direction
	compute_fluid_acceleration (
	    BL_TO_FORTRAN_BOX(vbx),  
	    BL_TO_FORTRAN_ANYD((*vacc[lev])[mfi]),
	    (*slopes_v[lev])[mfi].dataPtr (),
	    BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*rop_g[lev])[mfi].dataPtr (),
	    geom[lev].CellSize (), &ydir );

	// z direction
	compute_fluid_acceleration (
	    BL_TO_FORTRAN_BOX(wbx),  
	    BL_TO_FORTRAN_ANYD((*wacc[lev])[mfi]),
	    (*slopes_w[lev])[mfi].dataPtr (),
	    BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*rop_g[lev])[mfi].dataPtr (),
	    geom[lev].CellSize (), &zdir );
    }
}




void
mfix_level::mfix_apply_forcing_terms (int lev, amrex::Real dt)
{
    BL_PROFILE("mfix_level::mfix_apply_forcing_terms");

#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi)
    {
	// Boxes for staggered components
	Box ubx = amrex::convert ( mfi.tilebox(), e_x );
	Box vbx = amrex::convert ( mfi.tilebox(), e_y );
	Box wbx = amrex::convert ( mfi.tilebox(), e_z );

	// To be filled
	
    }
}

//
// Compute the slopes of each velocity component in the
// three directions.
// 
void
mfix_level::mfix_compute_velocity_slopes (int lev)
{
    BL_PROFILE("mfix_level::mfix_compute_velocity_slopes");

#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi)
    {
	// Boxes for staggered components
	Box domain(geom[lev].Domain());
	Box ubx = amrex::convert ( mfi.tilebox(), e_x );
	Box vbx = amrex::convert ( mfi.tilebox(), e_y );
	Box wbx = amrex::convert ( mfi.tilebox(), e_z );

	compute_u_slopes ( BL_TO_FORTRAN_BOX(ubx),
			   BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
			   (*slopes_u[lev])[mfi].dataPtr (),
			   domain.loVect (), domain.hiVect (),
			   bc_ilo.dataPtr(), bc_ihi.dataPtr() );

	compute_v_slopes ( BL_TO_FORTRAN_BOX(vbx),
			   BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
			   (*slopes_v[lev])[mfi].dataPtr (),
			   domain.loVect (), domain.hiVect (),
			   bc_jlo.dataPtr(), bc_jhi.dataPtr() );

	compute_w_slopes ( BL_TO_FORTRAN_BOX(wbx),
			   BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
			   (*slopes_w[lev])[mfi].dataPtr (),
			   domain.loVect (), domain.hiVect (),
			   bc_klo.dataPtr(), bc_khi.dataPtr() );
    }

    // Fill halo cells
    slopes_u[lev] -> FillBoundary(geom[lev].periodicity());
    slopes_v[lev] -> FillBoundary(geom[lev].periodicity());
    slopes_w[lev] -> FillBoundary(geom[lev].periodicity());
}



void 
mfix_level::mfix_apply_projection ( int lev, amrex::Real dt )
{
   BL_PROFILE("mfix_level::mfix_apply_projection");

   // Compute right hand side, AKA div(u)
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*trD_g[lev],true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	compute_divu ( BL_TO_FORTRAN_BOX(bx),
		       BL_TO_FORTRAN_ANYD((*trD_g[lev])[mfi]),
		       BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
		       BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
		       BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
		       geom[lev].CellSize() );		     

    }

    // RHS BCs (probably we do not need this)
    trD_g[lev] -> FillBoundary(geom[lev].periodicity());

    // Compute the PPE coefficients
    mfix_compute_oro_g ( lev );
    // For the time being set oro_g to 1
    
    
    // Solve PPE
    solve_poisson_equation ( lev, oro_g, p_g, trD_g );
    p_g[lev] -> FillBoundary(geom[lev].periodicity());
    
    // Apply pressure correction
    //  p_g is now p_g = p^{n+1}*dt/2
    Real coeff = -1.0;
    int  xdir  = 1;
    int  ydir  = 2;
    int  zdir  = 3;
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*phi[lev],true); mfi.isValid(); ++mfi)
    {
	Box ubx = amrex::convert ( mfi.tilebox(), e_x );
	Box vbx = amrex::convert ( mfi.tilebox(), e_y );
	Box wbx = amrex::convert ( mfi.tilebox(), e_z );

	add_gradient (
	    BL_TO_FORTRAN_BOX(ubx),  
	    BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
	    (*(oro_g[lev][0]))[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*p_g[lev])[mfi]),
	    geom[lev].CellSize(), &coeff, &xdir );
	
	add_gradient (
	    BL_TO_FORTRAN_BOX(vbx),  
	    BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
	    (*(oro_g[lev][1]))[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*p_g[lev])[mfi]),
	    geom[lev].CellSize(), &coeff, &ydir );

	add_gradient (
	    BL_TO_FORTRAN_BOX(wbx),  
	    BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
	    (*(oro_g[lev][2]))[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*p_g[lev])[mfi]),
	    geom[lev].CellSize(), &coeff, &zdir );

    }


    // 
    // Rescale pressure p_g = p_g * 2 / dt
    // 
    int nghost = p_g[lev] -> nGrow ();
    p_g[lev] -> mult ( 2.0/dt, nghost );
    p_g[lev] -> FillBoundary(geom[lev].periodicity());
    
 }


//
// Solves   div ( b grad(phi) ) = rhs   via Multigrid
//
// phi and rhs are cell-centered
// b           is  face-centered
//
void
mfix_level::solve_poisson_equation (  int lev,
				      Vector< Vector< std::unique_ptr<MultiFab> > >& b,
				      Vector< std::unique_ptr<MultiFab> >& phi,
				      Vector< std::unique_ptr<MultiFab> >& rhs )
{
    BL_PROFILE("mfix_level::solve_poisson_equation");

    // Multigrid inputs
    Vector<int>                         bc(2*AMREX_SPACEDIM, -1); // Periodic boundaries
    int                                 stencil =  amrex::CC_CROSS_STENCIL;
    int                                 verbose = 0;
    Real                                rel_tol = 1.0e-13;
    Real                                abs_tol = 1.0e-14;
    amrex::FMultiGrid                   solver(geom[lev]);

    solver.set_stencil (stencil);
    solver.set_verbose (verbose);
    solver.set_bc (bc.dataPtr());

    
    // bool has_nans = rhs[lev] -> contains_nan ();

    // amrex::Print() << " HAS NANs = " << has_nans << "\n";
   
    solver.set_gravity_coeffs ( amrex::GetVecOfPtrs ( b[lev] ) );
    phi[lev] -> setVal (0.);
 	
    solver.solve ( *phi[lev], *rhs[lev], rel_tol, abs_tol, 0, 0, 1 );

}




//
// Computes 1/ro_g = ep_g/rop_g at the faces of the scalar cells
// 
void
mfix_level::mfix_compute_oro_g (int lev)
{
    BL_PROFILE("mfix_level::mfix_compute_oro_g");

    // Directions
    int xdir = 1;
    int ydir = 2;
    int zdir = 3;

    // For now, set everything to 1
    oro_g[lev][0] -> setVal (1.0);
    oro_g[lev][1] -> setVal (1.0);
    oro_g[lev][2] -> setVal (1.0);
    
// #ifdef _OPENMP
// #pragma omp parallel 
// #endif
//     for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi)
//     {
// 	// Boxes for staggered components
// 	Box ubx = amrex::convert ( mfi.tilebox(), e_x );
// 	Box vbx = amrex::convert ( mfi.tilebox(), e_y );
// 	Box wbx = amrex::convert ( mfi.tilebox(), e_z );

// 	// X direction
// 	compute_oro_g (BL_TO_FORTRAN_BOX(ubx),
// 		       BL_TO_FORTRAN_ANYD((*(oro_g[lev][0]))[mfi]),
// 		       BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
// 		       (*ep_g[lev])[mfi].dataPtr(), &xdir );

// 	// Y direction
// 	compute_oro_g (BL_TO_FORTRAN_BOX(vbx),
// 		       BL_TO_FORTRAN_ANYD((*(oro_g[lev][1]))[mfi]),
// 		       BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
// 		       (*ep_g[lev])[mfi].dataPtr(), &ydir );

// 	// Z direction
// 	compute_oro_g (BL_TO_FORTRAN_BOX(wbx),
// 		       BL_TO_FORTRAN_ANYD((*(oro_g[lev][2]))[mfi]),
// 		       BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
// 		       (*ep_g[lev])[mfi].dataPtr(), &zdir );
	
	
//     }

    oro_g[lev][0] -> FillBoundary(geom[lev].periodicity());
    oro_g[lev][1] -> FillBoundary(geom[lev].periodicity());
    oro_g[lev][2] -> FillBoundary(geom[lev].periodicity());

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
// Estimate pressure by applying projection for one time step
// and then rolling back velocities:
//
//           u_g = u_g + dt * R_u - dt * (dp/dx) / ro_g
//           v_g = v_g + dt * R_v - dt * (dp/dy) / ro_g
//           w_g = w_g + dt * R_w - dt * (dp/dz) / ro_g 
//
//  This is the prediction step of the Heun's integration
//  scheme, AKA Predictor-Corrector Method (PCM)
// 
// void
// mfix_level::mfix_compute_initial_pressure (int lev, amrex::Real dt)
// {
//     BL_PROFILE("mfix_level::mfix_compute_initial_pressure");


//     // First add the fluid acceleration
//     MultiFab::Saxpy (*u_g[lev], dt, *uacc[lev], 0, 0, 1, 0);
//     MultiFab::Saxpy (*v_g[lev], dt, *vacc[lev], 0, 0, 1, 0);
//     MultiFab::Saxpy (*w_g[lev], dt, *wacc[lev], 0, 0, 1, 0);

//     // The add the pressure gradient
//     mfix_add_pressure_gradient ( lev, -dt ); 

// }
