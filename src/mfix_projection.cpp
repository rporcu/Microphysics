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
    //mfix_set_bc1(lev);

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

   

    // Compute the divergence of the velocity field, div(u).
    // div(u) is needed to compute the volumetric term in the
    // stress tensor
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
    fill_mf_bc ( lev, *trD_g[lev] );

    amrex::Print() << " Initial max(abs(div(u)))  = " << trD_g[lev] ->  norm0 () << "\n"; 

 
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
    std::cout << "Before predictor step :\n";
    std::cout << "max(abs(u))  = " << u_g[lev] -> norm0 () << "\n";
    std::cout << "max(abs(v))  = " << v_g[lev] -> norm0 () << "\n";	
    std::cout << "max(abs(w))  = " << w_g[lev] -> norm0 () << "\n";

    //
    // Time integration step
    //
    
    // Step 1: compute u* (predictor step) and store it in u_gt,
    // v_gt, and w_gt
    mfix_apply_pcm_prediction ( lev, dt );

    // Exchange halo nodes and apply BCs
    u_gt[lev] -> FillBoundary (geom[lev].periodicity());
    v_gt[lev] -> FillBoundary (geom[lev].periodicity());
    w_gt[lev] -> FillBoundary (geom[lev].periodicity());
  
    // Here we should call set_bc1 () but it works for
    // u)g, v_g, w_g only. Gotta take care of that

    // Step 2: compute u** (corrector step) and store it in
    //  u_g, v_g, and w_g
    mfix_apply_pcm_correction ( lev, dt );
    
    // Add forcing terms ( previous pressure, gravity, and momentum
    // exchange with particles )
    mfix_apply_forcing_terms ( lev, dt );

    // Fill ghost cells and reimpose boundary conditions
    u_g[lev] -> FillBoundary (geom[lev].periodicity());
    v_g[lev] -> FillBoundary (geom[lev].periodicity());
    w_g[lev] -> FillBoundary (geom[lev].periodicity());
    // mfix_set_bc1(lev);

    //
    std::cout << "\nAfter predictor step :\n";
    std::cout << "max(abs(u))  = " << u_g[lev] -> norm0 () << "\n";
    std::cout << "max(abs(v))  = " << v_g[lev] -> norm0 () << "\n";	
    std::cout << "max(abs(w))  = " << w_g[lev] -> norm0 () << "\n\n";

    
    check_for_nans (lev);

    // 
    //  Projection Step
    // 
    mfix_apply_projection ( lev, dt );

    //
    std::cout << "\nAfter projection step :\n";
    std::cout << "max(abs(u))  = " << u_g[lev] -> norm0 () << "\n";
    std::cout << "max(abs(v))  = " << v_g[lev] -> norm0 () << "\n";	
    std::cout << "max(abs(w))  = " << w_g[lev] -> norm0 () << "\n \n";

    
    
    check_for_nans (lev);
    
    // // Calculate transport coefficients
    // mfix_physical_prop(lev,0);

    // // Update fluid density
    // mfix_physical_prop(lev,0);



    // Compute the divergence of the velocity field, div(u).
    // div(u) is needed to compute the volumetric term in the
    // stress tensor
    u_g[lev] -> FillBoundary (geom[lev].periodicity());
    v_g[lev] -> FillBoundary (geom[lev].periodicity());
    w_g[lev] -> FillBoundary (geom[lev].periodicity());
   
    // mfix_set_bc1(lev);
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

    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*p_go[lev],true); mfi.isValid(); ++mfi)
    {
	Box ubx = amrex::convert(mfi.tilebox(),e_x);
	Box vbx = amrex::convert(mfi.tilebox(),e_y);
	Box wbx = amrex::convert(mfi.tilebox(),e_z);

	init_periodic_vorteces ( BL_TO_FORTRAN_BOX(ubx),  
				 BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
				 BL_TO_FORTRAN_BOX(vbx),  
				 BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
				 BL_TO_FORTRAN_BOX(wbx),  
				 BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
				 geom[lev].CellSize (),
				 (geom[lev].Domain ()).loVect ()); 
	
    }

    std::cout << "Entering the init routine\n";


    u_g[lev] -> FillBoundary (geom[lev].periodicity());
    v_g[lev] -> FillBoundary (geom[lev].periodicity());
    w_g[lev] -> FillBoundary (geom[lev].periodicity());

    //mfix_set_bc1 (lev);

    first_access = false;
	
}


void
mfix_level::mfix_apply_pcm_prediction (int lev, amrex::Real dt)
{
    BL_PROFILE("mfix_level::mfix_apply_pcm_prediction");

#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi)
    {
	// Boxes for staggered components
	Box ubx = amrex::convert ( mfi.tilebox(), e_x );
	Box vbx = amrex::convert ( mfi.tilebox(), e_y );
	Box wbx = amrex::convert ( mfi.tilebox(), e_z );
	
        // U velocity ( u* stored in u_gt )
	int dir = 1;
	
	apply_pcm_prediction (
	    BL_TO_FORTRAN_BOX(ubx),  
	    BL_TO_FORTRAN_ANYD((*u_gt[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*u_go[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v_go[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w_go[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*rop_g[lev])[mfi].dataPtr (),
	    geom[lev].CellSize (), &dt, &dir );

	// V velocity ( v* stored in v_gt )
	dir = 2;
	
	apply_pcm_prediction (
	    BL_TO_FORTRAN_BOX(vbx),  
	    BL_TO_FORTRAN_ANYD((*v_gt[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*u_go[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v_go[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w_go[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*rop_g[lev])[mfi].dataPtr (),
	    geom[lev].CellSize (), &dt, &dir );

	// W velocity ( w* stored in w_gt )
	dir = 3;
	
	apply_pcm_prediction (
	    BL_TO_FORTRAN_BOX(wbx),  
	    BL_TO_FORTRAN_ANYD((*w_gt[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*u_go[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v_go[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w_go[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*rop_g[lev])[mfi].dataPtr (),
	    geom[lev].CellSize (), &dt, &dir );
    }
}

void
mfix_level::mfix_apply_pcm_correction (int lev, amrex::Real dt)
{
    BL_PROFILE("mfix_level::mfix_apply_pcm_correction");

#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi)
    {
	// Boxes for staggered components
	Box ubx = amrex::convert ( mfi.tilebox(), e_x );
	Box vbx = amrex::convert ( mfi.tilebox(), e_y );
	Box wbx = amrex::convert ( mfi.tilebox(), e_z );
	
        // U velocity ( u* stored in u_gt )
	int dir = 1;
	
	apply_pcm_correction (
	    BL_TO_FORTRAN_BOX(ubx),  
	    BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
	    (*u_go[lev])[mfi].dataPtr (),
	    BL_TO_FORTRAN_ANYD((*u_gt[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v_gt[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w_gt[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*rop_g[lev])[mfi].dataPtr (),
	    geom[lev].CellSize (), &dt, &dir );

        // V velocity ( v* stored in v_gt )
	dir = 2;
	
	apply_pcm_correction (
	    BL_TO_FORTRAN_BOX(vbx),  
	    BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
	    (*v_go[lev])[mfi].dataPtr (),
	    BL_TO_FORTRAN_ANYD((*u_gt[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v_gt[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w_gt[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*rop_g[lev])[mfi].dataPtr (),
	    geom[lev].CellSize (), &dt, &dir );
	
        // V velocity ( v* stored in v_gt )
	dir = 3;
	
	apply_pcm_correction (
	    BL_TO_FORTRAN_BOX(wbx),  
	    BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
	    (*w_go[lev])[mfi].dataPtr (),
	    BL_TO_FORTRAN_ANYD((*u_gt[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*v_gt[lev])[mfi]),
	    BL_TO_FORTRAN_ANYD((*w_gt[lev])[mfi]),
            BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
            (*rop_g[lev])[mfi].dataPtr (),
	    geom[lev].CellSize (), &dt, &dir );

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

	amrex::Real mdto2 = - dt*0.5; 

	int dir = 1;

	add_gradient (
	    BL_TO_FORTRAN_BOX(ubx),  
	    BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
	    (*(oro_g[lev][0]))[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*p_g[lev])[mfi]),
	    geom[lev].CellSize(), &mdto2, &dir );

	dir = 2 ;
	
	add_gradient (
	    BL_TO_FORTRAN_BOX(vbx),  
	    BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
	    (*(oro_g[lev][1]))[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*p_g[lev])[mfi]),
	    geom[lev].CellSize(), &mdto2, &dir );

	dir = 3;
	
	add_gradient (
	    BL_TO_FORTRAN_BOX(wbx),  
	    BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
	    (*(oro_g[lev][2]))[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*p_g[lev])[mfi]),
	    geom[lev].CellSize(), &mdto2, &dir );

	
    }
}
	


void
mfix_level::mfix_compute_velocity_slopes (int lev)
{
    BL_PROFILE("mfix_level::mfix_velocity_slopes");

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

   // Compute right hand side
   amrex::Real offset;
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

//    amrex::ParallelDescriptor::ReduceRealSum (offset);
//    trD_g[lev] -> plus ( -offset, 1 );

    // RHS BCs (probably we do not need this)
    trD_g[lev] -> FillBoundary(geom[lev].periodicity());

    // Compute the PPE coefficients
    compute_oro_g ( lev );
    
    // Solve PPE
    solve_poisson_equation ( lev, oro_g, phi, trD_g );
    phi[lev] -> FillBoundary(geom[lev].periodicity());
    p_g[lev] -> FillBoundary(geom[lev].periodicity());

    // Update pressure
    amrex::Real scale = dt * 0.5;
	
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox ();

	update_pressure (
	    BL_TO_FORTRAN_BOX(bx),  
	    BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
	    (*p_g[lev])[mfi].dataPtr (), &scale );	    
	
    }
    p_g[lev] -> FillBoundary(geom[lev].periodicity());

    std::cout << "max(abs(p_g)) = " << p_g[lev] -> norm0 () << "\n";


    
    // Apply pressure correction
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
	    BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
	    geom[lev].CellSize(), &coeff, &xdir );

	
	add_gradient (
	    BL_TO_FORTRAN_BOX(vbx),  
	    BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
	    (*(oro_g[lev][1]))[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
	    geom[lev].CellSize(), &coeff, &ydir );

	add_gradient (
	    BL_TO_FORTRAN_BOX(wbx),  
	    BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
	    (*(oro_g[lev][2]))[mfi].dataPtr(),
	    BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
	    geom[lev].CellSize(), &coeff, &zdir );


	
	// apply_pressure_correction ( BL_TO_FORTRAN_BOX(ubx),  
	// 			    BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
	// 			    (*(oro_g[lev][0]))[mfi].dataPtr(),
	// 			    BL_TO_FORTRAN_BOX(vbx),  
	// 			    BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
	// 			    (*(oro_g[lev][1]))[mfi].dataPtr(),
	// 			    BL_TO_FORTRAN_BOX(wbx),  
	// 			    BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
	// 			    (*(oro_g[lev][2]))[mfi].dataPtr(),
	// 			    BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
	// 			    geom[lev].CellSize (),
	// 			    &dt );
    }

}


//
// Solves - div ( b grad(phi) ) = rhs via Multigrid
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
    bool                                nodal = false;
    int                                 stencil =  amrex::CC_CROSS_STENCIL;
    bool                                have_rhcc = false;
    int                                 nc = 0; // Don't know what it is but it should not
                                                // make any difference in the solve 
    int                                 verbose = 0;
    Real                                rel_tol = 1.0e-13;
    Real                                abs_tol = 1.0e-14;
    amrex::FMultiGrid                   solver(geom[lev]);

    solver.set_stencil (stencil);
    solver.set_verbose (verbose);
    solver.set_bc (bc.dataPtr());
    // solver.set_maxorder ();
    
    amrex::Print() << "RHS HAS NANS = " << rhs[lev] -> contains_nan () << "\n";
    amrex::Print() << "norm0(ppe rhs) = " << rhs[lev]  -> norm0 () << "\n";
    // amrex::Print() << "norm0(b_x) = "     << b[lev][0] -> norm0 () << "\n";
    // amrex::Print() << "norm0(b_y) = "     << b[lev][1] -> norm0 () << "\n";
    // amrex::Print() << "norm0(b_z) = "     << b[lev][2] -> norm0 () << "\n";

    bool has_nans = rhs[lev] -> contains_nan ();

    amrex::Print() << " HAS NANs = " << has_nans << "\n";
    
    solver.set_mac_coeffs ( amrex::GetVecOfPtrs ( b[lev] ) );
    phi[lev] -> setVal (0.);
    
//    VisMF::Write(*rhs[lev], "rhs");
	
    solver.solve ( *phi[lev], *rhs[lev], rel_tol, abs_tol, 0, 0, 1 );

}


void
mfix_level::compute_oro_g (int lev)
{
    BL_PROFILE("mfix_level::compute_oro_g");

    // Compute the PPE coefficients
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*u_g[lev],true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	compute_oro_g_x (BL_TO_FORTRAN_BOX(bx),
			 BL_TO_FORTRAN_ANYD((*(oro_g[lev][0]))[mfi]),
			 BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
			 (*ep_g[lev])[mfi].dataPtr() );
    }

    oro_g[lev][0] -> FillBoundary(geom[lev].periodicity());
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*v_g[lev],true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	compute_oro_g_y (BL_TO_FORTRAN_BOX(bx),
			 BL_TO_FORTRAN_ANYD((*(oro_g[lev][1]))[mfi]),
			 BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
			 (*ep_g[lev])[mfi].dataPtr() );
    }

    oro_g[lev][1] -> FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*w_g[lev],true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	compute_oro_g_z (BL_TO_FORTRAN_BOX(bx),
			 BL_TO_FORTRAN_ANYD((*(oro_g[lev][2]))[mfi]),
			 BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
			 (*ep_g[lev])[mfi].dataPtr() );
    }
    
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
