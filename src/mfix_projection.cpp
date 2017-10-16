#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>


// For multigrid
#include <AMReX_FMultiGrid.H>
#include <AMReX_stencil_types.H>






void
mfix_level::EvolveFluidProjection(int lev, int nstep, int set_normg,
				  Real dt, Real& prev_dt, Real time, Real normg)
{
    // Reimpose boundary conditions -- make sure to do this before we compute tau
    mfix_set_bc1(lev);

    // Calculate transport coefficients
    int calc_flag = 2;
    mfix_calc_coeffs(lev,calc_flag);
  
    // Back up field
    // Backup field variable to old
    int nghost = ep_go[lev]->nGrow();

    MultiFab::Copy(*ep_go[lev],  *ep_g[lev],  0, 0, 1, nghost);
    MultiFab::Copy( *p_go[lev],   *p_g[lev],  0, 0, 1, nghost);
    MultiFab::Copy(*ro_go[lev],  *ro_g[lev],  0, 0, 1, nghost);
    MultiFab::Copy(*rop_go[lev], *rop_g[lev], 0, 0, 1, nghost);
    MultiFab::Copy(*u_go[lev],   *u_g[lev],   0, 0, 1, nghost);
    MultiFab::Copy(*v_go[lev],   *v_g[lev],   0, 0, 1, nghost);
    MultiFab::Copy(*w_go[lev],   *w_g[lev],   0, 0, 1, nghost);


    // Compute the divergence of the velocity field, div(u).
    // div(u) is needed to compute the volumetric term in the
    // stress tensor
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*ep_g[lev],true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	compute_divu ( BL_TO_FORTRAN_BOX(bx),
		       BL_TO_FORTRAN_ANYD((*trD_g[lev])[mfi]),
		       BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
		       BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
		       BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
		       geom[lev].CellSize() );		     

    }
  
    fill_mf_bc(lev,*trD_g[lev]);



    // User hooks
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*ep_g[lev], true); mfi.isValid(); ++mfi)
	mfix_usr2();

    // Calculate drag coefficient
    if (solve_dem)
	mfix_calc_drag_fluid(lev);
  

    // Here we should check the CFL condition

  
    // Compute intermediate velocity
    mfix_compute_u_star ( lev, dt );
    mfix_compute_v_star ( lev, dt );
    mfix_compute_w_star ( lev, dt );
  
    // Reimpose boundary conditions
    mfix_set_bc1(lev);

    //  Do projection HERE
    mfix_apply_projection( lev, dt );
    
// Loop over iterate for auto time-step size adjustment
    int reiterate;
    do {
	prev_dt = dt;

	// Calculate bulk density (epg*ro_g) at cell faces
	mfix_conv_rop(lev);

	// Calculate face mass fluxes
	mfix_calc_mflux(lev);

	int converged=0;
	int nit=0;          // number of iterations

	///////////////// ---- call to iterate -------- /////////////////
	do {
	    nit++;

	    Real residuals[2*8];
	    for (int i=0; i<2*8; ++i)
		residuals[i] = 0.0L;

	    // Solve momentum equations in a thread safe way.
	    Real num_u = 0.0L;
	    Real num_v = 0.0L;
	    Real num_w = 0.0L;

	    Real denom_u = 0.0L;
	    Real denom_v = 0.0L;
	    Real denom_w = 0.0L;

	    mfix_solve_for_u(lev, dt, num_u, denom_u);
	    mfix_solve_for_v(lev, dt, num_v, denom_v);
	    mfix_solve_for_w(lev, dt, num_w, denom_w);

	    residuals[1] = num_u;
	    residuals[2] = num_v;
	    residuals[3] = num_w;

	    residuals[9]  = denom_u;
	    residuals[10] = denom_v;
	    residuals[11] = denom_w;

	    //Called after each momentum equation variable solve
	    MultiFab::Copy(*u_g[lev], *u_gt[lev], 0, 0, 1, u_g[lev]->nGrow());
	    MultiFab::Copy(*v_g[lev], *v_gt[lev], 0, 0, 1, v_g[lev]->nGrow());
	    MultiFab::Copy(*w_g[lev], *w_gt[lev], 0, 0, 1, w_g[lev]->nGrow());

	    u_g[lev]->FillBoundary(geom[lev].periodicity());
	    v_g[lev]->FillBoundary(geom[lev].periodicity());
	    w_g[lev]->FillBoundary(geom[lev].periodicity());

	    // Calculate transport coefficients
	    mfix_physical_prop(lev,0);

	    // Reimpose boundary conditions
	    mfix_set_bc1(lev);

	    // Calculate bulk density (epg*ro_g) at cell faces
	    mfix_conv_rop(lev);

	    // Solve the pressure correction equation
	    Real   num_p = 0.0L;
	    Real denom_p = 0.0L;
	    mfix_solve_for_pp(lev,dt,num_p,denom_p);
	    residuals[0] = num_p;
	    residuals[8] = denom_p;

	    // Apply pressure correction to all Pg, Ug, Vg, Wg
	    mfix_correct_0(lev);

	    // Update fluid density
	    mfix_physical_prop(lev,0);

	    // Calculate face mass fluxes
	    mfix_calc_mflux(lev);

	    // Check for convergence
	    ParallelDescriptor::ReduceRealSum(residuals,16);
	    converged = check_convergence(&nit, residuals);

	    // Display current iteration residuals
	    if ( ParallelDescriptor::IOProcessor() )
		display_resid(&time, &dt, &nit, residuals);

	} while(converged==0 && nit<max_nit);

	// Adjust time step if iteration failed.
	reiterate = adjustdt(&converged, &nit, &dt);
	if(reiterate == 1) {

	    // Reset the field variables
	    MultiFab::Copy(*ep_g[lev],  *ep_go[lev],  0, 0, 1, nghost);
	    MultiFab::Copy(*p_g[lev],   *p_go[lev],   0, 0, 1, nghost);
	    MultiFab::Copy(*ro_g[lev],  *ro_go[lev],  0, 0, 1, nghost);
	    MultiFab::Copy(*rop_g[lev], *rop_go[lev], 0, 0, 1, nghost);
	    MultiFab::Copy(*u_g[lev],   *u_go[lev],   0, 0, 1, nghost);
	    MultiFab::Copy(*v_g[lev],   *v_go[lev],   0, 0, 1, nghost);
	    MultiFab::Copy(*w_g[lev],   *w_go[lev],   0, 0, 1, nghost);

	}
    } while (reiterate==1);

}




void
mfix_level::mfix_compute_u_star (int lev, amrex::Real dt)
{
    BL_PROFILE("mfix_level::mfix_compute_u_star");

#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*u_go[lev],true); mfi.isValid(); ++mfi)
    {
	const Box&  bx = mfi.tilebox();

	// Compute convection term
	compute_ugradu_x ( BL_TO_FORTRAN_BOX(bx),  
			   BL_TO_FORTRAN_ANYD((*u_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*v_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*w_go[lev])[mfi]),
			   (*ugradu_x[lev])[mfi].dataPtr(),
			   geom[lev].CellSize() );

	// Compute diffusion term
	compute_divtau_x ( BL_TO_FORTRAN_BOX(bx),  
			   BL_TO_FORTRAN_ANYD((*u_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*v_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*w_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
			   (*lambda_g[lev])[mfi].dataPtr(),
			   (*trD_g[lev])[mfi].dataPtr(),
			   (*divtau_x[lev])[mfi].dataPtr(),
			   geom[lev].CellSize() );

	int dir = 1;
	compute_intermediate_velocity ( BL_TO_FORTRAN_BOX(bx),  
					BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
					(*u_go[lev])[mfi].dataPtr(),
					(*ugradu_x[lev])[mfi].dataPtr(),
					(*divtau_x[lev])[mfi].dataPtr(),
					(*drag_u[lev])[mfi].dataPtr(),
					(*f_gds_u[lev])[mfi].dataPtr(),
					BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
					&dt, &dir);


    }

}


void
mfix_level::mfix_compute_v_star (int lev, amrex::Real dt)
{
    BL_PROFILE("mfix_level::mfix_compute_v_star");
   
    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*v_go[lev],true); mfi.isValid(); ++mfi)
    {
	const Box&  bx = mfi.tilebox();

	// Compute convection term
	compute_ugradu_y ( BL_TO_FORTRAN_BOX(bx),  
			   BL_TO_FORTRAN_ANYD((*u_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*v_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*w_go[lev])[mfi]),
			   (*ugradu_y[lev])[mfi].dataPtr(),
			   geom[lev].CellSize() );

	// Compute diffusion term
	compute_divtau_y ( BL_TO_FORTRAN_BOX(bx),  
			   BL_TO_FORTRAN_ANYD((*u_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*v_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*w_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
			   (*lambda_g[lev])[mfi].dataPtr(),
			   (*trD_g[lev])[mfi].dataPtr(),
			   (*divtau_y[lev])[mfi].dataPtr(),
			   geom[lev].CellSize() );

	int dir = 2;
	compute_intermediate_velocity ( BL_TO_FORTRAN_BOX(bx),  
					BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
					(*v_go[lev])[mfi].dataPtr(),
					(*ugradu_y[lev])[mfi].dataPtr(),
					(*divtau_y[lev])[mfi].dataPtr(),
					(*drag_u[lev])[mfi].dataPtr(),
					(*f_gds_u[lev])[mfi].dataPtr(),
					BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
					&dt, &dir);


    }

    
}


void
mfix_level::mfix_compute_w_star (int lev, amrex::Real dt)
{
    BL_PROFILE("mfix_level::mfix_compute_w_star");
  
    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*w_go[lev],true); mfi.isValid(); ++mfi)
    {
	const Box&  bx = mfi.tilebox();

	// Compute convection term
	compute_ugradu_z ( BL_TO_FORTRAN_BOX(bx),  
			   BL_TO_FORTRAN_ANYD((*u_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*v_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*w_go[lev])[mfi]),
			   (*ugradu_z[lev])[mfi].dataPtr(),
			   geom[lev].CellSize() );

	// Compute diffusion term
	compute_divtau_z ( BL_TO_FORTRAN_BOX(bx),  
			   BL_TO_FORTRAN_ANYD((*u_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*v_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*w_go[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*mu_g[lev])[mfi]),
			   (*lambda_g[lev])[mfi].dataPtr(),
			   (*trD_g[lev])[mfi].dataPtr(),
			   (*divtau_z[lev])[mfi].dataPtr(),
			   geom[lev].CellSize() );

	int dir = 3;
	compute_intermediate_velocity ( BL_TO_FORTRAN_BOX(bx),  
					BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
					(*w_go[lev])[mfi].dataPtr(),
					(*ugradu_z[lev])[mfi].dataPtr(),
					(*divtau_z[lev])[mfi].dataPtr(),
					(*drag_u[lev])[mfi].dataPtr(),
					(*f_gds_u[lev])[mfi].dataPtr(),
					BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
					&dt, &dir);


    }
    
}

void 
mfix_level::mfix_apply_projection ( int lev, amrex::Real dt )
{
   BL_PROFILE("mfix_level::mfix_apply_projection");

   // Compute right hand side
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*ep_g[lev],true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	compute_divu ( BL_TO_FORTRAN_BOX(bx),
		       BL_TO_FORTRAN_ANYD((*trD_g[lev])[mfi]),
		       BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
		       BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
		       BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
		       geom[lev].CellSize() );		     

    }
  
    fill_mf_bc(lev,*trD_g[lev]);

    
}


//
// Solves - div ( b grad(phi) ) = rhs via Multigrid
//
// phi and rhs are cell-centered
// b           is  face-centered
// 
void
mfix_level::solve_poisson_equation ( int lev, MultiFab& b, MultiFab& phi, MultiFab& rhs )
{
    BL_PROFILE("mfix_level::solve_poisson_equation");

    // Multigrid inputs
    Vector<int>                         bc(2*AMREX_SPACEDIM, 2); // Neumann boundaries
    bool                                nodal = false;
    int                                 stencil =  amrex::CC_CROSS_STENCIL;
    bool                                have_rhcc = false;;
    int                                 nc = 0; // Don't know what it is but it should not
                                                // make any difference in the solve 
    int                                 verbose = 1;
    Real                                rel_tol = 1.0e-14;
    Real                                abs_tol = 1.0e-14;
    amrex::FMultiGrid                   solver(geom[lev]);

    solver.set_stencil (stencil);
    solver.set_verbose (verbose);
    solver.set_bc (bc.dataPtr());
    solver.set_const_gravity_coeffs ();

    solver.solve ( phi, rhs, rel_tol, abs_tol, 0, 0, 1 );

}

