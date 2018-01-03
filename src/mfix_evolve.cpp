#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

// This subroutine is the driver for the whole time stepping (fluid + particles )
void
mfix_level::Evolve(int lev, int nstep, int set_normg, int steady_state,  Real &dt, Real& prev_dt,
                   Real time, Real normg)
{
  AMREX_ALWAYS_ASSERT(lev == 0);

  Real sum_vol;
  if (solve_dem && solve_fluid)
  {
    mfix_calc_volume_fraction(lev,sum_vol);
//  Print() << "Testing new sum_vol " << sum_vol << " against original sum_vol " << sum_vol_orig << std::endl;
    if (abs(sum_vol_orig - sum_vol) > 1.e-12 * sum_vol_orig)
    {
       amrex::Print() << "Original volume fraction " << sum_vol_orig << std::endl;
       amrex::Print() << "New      volume fraction " << sum_vol      << std::endl;
       // amrex::Abort("Volume fraction in domain has changed!");
    }
  }
  
  if (solve_fluid)  {
      if ( use_proj_method )	{
	  EvolveFluidProjection(lev,nstep,steady_state,dt,prev_dt,time);
      }	else {
	    EvolveFluid(lev,nstep,set_normg,dt,prev_dt,time,normg);
      }
  }
  
  // This returns the drag force on the particle
  if (solve_dem && solve_fluid)
      mfix_calc_drag_particle(lev);
  
  if (solve_dem)
      pc -> EvolveParticles( lev, nstep, dt, time, particle_ebfactory, particle_cost[lev], knapsack_weight_type, subdt_io);

}

void
mfix_level::EvolveFluid(int lev, int nstep, int set_normg,
                        Real dt, Real& prev_dt, Real time, Real normg)
{
    // Reimpose boundary conditions -- make sure to do this before we compute tau
    mfix_set_bc1(lev);

    // Calculate transport coefficients
    int calc_flag = 2;
    mfix_calc_coeffs(lev,calc_flag);

    // Calculate the stress tensor trace and cross terms for all phases.
    mfix_calc_trd_and_tau(lev);

    // Backup field variable to old
    int nghost = ep_go[lev]->nGrow();

    MultiFab::Copy(*ep_go[lev],  *ep_g[lev],  0, 0, 1, nghost);
    MultiFab::Copy( *p_go[lev],   *p_g[lev],  0, 0, 1, nghost);
    MultiFab::Copy(*ro_go[lev],  *ro_g[lev],  0, 0, 1, nghost);
    MultiFab::Copy(*rop_go[lev], *rop_g[lev], 0, 0, 1, nghost);
    MultiFab::Copy(*u_go[lev],   *u_g[lev],   0, 0, 1, nghost);
    MultiFab::Copy(*v_go[lev],   *v_g[lev],   0, 0, 1, nghost);
    MultiFab::Copy(*w_go[lev],   *w_g[lev],   0, 0, 1, nghost);

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

	    // User hooks
#ifdef _OPENMP
#pragma omp parallel
#endif
	    for (MFIter mfi(*ep_g[lev], true); mfi.isValid(); ++mfi)
		mfix_usr2();

	    // Calculate drag coefficient
	    if (solve_dem)
		mfix_calc_drag_fluid(lev);

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

