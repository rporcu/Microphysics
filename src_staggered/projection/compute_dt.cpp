#include <mfix_F.H>
#include <mfix.H>

void
mfix::mfix_compute_dt(int lev, Real time, Real stop_time, int steady_state, Real& dt)
{

    // DT is always computed even for fixed dt, so we can
    // issue a warning if fixed dt does not satisfy CFL condition.
    Real dt_new = dt;
    
    // Compute dt for this time step
    Real umax  =   u_g[lev] -> norm0 ();
    Real vmax  =   v_g[lev] -> norm0 ();
    Real wmax  =   w_g[lev] -> norm0 ();
    Real romin = rop_g[lev] -> min   (0);
    Real mumax =  mu_g[lev] -> max   (0);
    
    Real gradp0max[3];
    
    for (MFIter mfi(*p0_g[lev], true); mfi.isValid(); ++mfi) 
    {
	// Tilebox
	Box bx = mfi.tilebox();
	
	compute_gradp0_max (
               bx.loVect(), bx.hiVect(),
               BL_TO_FORTRAN_ANYD((*p0_g[lev])[mfi]),
               gradp0max, geom[lev].CellSize());
    }

    ParallelDescriptor::ReduceRealMax(gradp0max[0]);
    ParallelDescriptor::ReduceRealMax(gradp0max[1]);
    ParallelDescriptor::ReduceRealMax(gradp0max[2]);

    compute_new_dt ( &umax, &vmax, &wmax, &romin, &mumax, 
		     gradp0max, geom[lev].CellSize(), &cfl, 
		     &steady_state, &time, &stop_time, &dt_new );

    if ( fixed_dt )
    {
	if ( dt_new < dt )
	    amrex::Print() << "WARNING: current dt does not satisfy CFL condition: " 
			   << "min dt required = " << dt_new
			   << " <  current dt = "  << dt
			   << endl;
    }
    else
    {	    
	dt = min( dt_new, dt_max );

	if ( dt < dt_min )
	    amrex::Abort ("Current dt is smaller than dt_min");
	
    }
}
