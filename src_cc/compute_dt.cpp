#include "mfix_proj_F.H"
#include <mfix.H>

void
mfix::mfix_compute_dt(Real time, Real stop_time, int steady_state, Real& dt)
{
    // DT is always computed even for fixed dt, so we can
    // issue a warning if fixed dt does not satisfy CFL condition.
    Real dt_new = dt;

    Real umax = -1.e20;
    Real vmax = -1.e20;
    Real wmax = -1.e20;
    Real romin = 1.e20;
    Real mumax = 0.0;
    
    // We only compute gp0max on the coarset level because it is the same at all levels
    Real gp0max[3];
    gp0max[0]  = mfix_norm0 ( gp0, 0, 0 );
    gp0max[1]  = mfix_norm0 ( gp0, 0, 1 );
    gp0max[2]  = mfix_norm0 ( gp0, 0, 2 );

    for (int lev = 0; lev < nlev; lev++)
    {
       // Compute dt for this time step
       umax  = max(umax,mfix_norm0 ( vel_g, lev, 0 ));
       vmax  = max(vmax,mfix_norm0 ( vel_g, lev, 1 ));
       wmax  = max(wmax,mfix_norm0 ( vel_g, lev, 2 ));
       romin = min(romin,mfix_norm0( rop_g, lev, 0 ));
       mumax = max(mumax,mfix_norm0( mu_g,  lev, 0 ));
    }

    compute_new_dt ( &umax, &vmax, &wmax, &romin, &mumax, 
                     gp0max, geom[finest_level].CellSize(), &cfl, 
                     &steady_state, &time, &stop_time, &dt_new);

    if ( fixed_dt )
    {
        if ( dt_new < dt && cfl > 0)
        {
            amrex::Print() << "WARNING: fixed dt does not satisfy CFL condition: "
                           << " fixed dt = "  << dt
                           << " > dt based on cfl = " << dt_new
                           << endl;
            amrex::Abort ("Fixed dt is too large for fluid solve");
        }
    }
    else
    {
        dt = min( dt_new, dt_max );

        if ( dt < dt_min )
            amrex::Abort ("Current dt is smaller than dt_min");

    }
}
