#include "mfix_proj_F.H"
#include <mfix.H>

#include <cmath>
#include <limits>

void
mfix::mfix_compute_dt(Real time, Real stop_time, Real& dt)
{
    // dt is always computed even when fixed_dt is set, 
    // so we can issue a warning if the value of fixed dt does not satisfy the CFL condition.
    Real dt_new = dt;

    /* 
       Compute new dt by using the formula derived in
       "A Boundary Condition Capturing Method for Multiphase Incompressible Flow"
       by Kang et al. (JCP).
      
       dt/2 * ( C+V + sqrt( (C+V)**2 + 4Fx/dx + 4Fy/dy + 4Fz/dz )
     
      where
      
      C = max(|U|)/dx + max(|V|)/dy + max(|W|)/dz    --> Convection
      
      V = 2 * max(mu/ro) * (1/dx^2 + 1/dy^2 +1/dz^2) --> Diffusion
      
      Fx, Fy, Fz = net acceleration due to external forces
    */

    Real umax = -1.e20;
    Real vmax = -1.e20;
    Real wmax = -1.e20;
    Real romin = 1.e20;
    Real mumax = 0.0;

    Real ope = 1.0 + 1.e-8;
    
    Real gp0max[3];
    gp0max[0]  = std::abs(gp0[0]);
    gp0max[1]  = std::abs(gp0[1]);
    gp0max[2]  = std::abs(gp0[2]);

    Real c_cfl  = 0.0;
    Real v_cfl  = 0.0;
    Real f_cfl  = 0.0;
    Real old_dt = dt;

    for (int lev = 0; lev < nlev; lev++)
    {
       // These take the min over un-covered cells
       umax  = amrex::max(umax,mfix_norm0 ( vel_g, lev, 0 ));
       vmax  = amrex::max(vmax,mfix_norm0 ( vel_g, lev, 1 ));
       wmax  = amrex::max(wmax,mfix_norm0 ( vel_g, lev, 2 ));
       mumax = amrex::max(mumax,mfix_norm0( mu_g,  lev, 0 ));

       // This takes the min of (ro_g * ep_g) over un-covered cells
       romin = amrex::min(romin, mfix_norm0( ro_g, ep_g, lev, 0, 0 ));
    }

    const Real* dx = geom[finest_level].CellSize();

    Real odx    = 1.0 / dx[0];
    Real ody    = 1.0 / dx[1];
    Real odz    = 1.0 / dx[2];

    // Convection
    c_cfl = std::max(std::max(umax*odx,vmax*ody), wmax*odz);

    // Viscous
    v_cfl = 2.0 * ( mumax / romin ) * ( odx*odx + ody*ody + odz*odz );

    // Gravity and/or gradient of p0
    f_cfl = std::abs(gravity[0]-gp0max[0]) * odx + 
            std::abs(gravity[1]-gp0max[1]) * ody + 
            std::abs(gravity[2]-gp0max[2]) * odz;

    // Put all together
    Real tmp_cfl = (c_cfl+v_cfl);
    Real tmp     = tmp_cfl + std::sqrt( tmp_cfl*tmp_cfl + 4.0*f_cfl );
    dt_new  = cfl * 2.0 / tmp;

    // Protect against tmp very small
    // This may happen, for example, when the initial velocity field
    // is zero for an inviscid flow with no external forcing
    Real eps = std::numeric_limits<Real>::epsilon();
    if ( tmp <= eps ) dt_new = 0.5 * old_dt;

    // Don't let the timestep grow by more than 1% per step.
    dt_new = std::min ( dt_new, 1.01*old_dt );

    // Don't overshoot the final time if not running to steady state
    if (steady_state == 0 && stop_time > 0.) 
    {
       if (time+dt_new > stop_time) dt_new = stop_time - time;
    }

    // dt_new is the step calculated with a cfl contraint; dt is the value set by fixed_dt
    // When the test was on dt > dt_new, there were cases where they were effectively equal 
    //   but (dt > dt_new) was being set to true due to precision issues.
    if ( fixed_dt )
    {
        if ( dt > dt_new*ope && cfl > 0)
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
        dt = amrex::min( dt_new, dt_max );

        if ( dt < dt_min )
            amrex::Abort ("Current dt is smaller than dt_min");

    }
}
