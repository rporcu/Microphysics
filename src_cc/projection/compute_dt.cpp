#include "mfix_proj_F.H"
#include <mfix_level.H>

void
mfix_level::mfix_compute_dt(int lev, Real time, Real stop_time, int steady_state, Real& dt, int nodal_pressure)
{
    if (!fixed_dt)
    {
       // Compute dt for this time step
       Real umax  = vel_g[lev] -> norm0 (0);
       Real vmax  = vel_g[lev] -> norm0 (1);
       Real wmax  = vel_g[lev] -> norm0 (2);
       Real romin = rop_g[lev] -> min   (0);
       Real mumax =  mu_g[lev] -> max   (0);
    
       Real gradp0max[3];

       for (MFIter mfi(*vel_g[lev], true); mfi.isValid(); ++mfi) 
       {
   	   // Cell-centered tilebox
	   Box bx = mfi.tilebox();

           compute_gradp0_max (
               bx.loVect(), bx.hiVect(),
               BL_TO_FORTRAN_ANYD((*p0_g[lev])[mfi]),
               gradp0max, geom[lev].CellSize(), &nodal_pressure);
       }

       compute_new_dt ( &umax, &vmax, &wmax, &romin, &mumax, 
   		        gradp0max, geom[lev].CellSize(), &cfl, 
                        &steady_state, &time, &stop_time, &dt);

   }
}
