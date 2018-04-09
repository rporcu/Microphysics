#include <mfix_F.H>
#include <mfix_level.H>

void
mfix_level::mfix_compute_dt(int lev, Real time, Real stop_time, int steady_state, Real& dt)
{
//  if (!fixed_dt)
    {
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

       compute_new_dt ( &umax, &vmax, &wmax, &romin, &mumax, 
   		        gradp0max, geom[lev].CellSize(), &cfl, 
                        &steady_state, &time, &stop_time, &dt );

    }
}
