#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

// This subroutine is the driver for the whole time stepping (fluid + particles )
void
mfix_level::Evolve(int lev, int nstep, int set_normg, int steady_state,  Real& dt, Real & prev_dt,
                   Real time, Real stop_time, Real normg)
{
    BL_PROFILE_REGION_START("mfix::Evolve");

    AMREX_ALWAYS_ASSERT(lev == 0);

    Real sum_vol;
    if (solve_dem && solve_fluid)
    {
        mfix_calc_volume_fraction(lev,sum_vol);
        //Print() << "Testing new sum_vol " << sum_vol << " against original sum_vol " << sum_vol_orig << std::endl;
        if (abs(sum_vol_orig - sum_vol) > 1.e-12 * sum_vol_orig)
        {
            amrex::Print() << "Original volume fraction " << sum_vol_orig << std::endl;
            amrex::Print() << "New      volume fraction " << sum_vol      << std::endl;
        }
    }

    if (solve_fluid)  {
        if ( use_proj_method ) {
            EvolveFluidProjection(lev,nstep,steady_state,dt,time,stop_time);
            prev_dt = dt;
        } else {
            EvolveFluidSimple(lev,nstep,set_normg,dt,prev_dt,time,normg);
        }
    }

    // This returns the drag force on the particle
    if (solve_dem && solve_fluid)
        mfix_calc_drag_particle(lev);

    if (solve_dem)
    {
        pc->EvolveParticles(lev, nstep, dt, time,
                            particle_ebfactory[lev].get(), eb_normals.get(),
                            level_set->get_data(), level_set->get_valid(), level_set->get_ls_ref(),
                            dummy.get(), particle_cost[lev].get(), knapsack_weight_type, subdt_io
	    );

        //  Compute Eulerian velocities in selected regions
        if ( ( avg_vel_int > 0) && ( nstep % avg_vel_int == 0 ) )
   	    pc -> ComputeAverageVelocities ( lev,
	 				     time,
					     avg_vel_file,
					     avg_region_x_w, avg_region_x_e,
					     avg_region_y_s, avg_region_y_n,
					     avg_region_z_b, avg_region_z_t );
    }
    BL_PROFILE_REGION_STOP("mfix::Evolve");
}
 
