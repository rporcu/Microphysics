#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

// This subroutine is the driver for the whole time stepping (fluid + particles )
void
mfix::Evolve(int nstep, int steady_state, Real & dt, Real & prev_dt, Real time, Real stop_time)
{
    BL_PROFILE_REGION_START("mfix::Evolve");

    Real sum_vol;
    if (solve_dem && solve_fluid)
    {
      mfix_calc_volume_fraction(sum_vol);
      if (abs(sum_vol_orig - sum_vol) > 1.e-12 * sum_vol_orig)
      {
         amrex::Print() << "Original volume fraction " << sum_vol_orig << std::endl;
         amrex::Print() << "New      volume fraction " << sum_vol      << std::endl;
      }
    }

    Real start_fluid = ParallelDescriptor::second();

    BL_PROFILE_VAR("FLUID SOLVE",fluidSolve);
    for (int lev = 0; lev < nlev; lev++)
    {
       if (solve_fluid)
       {
          EvolveFluid(nstep,steady_state,dt,time,stop_time);
          prev_dt = dt;
       }
    }
    BL_PROFILE_VAR_STOP(fluidSolve);

    // This returns the drag force on the particle
    Real new_time = time+dt;
    if (solve_dem && solve_fluid)
       mfix_calc_drag_particle(new_time);

    Real end_fluid = ParallelDescriptor::second() - start_fluid;
    ParallelDescriptor::ReduceRealMax(end_fluid, ParallelDescriptor::IOProcessorNumber());

    Real start_particles = ParallelDescriptor::second();

    BL_PROFILE_VAR("PARTICLES SOLVE",particlesSolve);
    if (solve_dem)
    {
        if (use_amr_ls) {
            // TODO: this is a bit of hack, but it works for now
            int lev_ebfactory = 0;
            for (int lev = 0; lev <= amr_level_set->finestLevel(); lev ++){
                const MultiFab * ls_lev = amr_level_set->getLevelSet(lev);

                //TODO: valid is defunct, passing it as a legacy
                iMultiFab ls_valid(ls_lev->boxArray(), ls_lev->DistributionMap(),
                                   ls_lev->nComp(), ls_lev->nGrow());
                ls_valid.setVal(1);
                amrex::Print() << "hi lev " << lev << std::endl;
                pc->EvolveParticles(lev, nstep, dt, time,
                                    particle_ebfactory[lev_ebfactory].get(), eb_normals.get(),
                                    ls_lev, & ls_valid, 1,
                                    dummy.get(),
                                    particle_cost[lev].get(), knapsack_weight_type,
                                    subdt_io
                    );
                amrex::Print() << "bye lev " << lev << std::endl;
            }
        } else {

            for (int lev = 0; lev < nlev; lev++)
                pc->EvolveParticles(lev, nstep, dt, time,
                                    particle_ebfactory[lev].get(), eb_normals.get(),
                                    level_set->get_data(),
                                    level_set->get_valid(),
                                    level_set->get_ls_ref(),
                                    dummy.get(),
                                    particle_cost[lev].get(), knapsack_weight_type,
                                    subdt_io
                    );
        }

        //  Compute Eulerian velocities in selected regions
        for (int lev = 0; lev < nlev; lev++)
           if ( ( avg_vel_int > 0) && ( nstep % avg_vel_int == 0 ) )
             pc -> ComputeAverageVelocities ( lev,
               time,
               avg_vel_file,
               avg_region_x_w, avg_region_x_e,
               avg_region_y_s, avg_region_y_n,
               avg_region_z_b, avg_region_z_t );
    }
    BL_PROFILE_VAR_STOP(particlesSolve);

    Real end_particles = ParallelDescriptor::second() - start_particles;
    ParallelDescriptor::ReduceRealMax(end_particles, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "   Time per fluid step   " << end_fluid << std::endl;
        std::cout << "Time per particle step   " << end_particles << std::endl;
    }

    BL_PROFILE_REGION_STOP("mfix::Evolve");
}
