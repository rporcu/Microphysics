#include <mfix.H>
#include <MFIX_FLUID_Parms.H>
#include <MFIX_DEM_Parms.H>

// This subroutine is the driver for the whole time stepping (fluid + particles )
void
mfix::Evolve (int nstep, Real & dt, Real & prev_dt, Real time, Real stop_time)
{
    BL_PROFILE_REGION_START("mfix::Evolve");

    Real coupling_timing;
    Real sum_vol;
    if (DEM::solve and FLUID::solve)
    {
      Real start_coupling = ParallelDescriptor::second();
      mfix_calc_volume_fraction(sum_vol);

      if (abs(sum_vol_orig - sum_vol) > 1.e-12 * sum_vol_orig)
        {
          amrex::Print() << "Original volume fraction " << sum_vol_orig << std::endl;
          amrex::Print() << "New      volume fraction " << sum_vol      << std::endl;
        }
      coupling_timing = ParallelDescriptor::second() - start_coupling;
    }

    Real start_fluid = ParallelDescriptor::second();
    Real drag_timing = 0.;

    BL_PROFILE_VAR("FLUID SOLVE",fluidSolve);
    for (int lev = 0; lev < nlev; lev++)
    {
       if (FLUID::solve)
       {
         EvolveFluid(nstep,dt,time,stop_time, drag_timing);
          prev_dt = dt;
       }
    }
    BL_PROFILE_VAR_STOP(fluidSolve);

    Real end_fluid = ParallelDescriptor::second() - start_fluid - drag_timing;
    ParallelDescriptor::ReduceRealMax(end_fluid, ParallelDescriptor::IOProcessorNumber());


    // This returns the drag force on the particle
    Real new_time = time+dt;
    if (DEM::solve and FLUID::solve){
      Real start_coupling = ParallelDescriptor::second();
      mfix_calc_drag_particle(new_time);
      coupling_timing += ParallelDescriptor::second() - start_coupling + drag_timing;
      ParallelDescriptor::ReduceRealMax(coupling_timing, ParallelDescriptor::IOProcessorNumber());
    }

    /****************************************************************************
     *                                                                          *
     * Evolve Particles (Using Particle MD)                                     *
     *                                                                          *
     ***************************************************************************/

    Real start_particles = ParallelDescriptor::second();

    BL_PROFILE_VAR("PARTICLES SOLVE", particlesSolve);

    int nsubsteps;

    if (DEM::solve)
    {
        if (nlev == 1)
        {
            //___________________________________________________________________
            // Single level case: the refined level-set is stored on level 1,
            // everything else lives on level 0
            int ilev = 0;

            const MultiFab* ls_data = level_sets[1];

            if (!test_tracer_conservation)
              pc->EvolveParticles(ilev, nstep, dt, time, mfix::gravity,
                                  particle_ebfactory[ilev], ls_data,
                                  levelset_refinement, particle_cost[ilev],
                                  knapsack_weight_type, nsubsteps);
        }
        else
        {
            //___________________________________________________________________
            // Multi-level case: each level is treated separately

            for (int lev = 0; lev < nlev; lev ++ )
            {
                const MultiFab* ls_data = level_sets[lev];

                if (!test_tracer_conservation)
                pc->EvolveParticles(lev, nstep, dt, time, mfix::gravity,
                                    particle_ebfactory[lev], ls_data, 1,
                                    particle_cost[lev], knapsack_weight_type,
                                    nsubsteps);
            }
        }
    }

    BL_PROFILE_VAR_STOP(particlesSolve);

    Real end_particles = ParallelDescriptor::second() - start_particles;
    ParallelDescriptor::ReduceRealMax(end_particles, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor()) {
      if(FLUID::solve) std::cout << "   Time per fluid step      " << end_fluid << std::endl;
      if(DEM::solve  ) std::cout << "   Time per " << nsubsteps << " particle steps " << end_particles << std::endl;
      if(DEM::solve && FLUID::solve) std::cout << "   Coupling time per step   " << coupling_timing << std::endl;
    }

    BL_PROFILE_REGION_STOP("mfix::Evolve");
}
