#include <mfix.H>
#include <mfix_fluid_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_reactions_parms.H>

// This subroutine is the driver for the whole time stepping (fluid + particles )
void
mfix::Evolve (int nstep, Real & dt, Real & prev_dt, Real time, Real stop_time)
{
    BL_PROFILE_REGION_START("mfix::Evolve");

    Real coupling_timing(0.);
    Real drag_timing(0.);
    Real sum_vol;

    if ((DEM::solve || PIC::solve) && fluid.solve)
    {
      Real start_coupling = ParallelDescriptor::second();
      mfix_calc_volume_fraction(sum_vol);
      //const IntVect min_epg_cell = mfix_print_min_epg();

      if (amrex::Math::abs(sum_vol_orig - sum_vol) > 1.e-12 * sum_vol_orig)
        {
          amrex::Print() << "Original volume fraction " << sum_vol_orig << std::endl;
          amrex::Print() << "New      volume fraction " << sum_vol      << std::endl;
        }
      coupling_timing = ParallelDescriptor::second() - start_coupling;
    }

    Real start_fluid = ParallelDescriptor::second();
    BL_PROFILE_VAR("FLUID SOLVE",fluidSolve);
    for (int lev = 0; lev <= finest_level; lev++)
    {
       if (fluid.solve)
       {
          EvolveFluid(nstep, dt, prev_dt, time, stop_time, drag_timing);
          prev_dt = dt;
       }
    }
    BL_PROFILE_VAR_STOP(fluidSolve);

    Real end_fluid = ParallelDescriptor::second() - start_fluid - drag_timing;
    ParallelDescriptor::ReduceRealMax(end_fluid, ParallelDescriptor::IOProcessorNumber());


    // This returns the drag force on the particle
    Real new_time = time+dt;
    if ( (DEM::solve || PIC::solve) && fluid.solve){
      Real start_coupling = ParallelDescriptor::second();

      mfix_calc_txfr_particle(new_time, get_vel_g(), get_gp(), get_T_g());

      if (reactions.solve)
        mfix_calc_chem_txfr(get_chem_txfr(), get_ep_g(), get_ro_g(), get_vel_g(),
                            get_T_g(), get_X_gk(), new_time);

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
        if (finest_level == 0)
        {
            //___________________________________________________________________
            // Single level case: the refined level-set is stored on level 1,
            // everything else lives on level 0
            int ilev = 0;

            const MultiFab* ls_data = level_sets[1].get();

            if (!test_tracer_conservation)
              pc->EvolveParticles(ilev, nstep, dt, time, mfix::gravity,
                                  particle_ebfactory[ilev].get(), ls_data,
                                  levelset_refinement,
                                  particle_cost[ilev],
                                  knapsack_weight_type, nsubsteps,
                                  advect_enthalpy, enthalpy_source,
                                  update_mass, update_momentum,
                                  update_enthalpy);
        }
        else
        {
            //___________________________________________________________________
            // Multi-level case: each level is treated separately

            for (int lev = 0; lev <= finest_level; lev ++ )
            {
                const MultiFab* ls_data = level_sets[lev].get();

                if (!test_tracer_conservation)
                pc->EvolveParticles(lev, nstep, dt, time, mfix::gravity,
                                    particle_ebfactory[lev].get(), ls_data, 1,
                                    particle_cost[lev],
                                    knapsack_weight_type, nsubsteps,
                                    advect_enthalpy, enthalpy_source,
                                    update_mass, update_momentum,
                                    update_enthalpy);
            }
        }
    }

    if (PIC::solve) {
        //const IntVect min_epg_cell = mfix_print_min_epg();
        EvolveParcels(dt, time, mfix::gravity, levelset_refinement,
                      particle_cost, knapsack_weight_type, advect_enthalpy,
                      enthalpy_source);
    }

    BL_PROFILE_VAR_STOP(particlesSolve);

    Real end_particles = ParallelDescriptor::second() - start_particles;
    ParallelDescriptor::ReduceRealMax(end_particles, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor()) {
      if(fluid.solve)
        std::cout << "   Time per fluid step      " << end_fluid << std::endl;

      if(DEM::solve)
        std::cout << "   Time per " << nsubsteps
                  << " particle steps " << end_particles << std::endl;

      if(PIC::solve)
        std::cout << "   Time per parcel step " << end_particles << std::endl;

      if((DEM::solve || PIC::solve) && fluid.solve) {
        std::cout << "   Coupling time per step   " << coupling_timing << std::endl;
      }
    }

    BL_PROFILE_REGION_STOP("mfix::Evolve");
}
