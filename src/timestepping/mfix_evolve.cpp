#include <mfix.H>
#include <mfix_fluid.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_reactions.H>

// This subroutine is the driver for the whole time stepping (fluid + particles )
void
mfix::Evolve (int nstep, Real & dt, Real & prev_dt, Real time, Real stop_time)
{
    BL_PROFILE_REGION_START("mfix::Evolve");

    Real coupling_timing(0.);
    Real drag_timing(0.);
    Real sum_vol;

    if ((m_dem.solve() || m_pic.solve()) && fluid.solve()) {

      //BL_PROFILE_REGION("CALC VOLUME FRACTION");

      Real start_coupling = ParallelDescriptor::second();
      mfix_calc_volume_fraction(sum_vol);
      //const IntVect min_epg_cell = mfixRW->mfix_print_min_epg();

      if (amrex::Math::abs(sum_vol_orig - sum_vol) > 1.e-12 * sum_vol_orig)
        {
          amrex::Print() << "Original volume fraction " << sum_vol_orig << std::endl;
          amrex::Print() << "New      volume fraction " << sum_vol      << std::endl;
        }
      coupling_timing = ParallelDescriptor::second() - start_coupling;
    }

    Real start_fluid = ParallelDescriptor::second();
    {
      BL_PROFILE_REGION("FLUID SOLVE");
      for (int lev = 0; lev <= finest_level; lev++)
      {
         if (fluid.solve())
         {
            EvolveFluid(nstep, dt, prev_dt, time, stop_time, drag_timing);
            prev_dt = dt;
         }
      }
    } // end fluid profile region

    Real end_fluid = ParallelDescriptor::second() - start_fluid - drag_timing;
    ParallelDescriptor::ReduceRealMax(end_fluid, ParallelDescriptor::IOProcessorNumber());


    // This returns the drag force on the particle
    Real new_time = time+dt;
    if ( (m_dem.solve() || m_pic.solve()) && fluid.solve()){

      //BL_PROFILE_REGION("CALC TXFR PARTICLE");

      Real start_coupling = ParallelDescriptor::second();
      for (int lev = 0; lev <= finest_level; lev++) {
        mfix_calc_txfr_particle(new_time, get_ep_g(), get_ro_g(), get_vel_g(),
                                get_T_g(), get_X_gk(), get_thermodynamic_p_g(),
                                get_gp());
      }
      coupling_timing += ParallelDescriptor::second() - start_coupling + drag_timing;

      ParallelDescriptor::ReduceRealMax(coupling_timing, ParallelDescriptor::IOProcessorNumber());
    }

    /****************************************************************************
     *                                                                          *
     * Evolve Particles (Using Particle MD)                                     *
     *                                                                          *
     ***************************************************************************/

    Real start_particles = ParallelDescriptor::second();


    int nsubsteps;

    if (m_dem.solve()) {

        BL_PROFILE_REGION("DEM PARTICLE SOLVE");

        if (finest_level == 0)
        {
            //___________________________________________________________________
            // Single level case: the refined level-set is stored on level 1,
            // everything else lives on level 0
            int ilev = 0;

            const MultiFab* ls_data = level_sets[1].get();

            if (!test_tracer_conservation)
              pc->EvolveParticles(ilev, nstep, dt, time, mfix::gravity,
                                  ebfactory[ilev].get(),
                                  particle_ebfactory[ilev].get(), ls_data,
                                  levelset_refinement,
                                  particle_cost[ilev],
                                  knapsack_weight_type, nsubsteps);
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
                                    ebfactory[lev].get(),
                                    particle_ebfactory[lev].get(), ls_data, 1,
                                    particle_cost[lev],
                                    knapsack_weight_type, nsubsteps);
            }
        }
    }

    if (m_pic.solve()) {

        BL_PROFILE_REGION("PIC PARTICLE SOLVE");
        //const IntVect min_epg_cell = mfixRW->mfix_print_min_epg();
        EvolveParcels(dt, time, mfix::gravity, levelset_refinement,
                      particle_cost, knapsack_weight_type);
    }

    Real end_particles = ParallelDescriptor::second() - start_particles;
    ParallelDescriptor::ReduceRealMax(end_particles, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor()) {
      if(fluid.solve())
        std::cout << "   Time per fluid step      " << end_fluid << std::endl;

      if(m_dem.solve())
        std::cout << "   Time per " << nsubsteps
                  << " particle steps " << end_particles << std::endl;

      if(m_pic.solve())
        std::cout << "   Time per parcel step " << end_particles << std::endl;

      if((m_dem.solve() || m_pic.solve()) && fluid.solve()) {
        std::cout << "   Coupling time per step   " << coupling_timing << std::endl;
      }
    }


    if (mfixRW->report_mass_balance and reactions.solve()) {
      mfixRW->ComputeMassProduction(dt, get_txfr_const());
    }

    BL_PROFILE_REGION_STOP("mfix::Evolve");
}
