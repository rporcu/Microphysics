#include <mfix.H>
#include <mfix_fluid.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_reactions.H>

// This subroutine is the driver for the whole time stepping (fluid + particles )
void
mfix::Evolve (int nstep,
              Real & dt,
              Real & prev_dt,
              const Real time,
              Real stop_time)
{
    BL_PROFILE_REGION_START("mfix::Evolve");

    Real coupling_timing(0.);
    Real drag_timing(0.);
    Real sum_vol;

    if (m_dem.solve() || m_pic.solve()) {

      // sort particles by cell, this can significantly improve the locality
      if (pc->sortNow(nstep)) {
        Print() << "   Sorting particles at step " << nstep << "\n";
        pc->SortParticlesByBin(pc->getSortingBinSizes());
      }

      if ( fluid.solve() ) {
        //BL_PROFILE_REGION("CALC VOLUME FRACTION");

        Real start_coupling = ParallelDescriptor::second();
        mfix_calc_volume_fraction(time, sum_vol);
        //const IntVect min_epg_cell = m_rw->mfix_print_min_epg();

        coupling_timing = ParallelDescriptor::second() - start_coupling;
      }
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

    if (m_run_type != RunType::PIC2DEM) {

      if ( (m_dem.solve() || m_pic.solve()) && fluid.solve()) {

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
                                    knapsack_weight_type, nsubsteps,
                                    m_rw->report_mass_balance);
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
                                      knapsack_weight_type, nsubsteps,
                                      m_rw->report_mass_balance);
              }
          }
      }

      if (m_pic.solve()) {

          BL_PROFILE_REGION("PIC PARTICLE SOLVE");
          //const IntVect min_epg_cell = m_rw->mfix_print_min_epg();
          EvolveParcels(dt, time, mfix::gravity, levelset_refinement,
                        particle_cost, knapsack_weight_type,
                        m_rw->report_mass_balance);
      }


      if (m_dem.solve() || m_pic.solve()) {
        if (pc->UseConstraint()) {
          for (int lev = 0; lev <= finest_level; lev ++ ) {
            pc->ImposeMean(lev);
          }
        }
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


      if (m_rw->report_mass_balance && reactions.solve()) {
        m_rw->ComputeMassProduction(dt, get_txfr_const());
      }
    }

    if (m_run_type == RunType::PIC2DEM) {

      for (int lev(0); lev < nlev; ++lev) {

        std::swap(m_leveldata[lev]->ro_g, m_leveldata[lev]->ro_go);
        std::swap(m_leveldata[lev]->trac, m_leveldata[lev]->trac_o);
        std::swap(m_leveldata[lev]->vel_g, m_leveldata[lev]->vel_go);

        if (fluid.solve_enthalpy()) {
          std::swap(m_leveldata[lev]->T_g, m_leveldata[lev]->T_go);
        }

        if (fluid.solve_species())
          std::swap(m_leveldata[lev]->X_gk, m_leveldata[lev]->X_gko);
      }
    }

    BL_PROFILE_REGION_STOP("mfix::Evolve");
}
