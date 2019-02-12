#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

// This subroutine is the driver for the whole time stepping (fluid + particles )
void
mfix::Evolve(int nstep, Real & dt, Real & prev_dt, Real time, Real stop_time)
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
          EvolveFluid(nstep,dt,time,stop_time);
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


    /****************************************************************************
     *                                                                          *
     * Evolve Particles (Using Particle MD)                                     *
     *                                                                          *
     ***************************************************************************/

    Real start_particles = ParallelDescriptor::second();

    BL_PROFILE_VAR("PARTICLES SOLVE",particlesSolve);
    if (solve_dem)
    {

        if (nlev == 1)
        {
            //___________________________________________________________________
            // Single level case: the refined level-set is stored on level 1,
            // everything else lives on level 0
            int ilev = 0;

            const MultiFab * ls_data = level_sets[1].get();
            iMultiFab ls_valid(ls_data->boxArray(), ls_data->DistributionMap(),
                               ls_data->nComp(), ls_data->nGrow());

            // Consider all values as valid
            ls_valid.setVal(1);
            ls_valid.FillBoundary(geom[ilev].periodicity());

            pc->EvolveParticles(ilev, nstep, dt, time, particle_ebfactory[ilev].get(),
                                ls_data, & ls_valid, levelset__refinement,
                                particle_cost[ilev].get(), knapsack_weight_type, subdt_io);
        }
        else
        {
            //___________________________________________________________________
            // Multi-level case: each level is treated separately

            for (int lev = 0; lev < nlev; lev ++ )
            {
                const MultiFab * ls_data = level_sets[lev].get();
                iMultiFab ls_valid(ls_data->boxArray(), ls_data->DistributionMap(),
                                   ls_data->nComp(), ls_data->nGrow());

                // Consider all values as valid
                ls_valid.setVal(1);
                ls_valid.FillBoundary(geom[lev].periodicity());

                pc->EvolveParticles(lev, nstep, dt, time, particle_ebfactory[lev].get(),
                                    ls_data, & ls_valid, 1,
                                    particle_cost[lev].get(), knapsack_weight_type, subdt_io);
            }
        }

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
