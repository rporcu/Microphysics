#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

#include <AMReX_EB_utils.H>

void
mfix::Regrid ()
{

    BL_PROFILE_REGION_START("mfix::Regrid()");

    int base_lev = 0;

    if (load_balance_type == "KDTree")  // KDTree load balancing type
    {
        if (solve_dem)
           AMREX_ALWAYS_ASSERT(particle_cost[0] == nullptr);
        if (solve_fluid)
           AMREX_ALWAYS_ASSERT(fluid_cost[0]    == nullptr);

       // This creates a new BA and new DM, re-defines the particle BA and DM to be these new ones,
       //      and calls Redistribute.  This doesn't touch the fluid grids.
       pc -> BalanceParticleLoad_KDTree ();

       if (!dual_grid)
       {
           bool ba_changed = (pc->ParticleBoxArray(base_lev)        != grids[base_lev]);
           bool dm_changed = (pc->ParticleDistributionMap(base_lev) !=  dmap[base_lev]);

           SetBoxArray       (base_lev, pc->ParticleBoxArray(base_lev));
           SetDistributionMap(base_lev, pc->ParticleDistributionMap(base_lev));

           // Since we have already allocated the fluid data we need to
           // re-define those arrays and copy from the old BoxArray to the new
           // one if the grids and/or dmap have changed.  Note that the
           // SetBoxArray and SetDistributionMap calls above have re-defined
           // grids and dmap to be the new ones.
           if (solve_fluid)
               RegridArrays(base_lev);
       }

       if (solve_fluid)
       {
           mfix_set_p0();
           mfix_set_bc0();
       }

       // This calls re-creates a proper particle_ebfactories and regrids all
       // the multifab that depend on it
       if (solve_dem)
           RegridLevelSetArray(base_lev);

    }
    else if (load_balance_type == "KnapSack" || load_balance_type == "SFC") // Knapsack and SFC
    {

        amrex::Print() << "Load balancing using " << load_balance_type << std::endl;

        if (solve_dem)
           AMREX_ALWAYS_ASSERT(particle_cost[0] != nullptr);
        if (solve_fluid)
           AMREX_ALWAYS_ASSERT(fluid_cost[0]    != nullptr);

        if (ParallelDescriptor::NProcs() == 1) return;

        if (dual_grid)  //  Beginning of dual grid regridding
        {
            AMREX_ALWAYS_ASSERT(solve_fluid);

            if (load_balance_fluid > 0)
            {
                for (int lev = base_lev; lev <= finestLevel(); ++lev)
                {

                    DistributionMapping new_fluid_dm;

                    if ( load_balance_type == "KnapSack" )
                    {
                        new_fluid_dm = DistributionMapping::makeKnapSack(*fluid_cost[lev],knapsack_nmax);
                    }
                    else
                    {
                        new_fluid_dm = DistributionMapping::makeSFC(*fluid_cost[lev],false);
                    }

                    SetDistributionMap(lev, new_fluid_dm);

                    RegridArrays(lev);

                    fluid_cost[lev].reset(new MultiFab(grids[lev], new_fluid_dm, 1, 0));
                    fluid_cost[lev]->setVal(0.0);
                }
            }

            mfix_set_p0();
            mfix_set_bc0();

            for (int lev = base_lev; lev <= finestLevel(); ++lev)
            {
                DistributionMapping new_particle_dm;

                if ( load_balance_type == "KnapSack" )
                {
                    new_particle_dm = DistributionMapping::makeKnapSack(*particle_cost[lev],knapsack_nmax);
                }
                else
                {
                    new_particle_dm = DistributionMapping::makeSFC(*particle_cost[lev],false);
                }

                pc->Regrid(new_particle_dm, pc->ParticleBoxArray(lev), lev);

                particle_cost[lev].reset(new MultiFab(pc->ParticleBoxArray(lev),
                                                      new_particle_dm, 1, 0));
                particle_cost[lev]->setVal(0.0);

                // This calls re-creates a proper particle_ebfactories
                //  and regrids all the multifab that depend on it
                if (solve_dem)
                    RegridLevelSetArray(lev);
            }

        }
        else  // Single-grid regridding
        {

            //NOTE: why are particle costs defined on fluid grids here? Or am I
            //not supposed to care because the grids are the same? (this might
            //break if there are more particle grid levels).

            MultiFab costs(grids[base_lev], dmap[base_lev], 1, 0);
            costs.setVal(0.0);
            if (solve_dem) {
                // costs.plus(* particle_cost[base_lev], 0, 1, 0);
                MultiFab particle_cost_loc(grids[base_lev], dmap[base_lev], 1, 0);
                particle_cost_loc.copy(* particle_cost[base_lev], 0, 0, 1);
                costs.plus(particle_cost_loc, 0, 1, 0);
            }
            if (solve_fluid) {
                // costs.plus(* fluid_cost[base_lev], 0, 1, 0);
                MultiFab fluid_cost_loc(grids[base_lev], dmap[base_lev], 1, 0);
                fluid_cost_loc.copy(* fluid_cost[base_lev], 0, 0, 1);
                costs.plus(fluid_cost_loc, 0, 1, 0);
            }

            DistributionMapping newdm = DistributionMapping::makeKnapSack(costs,knapsack_nmax);

            SetDistributionMap(base_lev, newdm);

            if (solve_fluid)
                RegridArrays(base_lev);

            if (solve_fluid)
            {
               fluid_cost[base_lev].reset(new MultiFab(grids[base_lev], newdm, 1, 0));
               fluid_cost[base_lev]->setVal(0.0);
            }

            if (solve_dem)
            {
               particle_cost[base_lev].reset(new MultiFab(grids[base_lev], newdm, 1, 0));
               particle_cost[base_lev]->setVal(0.0);
            }


            if (solve_dem){
                pc->Regrid(dmap[base_lev], grids[base_lev], base_lev);
            }

            if (solve_fluid) mfix_set_bc0();

            // This calls re-creates a proper particles_ebfactory and regrids
            // all the multifab that depend on it
            if (solve_dem)
                RegridLevelSetArray(base_lev);

        }
    }

    if (solve_dem)
        for (int i_lev = base_lev; i_lev < nlev; i_lev++)
        {
            // This calls re-creates a proper particle_ebfactories and regrids
            //  all the multifab that depend on it
            RegridLevelSetArray(i_lev);
        }

    if (use_amr_ls)
        for (int i_lev = 0; i_lev < pc->finestLevel(); i_lev ++)
            amr_level_set->UpdateGrids(i_lev, pc->ParticleBoxArray(i_lev),
                                       pc->ParticleDistributionMap(i_lev))

    BL_PROFILE_REGION_STOP("mfix::Regrid()");
}
