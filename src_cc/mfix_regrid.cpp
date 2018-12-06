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
           if (solve_fluid && (ba_changed || dm_changed) )
               RegridArrays(base_lev);
       }

       if (solve_fluid)
       {
           mfix_set_p0();
           mfix_set_bc0();
           mfix_extrap_pressure(base_lev,p0_g[base_lev]);
       }


       if (particle_ebfactory[base_lev]) {

           particle_ebfactory[base_lev].reset(
               new EBFArrayBoxFactory(* eb_level_particles, geom[base_lev],
                                      pc->ParticleBoxArray(base_lev),
                                      pc->ParticleDistributionMap(base_lev),
                                      {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                       m_eb_full_grow_cells}, m_eb_support_level )
               );

           // eb_normals is a legacy of the old collision algorithm -> deprecated
           // eb_normals[base_lev] = pc->EBNormals(
           //     base_lev, particle_ebfactory[base_lev].get(),dummy[base_lev].get());

           // eb_normals[base_lev]->define(grids[base_lev], dmap[base_lev], 3, 2,
           //                              MFInfo(), * particle_ebfactory[base_lev]);
           eb_normals[base_lev]->define(pc->ParticleBoxArray(base_lev),
                                        pc->ParticleDistributionMap(base_lev),
                                        3, 2, MFInfo(), * particle_ebfactory[base_lev]);

           amrex::FillEBNormals( * eb_normals[base_lev], * particle_ebfactory[base_lev],
                                 geom[base_lev]);

           dummy[base_lev]->define(pc->ParticleBoxArray(base_lev),
                                   pc->ParticleDistributionMap(base_lev),
                                   3, 2, MFInfo(), * particle_ebfactory[base_lev]);

       }

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
                        new_fluid_dm = DistributionMapping::makeKnapSack(*fluid_cost[lev]);
                    }
                    else
                    {
                        new_fluid_dm = DistributionMapping::makeSFC(*fluid_cost[lev],false);
                    }

                    bool dm_changed = (new_fluid_dm !=  dmap[lev]);

                    SetDistributionMap(lev, new_fluid_dm);

                    if (dm_changed)
                        RegridArrays(lev);

                    fluid_cost[lev].reset(new MultiFab(grids[lev], new_fluid_dm, 1, 0));
                    fluid_cost[lev]->setVal(0.0);
                }
            }

            mfix_set_p0();
            mfix_set_bc0();

            for (int lev = base_lev; lev <= finestLevel(); ++lev)
              mfix_extrap_pressure(lev,p0_g[lev]);

            for (int lev = base_lev; lev <= finestLevel(); ++lev)
            {
                DistributionMapping new_particle_dm;

                if ( load_balance_type == "KnapSack" )
                {
                    new_particle_dm = DistributionMapping::makeKnapSack(*particle_cost[lev]);
                }
                else
                {
                    new_particle_dm = DistributionMapping::makeSFC(*particle_cost[lev],false);
                }

                pc->Regrid(new_particle_dm, pc->ParticleBoxArray(lev));

                particle_cost[lev].reset(new MultiFab(pc->ParticleBoxArray(lev),
                                                      new_particle_dm, 1, 0));
                particle_cost[lev]->setVal(0.0);

                if (particle_ebfactory[lev]) {

                    particle_ebfactory[lev].reset(
                        new EBFArrayBoxFactory(* eb_level_particles, geom[lev],
                                               pc->ParticleBoxArray(lev),
                                               pc->ParticleDistributionMap(lev),
                                               {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                                m_eb_full_grow_cells}, m_eb_support_level)
                        );

                    // eb_normals is a legacy of the old collision algorithm -> deprecated
                    // eb_normals[lev] = pc->EBNormals(
                    //     lev, particle_ebfactory[lev].get(), dummy[lev].get());

                    // eb_normals[lev]->define(grids[lev], dmap[lev], 3, 2,
                    //                         MFInfo(), *particle_ebfactory[lev]);
                    eb_normals[lev]->define(pc->ParticleBoxArray(lev),
                                                 pc->ParticleDistributionMap(lev),
                                                 3, 2, MFInfo(), * particle_ebfactory[lev]);

                    amrex::FillEBNormals( * eb_normals[lev], * particle_ebfactory[lev], geom[lev]);

                    dummy[lev]->define(pc->ParticleBoxArray(lev),
                                            pc->ParticleDistributionMap(lev),
                                            3, 2, MFInfo(), * particle_ebfactory[lev]);
                }
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

            DistributionMapping newdm = DistributionMapping::makeKnapSack(costs);

            bool dm_changed = (newdm !=  dmap[base_lev]);

            SetDistributionMap(base_lev, newdm);

            if (solve_fluid && dm_changed)
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
                pc->Regrid(dmap[base_lev], grids[base_lev]);
            }

            if (solve_fluid) mfix_set_bc0();

            if (particle_ebfactory[base_lev]) {

                particle_ebfactory[base_lev].reset(
                    new EBFArrayBoxFactory(* eb_level_particles, geom[base_lev],
                                           pc->ParticleBoxArray(base_lev),
                                           pc->ParticleDistributionMap(base_lev),
                                           {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                            m_eb_full_grow_cells}, m_eb_support_level)
                    );

                // eb_normals is a legacy of the old collision algorithm -> deprecated
                // eb_normals[base_lev] = pc->EBNormals(
                //     base_lev, particle_ebfactory[base_lev].get(), dummy[base_lev].get());

                // eb_normals[base_lev]->define(grids[base_lev], dmap[base_lev], 3, 2,
                //                              MFInfo(), * particle_ebfactory[base_lev]);
                eb_normals[base_lev]->define(pc->ParticleBoxArray(base_lev),
                                             pc->ParticleDistributionMap(base_lev),
                                             3, 2, MFInfo(), * particle_ebfactory[base_lev]);

                amrex::FillEBNormals( * eb_normals[base_lev], * particle_ebfactory[base_lev],
                                      geom[base_lev]);

                dummy[base_lev]->define(pc->ParticleBoxArray(base_lev),
                                        pc->ParticleDistributionMap(base_lev),
                                        3, 2, MFInfo(), * particle_ebfactory[base_lev]);
            }
        }
    }


    // Note that this is still being done here (instead of mfix::RegridArrays,
    // which only acts on the fluid grid) because of the dual grid.

    if (solve_dem)
        level_set->regrid(pc->ParticleBoxArray(base_lev),
                          pc->ParticleDistributionMap(base_lev));


    if (use_amr_ls)
        for (int i_lev = 0; i_lev < pc->finestLevel(); i_lev ++)
            amr_level_set->UpdateGrids(i_lev, pc->ParticleBoxArray(i_lev),
                                       pc->ParticleDistributionMap(i_lev))

    BL_PROFILE_REGION_STOP("mfix::Regrid()");
}
