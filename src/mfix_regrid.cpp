#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_EBFabFactory.H>

void
mfix_level::Regrid (int base_lev, int nstep)
{
    BL_PROFILE_REGION_START("mfix::Regrid()");

    amrex::Print() << "In Regrid at step " << nstep << std::endl;

    if (load_balance_type == "KDTree") {
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
               RegridArrays(base_lev,grids[base_lev],dmap[base_lev]);
       }

       if (solve_fluid)
       {
           mfix_set_p0(base_lev);
           mfix_set_bc0(base_lev);
           mfix_extrap_pressure(base_lev,p0_g[base_lev]);
       }

       if (ebfactory[base_lev]) {
           const EB2::IndexSpace& index_space = EB2::IndexSpace::top();
           const EB2::Level& eb_level = index_space.getLevel(geom[base_lev]);
           ebfactory[base_lev].reset(new EBFArrayBoxFactory(eb_level, 
                                            geom[base_lev], grids[base_lev], dmap[base_lev],
                                            {m_eb_basic_grow_cells,
                                             m_eb_volume_grow_cells,
                                             m_eb_full_grow_cells}, m_eb_support_level));
       }

       if (particle_ebfactory) {
           const EB2::IndexSpace& index_space = EB2::IndexSpace::top();
           const EB2::Level& eb_level = index_space.getLevel(geom[base_lev]);
           particle_ebfactory.reset(new EBFArrayBoxFactory(eb_level,
                                           geom[base_lev], pc->ParticleBoxArray(base_lev),
                                           pc->ParticleDistributionMap(base_lev),
                                           {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                            m_eb_full_grow_cells}, m_eb_support_level));

           // eb_normals is a legacy of the old collision algorithm -> deprecated
           eb_normals   = pc->EBNormals(base_lev, particle_ebfactory.get(), dummy.get());
       }

    } else if (load_balance_type == "KnapSack") {

        amrex::Print() << "Load balancing using KnapSack " << std::endl;

        if (solve_dem)
           AMREX_ALWAYS_ASSERT(particle_cost[0] != nullptr);
        if (solve_fluid)
           AMREX_ALWAYS_ASSERT(fluid_cost[0]    != nullptr);

        if (ParallelDescriptor::NProcs() == 1) return;

        if (dual_grid) {
            AMREX_ALWAYS_ASSERT(solve_fluid);

            for (int lev = base_lev; lev <= finestLevel(); ++lev)
            {
                DistributionMapping new_fluid_dm = DistributionMapping::makeKnapSack(*fluid_cost[lev]);

                bool dm_changed = (new_fluid_dm !=  dmap[lev]);

                SetDistributionMap(lev, new_fluid_dm);

                if (dm_changed)
                    RegridArrays(lev, grids[lev], new_fluid_dm);

                fluid_cost[lev].reset(new MultiFab(grids[lev], new_fluid_dm, 1, 0));
                fluid_cost[lev]->setVal(0.0);

                if (ebfactory[lev]) {
                    const EB2::IndexSpace& index_space = EB2::IndexSpace::top();
                    const EB2::Level& eb_level = index_space.getLevel(geom[base_lev]);
                    ebfactory[lev].reset(new EBFArrayBoxFactory(eb_level,
                                                geom[lev], grids[lev], dmap[lev],
                                                {m_eb_basic_grow_cells,
                                                 m_eb_volume_grow_cells,
                                                 m_eb_full_grow_cells}, m_eb_support_level));
                }

                {
                    mfix_set_p0(lev);
                    mfix_set_bc0(lev);
                    mfix_extrap_pressure(lev,p0_g[lev]);
                }

                DistributionMapping new_particle_dm = DistributionMapping::makeKnapSack(*particle_cost[lev]);

                pc->Regrid(new_particle_dm, pc->ParticleBoxArray(lev));

                particle_cost[lev].reset(new MultiFab(pc->ParticleBoxArray(lev),
                                                      new_particle_dm, 1, 0));
                particle_cost[lev]->setVal(0.0);

                if (particle_ebfactory) {
                    const EB2::IndexSpace& index_space = EB2::IndexSpace::top();
                    const EB2::Level& eb_level = index_space.getLevel(geom[base_lev]);
                    particle_ebfactory.reset(new EBFArrayBoxFactory(eb_level,
                                                    geom[lev], pc->ParticleBoxArray(lev),
                                                    pc->ParticleDistributionMap(lev),
                                                    {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                                     m_eb_full_grow_cells}, m_eb_support_level));

                    // eb_normals is a legacy of the old collision algorithm -> deprecated
                    eb_normals   = pc->EBNormals(lev, particle_ebfactory.get(), dummy.get());
                }
            }

        } else {

            MultiFab costs(grids[base_lev], dmap[base_lev], 1, 0);
            costs.setVal(0.0);
            if (solve_dem)
               costs.plus(*particle_cost[base_lev], 0, 1, 0);
            if (solve_fluid)
                costs.plus(*fluid_cost[base_lev], 0, 1, 0);

            DistributionMapping newdm = DistributionMapping::makeKnapSack(costs);

            bool dm_changed = (newdm !=  dmap[base_lev]);

            SetDistributionMap(base_lev, newdm);

            if (solve_fluid && dm_changed)
                RegridArrays(base_lev, grids[base_lev], newdm);

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

            if (solve_dem)   pc->Regrid(dmap[base_lev], grids[base_lev]);
            if (solve_fluid) mfix_set_bc0(base_lev);

            if (ebfactory[base_lev]) {
                const EB2::IndexSpace& index_space = EB2::IndexSpace::top();
                const EB2::Level& eb_level = index_space.getLevel(geom[base_lev]);
                ebfactory[base_lev].reset(new EBFArrayBoxFactory(eb_level,
                                                 geom[base_lev], grids[base_lev], dmap[base_lev],
                                                 {m_eb_basic_grow_cells,
                                                  m_eb_volume_grow_cells,
                                                  m_eb_full_grow_cells}, m_eb_support_level));
            }

            if (particle_ebfactory) {
                const EB2::IndexSpace& index_space = EB2::IndexSpace::top();
                const EB2::Level& eb_level = index_space.getLevel(geom[base_lev]);
                particle_ebfactory.reset(new EBFArrayBoxFactory(eb_level,
                                                geom[base_lev], pc->ParticleBoxArray(base_lev),
                                                pc->ParticleDistributionMap(base_lev),
                                                {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                                 m_eb_full_grow_cells}, m_eb_support_level));

                // eb_normals is a legacy of the old collision algorithm -> deprecated
                eb_normals  = pc->EBNormals(base_lev, particle_ebfactory.get(), dummy.get());
            }
        }
    }

    // Note that this is still being done here (instead of
    // mfix_level::RegridArrays, which only acts on the fluid grid) because of
    // a dual grid: the level-set factory object regrids using the
    // ParticleDistributionMap.
    level_set->regrid(pc->ParticleBoxArray(base_lev), pc->ParticleDistributionMap(base_lev));

    BL_PROFILE_REGION_STOP("mfix::Regrid()");
}
