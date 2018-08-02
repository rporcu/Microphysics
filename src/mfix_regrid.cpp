#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

void
mfix_level::Regrid (int lev, int nstep)
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
           bool ba_changed = (pc->ParticleBoxArray(lev)        != grids[lev]);
           bool dm_changed = (pc->ParticleDistributionMap(lev) !=  dmap[lev]);

           SetBoxArray       (lev, pc->ParticleBoxArray(lev));
           SetDistributionMap(lev, pc->ParticleDistributionMap(lev));

           // Since we have already allocated the fluid data we need to
           // re-define those arrays and copy from the old BoxArray to the new
           // one if the grids and/or dmap have changed.  Note that the
           // SetBoxArray and SetDistributionMap calls above have re-defined
           // grids and dmap to be the new ones.
           if (solve_fluid && (ba_changed || dm_changed) )
               RegridArrays(lev,grids[lev],dmap[lev]);
       }

       if (solve_fluid)
       {
           mfix_set_p0(lev);
           mfix_set_bc0(lev);
           mfix_extrap_pressure(lev,p0_g[lev]);
       }

       if (ebfactory) {
           ebfactory.reset(new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                                                  {m_eb_basic_grow_cells,
                                                   m_eb_volume_grow_cells,
                                                   m_eb_full_grow_cells}, m_eb_support_level));
       }

       if (particle_ebfactory) {
           particle_ebfactory.reset(new EBFArrayBoxFactory(geom[lev], pc->ParticleBoxArray(lev),
                                                           pc->ParticleDistributionMap(lev),
                                                           {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                                            m_eb_full_grow_cells}, m_eb_support_level));

           // eb_normals is a legacy of the old collision algorithm -> deprecated
           eb_normals   = pc->EBNormals(lev, particle_ebfactory.get(), dummy.get());
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

            for (int lev = 0; lev <= finestLevel(); ++lev)
            {
                DistributionMapping new_fluid_dm = DistributionMapping::makeKnapSack(*fluid_cost[lev]);

                bool dm_changed = (new_fluid_dm !=  dmap[lev]);

                SetDistributionMap(lev, new_fluid_dm);

                if (dm_changed)
                    RegridArrays(lev, grids[lev], new_fluid_dm);

                fluid_cost[lev].reset(new MultiFab(grids[lev], new_fluid_dm, 1, 0));
                fluid_cost[lev]->setVal(0.0);

                if (ebfactory) {
                    ebfactory.reset(new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
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
                    particle_ebfactory.reset(new EBFArrayBoxFactory(geom[lev], pc->ParticleBoxArray(lev),
                                                                    pc->ParticleDistributionMap(lev),
                                                                    {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                                                     m_eb_full_grow_cells}, m_eb_support_level));

                    // eb_normals is a legacy of the old collision algorithm -> deprecated
                    eb_normals   = pc->EBNormals(lev, particle_ebfactory.get(), dummy.get());
                }
            }

        } else {

            MultiFab costs(grids[lev], dmap[lev], 1, 0);
            costs.setVal(0.0);
            if (solve_dem)
               costs.plus(*particle_cost[lev], 0, 1, 0);
            if (solve_fluid)
                costs.plus(*fluid_cost[lev], 0, 1, 0);

            DistributionMapping newdm = DistributionMapping::makeKnapSack(costs);

            bool dm_changed = (newdm !=  dmap[lev]);

            SetDistributionMap(lev, newdm);

            if (solve_fluid && dm_changed)
                RegridArrays(lev, grids[lev], newdm);

            if (solve_fluid)
            {
               fluid_cost[lev].reset(new MultiFab(grids[lev], newdm, 1, 0));
               fluid_cost[lev]->setVal(0.0);
            }

            if (solve_dem)
            {
               particle_cost[lev].reset(new MultiFab(grids[lev], newdm, 1, 0));
               particle_cost[lev]->setVal(0.0);
            }

            if (solve_dem)   pc->Regrid(dmap[lev], grids[lev]);
            if (solve_fluid) mfix_set_bc0(lev);

            if (ebfactory) {
                ebfactory.reset(new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                                                       {m_eb_basic_grow_cells,
                                                        m_eb_volume_grow_cells,
                                                        m_eb_full_grow_cells}, m_eb_support_level));
            }

            if (particle_ebfactory) {
                particle_ebfactory.reset(new EBFArrayBoxFactory(geom[lev], pc->ParticleBoxArray(lev),
                                                                pc->ParticleDistributionMap(lev),
                                                                {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                                                 m_eb_full_grow_cells}, m_eb_support_level));

                // eb_normals is a legacy of the old collision algorithm -> deprecated
                eb_normals  = pc->EBNormals(lev, particle_ebfactory.get(), dummy.get());
            }
        }
    }

    // Note that this is still being done here (instead of
    // mfix_level::RegridArrays, which only acts on the fluid grid) because of
    // a dual grid: the level-set factory object regrids using the
    // ParticleDistributionMap.
    level_set->regrid(pc->ParticleBoxArray(lev), pc->ParticleDistributionMap(lev));

    BL_PROFILE_REGION_STOP("mfix::Regrid()");
}

void
mfix_level::RegridOnStart (int lev)
{

  if (load_balance_type == "FixedSize" || load_balance_type == "KnapSack") {

    ParmParse pp("particles");
    pp.query("max_grid_size_x", particle_max_grid_size_x);
    pp.query("max_grid_size_y", particle_max_grid_size_y);
    pp.query("max_grid_size_z", particle_max_grid_size_z);

    Regrid(lev,0);

  }
}
