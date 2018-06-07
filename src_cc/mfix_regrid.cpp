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

    if (load_balance_type == "KDTree")
    {
        if (solve_dem)
           AMREX_ALWAYS_ASSERT(particle_cost[0] == nullptr);
        if (solve_fluid)
           AMREX_ALWAYS_ASSERT(fluid_cost[0]    == nullptr);

       // This creates a new BA and new DM, re-defines the particle BA and DM to be these new ones,
       //      and calls Redistribute.  This doesn't touch the fluid grids.
       pc -> BalanceParticleLoad_KDTree ();

       if (dual_grid == 0)
       {
           SetBoxArray(lev, pc->ParticleBoxArray(lev));
           SetDistributionMap(lev, pc->ParticleDistributionMap(lev));
           
           // Since we have already allocated the fluid data we need to re-define those arrays
           //   and copy from the old BoxArray to the new one.  Note that the SetBoxArray and
           //   SetDistributionMap calls above have re-defined grids and dmap to be the new ones.
           if (solve_fluid) 
               RegridArrays(lev,grids[lev],dmap[lev]);
       }
       
       if (solve_fluid) 
           mfix_set_bc0(lev);
       
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

           // eb_normals is a legacy of the old collision algorithm -> depricated
           eb_normals   = pc -> EBNormals(lev, particle_ebfactory.get(), dummy.get());
       }
    }
    else if (load_balance_type == "KnapSack") {
        
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
                RegridArrays(lev, grids[lev], new_fluid_dm);
                SetDistributionMap(lev, new_fluid_dm);
                fluid_cost[lev].reset(new MultiFab(grids[lev], new_fluid_dm, 1, 0));
                fluid_cost[lev]->setVal(0.0);

                if (ebfactory) {
                    ebfactory.reset(new EBFArrayBoxFactory(geom[lev], grids[lev], dmap[lev],
                                                           {m_eb_basic_grow_cells,
                                                            m_eb_volume_grow_cells,
                                                            m_eb_full_grow_cells}, m_eb_support_level));                    
                }

                mfix_set_bc0(lev);

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

                    // eb_normals is a legacy of the old collision algorithm -> depricated
                    eb_normals   = pc -> EBNormals(lev, particle_ebfactory.get(), dummy.get());
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
            if (solve_fluid) 
                RegridArrays(lev, grids[lev], newdm);
            SetDistributionMap(lev, newdm);

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

                // eb_normals is a legacy of the old collision algorithm -> depricated
                eb_normals   = pc -> EBNormals(lev, particle_ebfactory.get(), dummy.get());
            }
        }
    }

    // Note: this might not be necessary anymore if the level-set data is
    // managed by mfix_level
    level_set->regrid();

    BL_PROFILE_REGION_STOP("mfix::Regrid()");
}

void
mfix_level::RegridOnRestart (int lev)
{
    amrex::Print() << "In RegridOnRestart " << std::endl;

    if (load_balance_type == "FixedSize" || load_balance_type == "KnapSack")
    {
       // We hold on to the old_ba so that we can test the new BoxArray against it
       //   to see if the grids have changed
       BoxArray old_ba(grids[lev]);

       // This creates a new BoxArray (based on the new max_grid_size)
       const BoxArray& ba = MakeBaseGrids();

       // This creates the associated Distribution Mapping
       DistributionMapping dm(ba, ParallelDescriptor::NProcs());

       // This sets grids[lev] = ba 
       SetBoxArray(lev, ba);

       bool ba_changed = (old_ba != grids[lev]);
       bool dm_changed = (dm     !=  dmap[lev]);

       // This sets dmap[lev] = dm
       SetDistributionMap(lev, dm);

       // If the grids have changed, we need to re-define those arrays
       // and copy from the old BoxArray to the new one since we have already read in the old fluid data
       if ( solve_fluid )
       {
         if ( ba_changed || dm_changed )
            RegridArrays(lev,grids[lev],dmap[lev]);
       }

       if (solve_fluid)
          mfix_set_bc0(lev);

       pc->Redistribute();

    }
}
void
mfix_level::RegridArrays (int lev, BoxArray& new_grids, DistributionMapping& new_dmap)
{
    // ********************************************************************************
    // Cell-based arrays
    // ********************************************************************************

    // Void fraction
    int ng = ep_g[lev]->nGrow();
    std::unique_ptr<MultiFab> ep_g_new(new MultiFab(new_grids,new_dmap,1,ng));
    ep_g_new->copy(*ep_g[lev],0,0,1,ng,ng);
    ep_g_new->FillBoundary(geom[lev].periodicity());
    ep_g[lev] = std::move(ep_g_new);

    // Old void fraction
    ng = ep_go[lev]->nGrow();
    std::unique_ptr<MultiFab> ep_go_new(new MultiFab(new_grids,new_dmap,1,ep_go[lev]->nGrow()));
    ep_go_new->copy(*ep_go[lev],0,0,1,ng,ng);
    ep_go_new->FillBoundary(geom[lev].periodicity());
    ep_go[lev] = std::move(ep_go_new);

    // Gas density
    ng = ro_g[lev]->nGrow();
    std::unique_ptr<MultiFab> ro_g_new(new MultiFab(new_grids,new_dmap,1,ro_g[lev]->nGrow()));
    ro_g_new->copy(*ro_g[lev],0,0,1,ng,ng);
    ro_g_new->FillBoundary(geom[lev].periodicity());
    ro_g[lev] = std::move(ro_g_new);

    // Old gas density
    ng = ro_go[lev]->nGrow();
    std::unique_ptr<MultiFab> ro_go_new(new MultiFab(new_grids,new_dmap,1,ro_go[lev]->nGrow()));
    ro_go_new->copy(*ro_go[lev],0,0,1,ng,ng);
    ro_go_new->FillBoundary(geom[lev].periodicity());
    ro_go[lev] = std::move(ro_go_new);

    // Gas bulk density
    ng = rop_g[lev]->nGrow();
    std::unique_ptr<MultiFab> rop_g_new(new MultiFab(new_grids,new_dmap,1,rop_g[lev]->nGrow()));
    rop_g_new->copy(*rop_g[lev],0,0,1,ng,ng);
    rop_g_new->FillBoundary(geom[lev].periodicity());
    rop_g[lev] = std::move(rop_g_new);

    // Old gas bulk density
    ng = rop_go[lev]->nGrow();
    std::unique_ptr<MultiFab> rop_go_new(new MultiFab(new_grids,new_dmap,1,rop_go[lev]->nGrow()));
    rop_go_new->copy(*rop_go[lev],0,0,1,ng,ng);
    rop_go_new->FillBoundary(geom[lev].periodicity());
    rop_go[lev] = std::move(rop_go_new);

    if (nodal_pressure)
    {
       const BoxArray & nd_grids = amrex::convert(new_grids, IntVect{1,1,1});

       ng = p_g[lev]->nGrow();
       std::unique_ptr<MultiFab> p_g_new(new MultiFab(nd_grids,new_dmap,1,p_g[lev]->nGrow()));
       p_g_new->copy(*p_g[lev],0,0,1,ng,ng);
       p_g_new->FillBoundary(geom[lev].periodicity());
       p_g[lev] = std::move(p_g_new);

       ng = p_go[lev]->nGrow();
       std::unique_ptr<MultiFab> p_go_new(new MultiFab(nd_grids,new_dmap,1,p_go[lev]->nGrow()));
       p_go_new->copy(*p_go[lev],0,0,1,ng,ng);
       p_go_new->FillBoundary(geom[lev].periodicity());
       p_go[lev] = std::move(p_go_new);

       ng = p0_g[lev]->nGrow();
       std::unique_ptr<MultiFab> p0_g_new(new MultiFab(nd_grids,new_dmap,1,ng));
       p0_g_new->copy(*p0_g[lev],0,0,1,ng,ng);
       p0_g_new->FillBoundary(p0_periodicity);
       p0_g[lev] = std::move(p0_g_new);

       ng = pp_g[lev]->nGrow();
       std::unique_ptr<MultiFab> pp_g_new(new MultiFab(nd_grids,new_dmap,1,pp_g[lev]->nGrow()));
       pp_g_new->copy(*pp_g[lev],0,0,1,ng,ng);
       pp_g_new->FillBoundary(geom[lev].periodicity());
       pp_g[lev] = std::move(pp_g_new);

       ng = diveu[lev]->nGrow();
       std::unique_ptr<MultiFab> diveu_new(new MultiFab(nd_grids,new_dmap,1,diveu[lev]->nGrow()));
       diveu[lev] = std::move(diveu_new);
       diveu[lev]->setVal(0.);

    } else {

       ng = p_g[lev]->nGrow();
       std::unique_ptr<MultiFab> p_g_new(new MultiFab(new_grids,new_dmap,1,p_g[lev]->nGrow()));
       p_g_new->copy(*p_g[lev],0,0,1,ng,ng);
       p_g_new->FillBoundary(geom[lev].periodicity());
       p_g[lev] = std::move(p_g_new);

       ng = p_go[lev]->nGrow();
       std::unique_ptr<MultiFab> p_go_new(new MultiFab(new_grids,new_dmap,1,p_go[lev]->nGrow()));
       p_go_new->copy(*p_go[lev],0,0,1,ng,ng);
       p_go_new->FillBoundary(geom[lev].periodicity());
       p_go[lev] = std::move(p_go_new);

       ng = p0_g[lev]->nGrow();
       std::unique_ptr<MultiFab> p0_g_new(new MultiFab(new_grids,new_dmap,1,p0_g[lev]->nGrow()));
       p0_g_new->copy(*p0_g[lev],0,0,1,ng,ng);
       p0_g_new->FillBoundary(p0_periodicity);
       p0_g[lev] = std::move(p0_g_new);

       ng = pp_g[lev]->nGrow();
       std::unique_ptr<MultiFab> pp_g_new(new MultiFab(new_grids,new_dmap,1,pp_g[lev]->nGrow()));
       pp_g_new->copy(*pp_g[lev],0,0,1,ng,ng);
       pp_g_new->FillBoundary(geom[lev].periodicity());
       pp_g[lev] = std::move(pp_g_new);

       ng = diveu[lev]->nGrow();
       std::unique_ptr<MultiFab> diveu_new(new MultiFab(new_grids,new_dmap,1,diveu[lev]->nGrow()));
       diveu[lev] = std::move(diveu_new);
       diveu[lev]->setVal(0.);

    }

    // Molecular viscosity
    ng = mu_g[lev]->nGrow();
    std::unique_ptr<MultiFab> mu_g_new(new MultiFab(new_grids,new_dmap,1,mu_g[lev]->nGrow()));
    mu_g_new->copy(*mu_g[lev],0,0,1,ng,ng);
    mu_g_new->FillBoundary(geom[lev].periodicity());
    mu_g[lev] = std::move(mu_g_new);

    // Lambda
    ng = lambda_g[lev]->nGrow();
    std::unique_ptr<MultiFab> lambda_g_new(new MultiFab(new_grids,new_dmap,1,lambda_g[lev]->nGrow()));
    lambda_g_new->copy(*lambda_g[lev],0,0,1,ng,ng);
    lambda_g_new->FillBoundary(geom[lev].periodicity());
    lambda_g[lev] = std::move(lambda_g_new);

    // Gas velocity
    ng = vel_g[lev]->nGrow();
    std::unique_ptr<MultiFab> vel_g_new(new MultiFab(new_grids,new_dmap,vel_g[lev]->nComp(),ng));
    vel_g_new->copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),ng,ng);
    vel_g_new->FillBoundary(geom[lev].periodicity());
    vel_g[lev] = std::move(vel_g_new);

    // Old gas velocity
    ng = vel_go[lev]->nGrow();
    std::unique_ptr<MultiFab> vel_go_new(new MultiFab(new_grids,new_dmap,vel_go[lev]->nComp(),ng));
    vel_go_new->copy(*vel_go[lev],0,0,vel_go[lev]->nComp(),ng,ng);
    vel_go_new->FillBoundary(geom[lev].periodicity());
    vel_go[lev] = std::move(vel_go_new);

    // Trace(D)
    ng = trD_g[lev]->nGrow();
    std::unique_ptr<MultiFab> trD_g_new(new MultiFab(new_grids,new_dmap,1,trD_g[lev]->nGrow()));
    trD_g[lev] = std::move(trD_g_new);
    trD_g[lev]->setVal(0.);

    // Vorticity
    ng = vort[lev]->nGrow();
    std::unique_ptr<MultiFab> vort_new(new MultiFab(new_grids,new_dmap,1,vort[lev]->nGrow()));
    vort[lev] = std::move(vort_new);
    vort[lev]->setVal(0.);

    // Coefficient in drag
    ng = f_gds[lev]->nGrow();
    std::unique_ptr<MultiFab> f_gds_new(new MultiFab(new_grids,new_dmap,f_gds[lev]->nComp(),ng));
    f_gds[lev] = std::move(f_gds_new);
    f_gds[lev]->setVal(0.);

    // Particle/fluid drag
    ng = drag[lev]->nGrow();
    std::unique_ptr<MultiFab> drag_new(new MultiFab(new_grids,new_dmap,drag[lev]->nComp(),ng));
    drag[lev] = std::move(drag_new);
    drag[lev]->setVal(0.);

   /****************************************************************************
    * Nodal Arrays                                                             *
    ****************************************************************************/

    // Create a nodal BoxArray
    const BoxArray & new_nodal_grids = amrex::convert(new_grids, IntVect{1,1,1});

    // Level-set
    ng = ls[lev]->nGrow();
    std::unique_ptr<MultiFab> ls_new(new MultiFab(new_nodal_grids, new_dmap, 1, ls[lev]->nGrow()));
    ls_new->copy(*ls[lev],0,0,1,ng,ng);
    ls_new->FillBoundary(geom[lev].periodicity());
    ls[lev] = std::move(ls_new);


    // ********************************************************************************
    // Make sure we fill the ghost cells as appropriate -- this is copied from init_fluid
    // ********************************************************************************

    fill_mf_bc(lev,*ep_g[lev]);
    fill_mf_bc(lev,*ep_go[lev]);
    fill_mf_bc(lev,*p_g[lev]);
    fill_mf_bc(lev,*p_go[lev]);
    fill_mf_bc(lev,*ro_g[lev]);
    fill_mf_bc(lev,*ro_go[lev]);
    fill_mf_bc(lev,*rop_g[lev]);
    fill_mf_bc(lev,*rop_go[lev]);

    fill_mf_bc(lev,*mu_g[lev]);
    fill_mf_bc(lev,*lambda_g[lev]);
}
