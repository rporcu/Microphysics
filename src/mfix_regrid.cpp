#include <mfix.H>
#include <MFIX_FLUID_Parms.H>
#include <MFIX_DEM_Parms.H>

void
mfix::Regrid ()
{
  if (ooo_debug) amrex::Print() << "Regrid" << std::endl;
  BL_PROFILE_REGION_START("mfix::Regrid()");

  int base_lev = 0;

  if (load_balance_type == "KnapSack" || load_balance_type == "SFC") // Knapsack and SFC
  {
    amrex::Print() << "Load balancing using " << load_balance_type << std::endl;

    if (DEM::solve)
       AMREX_ALWAYS_ASSERT(m_leveldata[0]->particle_cost != nullptr);
    if (FLUID::solve)
       AMREX_ALWAYS_ASSERT(m_leveldata[0]->fluid_cost != nullptr);

    if (ParallelDescriptor::NProcs() == 1) return;

    if (dual_grid)  //  Beginning of dual grid regridding
    {
      AMREX_ALWAYS_ASSERT(FLUID::solve);

      if (load_balance_fluid > 0)
      {
        for (int lev = base_lev; lev <= finestLevel(); ++lev)
        {
          DistributionMapping new_fluid_dm;

          if ( load_balance_type == "KnapSack" )
          {
            new_fluid_dm = DistributionMapping::makeKnapSack(*m_leveldata[lev]->fluid_cost,
                                                             knapsack_nmax);
          }
          else
          {
            new_fluid_dm = DistributionMapping::makeSFC(*m_leveldata[lev]->fluid_cost, false);
          }

          SetDistributionMap(lev, new_fluid_dm);

          RegridArrays(lev);

          if (m_leveldata[lev]->fluid_cost != nullptr)
            delete m_leveldata[lev]->fluid_cost;

          m_leveldata[lev]->fluid_cost = new MultiFab(grids[lev], new_fluid_dm, 1, 0);
          m_leveldata[lev]->fluid_cost->setVal(0.0);
        }
      }

      mfix_set_p0();
      mfix_set_bc0();

      for (int lev = base_lev; lev <= finestLevel(); ++lev)
      {
        DistributionMapping new_particle_dm;

        if ( load_balance_type == "KnapSack" )
        {
          new_particle_dm = DistributionMapping::makeKnapSack(*m_leveldata[lev]->particle_cost,
                                                              knapsack_nmax);
        }
        else
        {
          new_particle_dm = DistributionMapping::makeSFC(*m_leveldata[lev]->particle_cost,
                                                         false);
        }

        pc->Regrid(new_particle_dm, pc->ParticleBoxArray(lev), lev);

        if (m_leveldata[lev]->particle_cost != nullptr)
          delete m_leveldata[lev]->particle_cost;

        m_leveldata[lev]->particle_cost = new MultiFab(pc->ParticleBoxArray(lev),
                                                       new_particle_dm, 1, 0);
        m_leveldata[lev]->particle_cost->setVal(0.0);

        // This calls re-creates a proper particle_ebfactories
        //  and regrids all the multifabs that depend on it
        if (DEM::solve)
          RegridLevelSetArray(lev);
      }

    }
    else  // Single-grid regridding
    {
      MultiFab costs(grids[base_lev], dmap[base_lev], 1, 0);
      costs.setVal(0.0);

      Print() << "grids = " << grids[base_lev] << std::endl;
      Print() << "costs ba = " << costs.boxArray() << std::endl;

      if(DEM::solve)
        Print() << "particle_cost ba = " 
                << m_leveldata[base_lev]->particle_cost->boxArray()
                << std::endl;

      //Print() << "fluid cost ba = " << fluid_cost[base_lev]->boxArray() << std::endl;

      if (DEM::solve) {
        // costs.plus(* particle_cost[base_lev], 0, 1, 0);

        // MultiFab particle_cost_loc(grids[base_lev], dmap[base_lev], 1, 0);
        // particle_cost_loc.copy(* particle_cost[base_lev], 0, 0, 1);
        MultiFab particle_cost_loc = MFUtil::regrid(grids[base_lev], dmap[base_lev],
                                                    *m_leveldata[base_lev]->particle_cost,
                                                    true);

        costs.plus(particle_cost_loc, 0, 1, 0);
      }
      if (FLUID::solve) {
        // costs.plus(* fluid_cost[base_lev], 0, 1, 0);

        // MultiFab fluid_cost_loc(grids[base_lev], dmap[base_lev], 1, 0);
        // fluid_cost_loc.copy(* fluid_cost[base_lev], 0, 0, 1);
        MultiFab fluid_cost_loc = MFUtil::regrid(grids[base_lev], dmap[base_lev],
                                                 *m_leveldata[base_lev]->fluid_cost,
                                                 true);

        costs.plus(fluid_cost_loc, 0, 1, 0);
      }

      DistributionMapping newdm = DistributionMapping::makeKnapSack(costs,knapsack_nmax);

      SetDistributionMap(base_lev, newdm);

      if (FLUID::solve)
        RegridArrays(base_lev);

      if (FLUID::solve)
      {
        if (m_leveldata[base_lev]->fluid_cost != nullptr)
          delete m_leveldata[base_lev]->fluid_cost;

        m_leveldata[base_lev]->fluid_cost = new MultiFab(grids[base_lev], newdm, 1, 0);
        m_leveldata[base_lev]->fluid_cost->setVal(0.0);
      }

      if (DEM::solve)
      {
        if (m_leveldata[base_lev]->particle_cost != nullptr)
          delete m_leveldata[base_lev]->particle_cost;

        m_leveldata[base_lev]->particle_cost = new MultiFab(grids[base_lev], newdm, 1, 0);
        m_leveldata[base_lev]->particle_cost->setVal(0.0);
      }

      if (DEM::solve){
        pc->Regrid(dmap[base_lev], grids[base_lev], base_lev);
      }

      if (FLUID::solve) mfix_set_bc0();

      // This calls re-creates a proper particles_ebfactory and regrids
      // all the multifab that depend on it
      if (DEM::solve)
        RegridLevelSetArray(base_lev);
      }
  } else {
      amrex::Abort("load_balance_type must be KnapSack or SFC");
  }

  if (DEM::solve)
    for (int i_lev = base_lev; i_lev < nlev; i_lev++)
    {
      // This calls re-creates a proper particle_ebfactories and regrids
      //  all the multifab that depend on it
      RegridLevelSetArray(i_lev);
    }

  // This call resets both the nodal and the diffusion solvers
  if (FLUID::solve)
    mfix_setup_solvers();

  BL_PROFILE_REGION_STOP("mfix::Regrid()");
}
