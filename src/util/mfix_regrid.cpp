#include <mfix.H>
#include <mfix_fluid_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>

namespace {
  void print_process_boxes(amrex::DistributionMapping& dmap)
  {
    const Vector<int>& pmap = dmap.ProcessorMap();
    Vector<Vector<int>> pbox(ParallelDescriptor::NProcs());

    for (unsigned int i=0; i<pmap.size(); ++i)
      pbox[pmap[i]].push_back(i);

    for (unsigned int i=0; i<pbox.size(); ++i) {
      Print() << "Process  " << i << ":";
      for (unsigned int j=0; j<pbox[i].size(); ++j)
        Print() << " " << pbox[i][j];
      Print() << "\n";
    }
  }
}

void
mfix::Regrid ()
{
  if (ooo_debug) amrex::Print() << "Regrid" << std::endl;
  BL_PROFILE_REGION_START("mfix::Regrid()");

  int base_lev = 0;

  if (load_balance_type == "KnapSack" || load_balance_type == "SFC") // Knapsack and SFC
  {
    amrex::Print() << "Load balancing using " << load_balance_type << std::endl;

    if (DEM::solve  || PIC::solve)
       AMREX_ALWAYS_ASSERT(particle_cost[0] != nullptr);

    if (fluid.solve)
       AMREX_ALWAYS_ASSERT(fluid_cost[0] != nullptr);

    if (ParallelDescriptor::NProcs() == 1) return;

    if (dual_grid)  //  Beginning of dual grid regridding
    {
      AMREX_ALWAYS_ASSERT(fluid.solve);

      if (load_balance_fluid > 0)
      {
        for (int lev = base_lev; lev <= finestLevel(); ++lev)
        {
          DistributionMapping new_fluid_dm;

          if ( load_balance_type == "KnapSack" )
          {
            new_fluid_dm = DistributionMapping::makeKnapSack(*fluid_cost[lev],
                                                             knapsack_nmax);
          }
          else
          {
            new_fluid_dm = DistributionMapping::makeSFC(*fluid_cost[lev], false);
          }

          SetDistributionMap(lev, new_fluid_dm);
          print_process_boxes(new_fluid_dm);

          macproj = std::make_unique<MacProjector>(Geom(0,finest_level),
                                         MLMG::Location::FaceCentroid,  // Location of mac_vec
                                         MLMG::Location::FaceCentroid,  // Location of beta
                                         MLMG::Location::CellCenter,    // Location of solution variable phi
                                         MLMG::Location::CellCentroid);// Location of MAC RHS

          RegridArrays(lev);

          if (fluid_cost[lev] != nullptr)
            delete fluid_cost[lev];

          fluid_cost[lev] = new MultiFab(grids[lev], new_fluid_dm, 1, 0);
          fluid_cost[lev]->setVal(0.0);

          if (fluid_proc[lev] != nullptr)
            delete fluid_proc[lev];
          
          const Real proc = static_cast<Real>(ParallelDescriptor::MyProc());
          fluid_proc[lev] = new MultiFab(grids[lev], new_fluid_dm, 1, 0);
          fluid_proc[lev]->setVal(proc);
        }
      }

      mfix_set_p0();
      mfix_set_bc0();

      for (int lev = base_lev; lev <= finestLevel(); ++lev)
      {
        // update new particle grid size
        IntVect max_grid_size = pc->MaxGridSize();
        IntVect new_grid_size = max_grid_size;
        pc->checkParticleBoxSize(lev, new_grid_size, 1.0);

        // if new size is smaller, then downsize particle grids
        // reset the particle cost by # particles
        if (new_grid_size < max_grid_size) {
          pc->downsizeParticleBoxes(lev, new_grid_size);
          pc->resetCostByCount(lev, particle_cost);
          pc->setMaxGridSize(new_grid_size);
        }

        DistributionMapping new_particle_dm;
        if ( load_balance_type == "KnapSack" )
        {
          new_particle_dm = DistributionMapping::makeKnapSack(*particle_cost[lev],
                                                              knapsack_nmax);
        }
        else
        {
          new_particle_dm = DistributionMapping::makeSFC(*particle_cost[lev],
                                                         false);
        }

        pc->Regrid(new_particle_dm, pc->ParticleBoxArray(lev), lev);
        pc->SortParticlesByCell();
        print_process_boxes(new_particle_dm);

        if (particle_cost[lev] != nullptr)
          delete particle_cost[lev];

        particle_cost[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                          new_particle_dm, 1, 0);
        particle_cost[lev]->setVal(0.0);

        // reset rank of particle grids
        if (particle_proc[lev] != nullptr)
          delete particle_proc[lev];
        particle_proc[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                          new_particle_dm, 1, 0);
        const Real proc = static_cast<Real>(ParallelDescriptor::MyProc());
        particle_proc[lev]->setVal(proc);

        // This calls re-creates a proper particle_ebfactories
        //  and regrids all the multifabs that depend on it
        if (DEM::solve || PIC::solve)
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
                << particle_cost[base_lev]->boxArray()
                << std::endl;

      //Print() << "fluid cost ba = " << fluid_cost[base_lev]->boxArray() << std::endl;

      if (DEM::solve) {
        // costs.plus(* particle_cost[base_lev], 0, 1, 0);

        // MultiFab particle_cost_loc(grids[base_lev], dmap[base_lev], 1, 0);
        // particle_cost_loc.copy(* particle_cost[base_lev], 0, 0, 1);
        MultiFab particle_cost_loc = MFUtil::regrid(grids[base_lev], dmap[base_lev],
                                                    *particle_cost[base_lev],
                                                    true);

        costs.plus(particle_cost_loc, 0, 1, 0);
      }
      if (fluid.solve) {
        // costs.plus(* fluid_cost[base_lev], 0, 1, 0);

        // MultiFab fluid_cost_loc(grids[base_lev], dmap[base_lev], 1, 0);
        // fluid_cost_loc.copy(* fluid_cost[base_lev], 0, 0, 1);
        MultiFab fluid_cost_loc = MFUtil::regrid(grids[base_lev], dmap[base_lev],
                                                 *fluid_cost[base_lev],
                                                 true);

        costs.plus(fluid_cost_loc, 0, 1, 0);
      }

      DistributionMapping newdm = DistributionMapping::makeKnapSack(costs,knapsack_nmax);

      SetDistributionMap(base_lev, newdm);

      macproj = std::make_unique<MacProjector>(Geom(0,finest_level),
                                     MLMG::Location::FaceCentroid,  // Location of mac_vec
                                     MLMG::Location::FaceCentroid,  // Location of beta
                                     MLMG::Location::CellCenter,    // Location of solution variable phi
                                     MLMG::Location::CellCentroid);// Location of MAC RHS

      if (fluid.solve)
        RegridArrays(base_lev);

      if (fluid.solve)
      {
        if (fluid_cost[base_lev] != nullptr)
          delete fluid_cost[base_lev];

        fluid_cost[base_lev] = new MultiFab(grids[base_lev], newdm, 1, 0);
        fluid_cost[base_lev]->setVal(0.0);

        if (fluid_proc[base_lev] != nullptr)
          delete fluid_proc[base_lev];

        const Real proc = static_cast<Real>(ParallelDescriptor::MyProc());
        fluid_proc[base_lev] = new MultiFab(grids[base_lev], newdm, 1, 0);
        fluid_proc[base_lev]->setVal(proc);
      }

      if (DEM::solve)
      {
        if (particle_cost[base_lev] != nullptr)
          delete particle_cost[base_lev];

        particle_cost[base_lev] = new MultiFab(grids[base_lev], newdm, 1, 0);
        particle_cost[base_lev]->setVal(0.0);

        if (particle_proc[base_lev] != nullptr)
          delete particle_proc[base_lev];

        const Real proc = static_cast<Real>(ParallelDescriptor::MyProc());
        particle_proc[base_lev] = new MultiFab(grids[base_lev], newdm, 1, 0);
        particle_proc[base_lev]->setVal(proc);
      }

      if (DEM::solve || PIC::solve){
        pc->Regrid(dmap[base_lev], grids[base_lev], base_lev);
      }

      if (fluid.solve) mfix_set_bc0();

      // This calls re-creates a proper particles_ebfactory and regrids
      // all the multifab that depend on it
      if (DEM::solve || PIC::solve)
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
  if (fluid.solve)
    mfix_setup_solvers();

  BL_PROFILE_REGION_STOP("mfix::Regrid()");
}
