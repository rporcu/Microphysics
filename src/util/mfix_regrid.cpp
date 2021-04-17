#include <mfix.H>
#include <mfix_fluid_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>

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

    if (FLUID::solve)
       AMREX_ALWAYS_ASSERT(fluid_cost[0] != nullptr);

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
            new_fluid_dm = DistributionMapping::makeKnapSack(*fluid_cost[lev],
                                                             knapsack_nmax);
          }
          else
          {
            new_fluid_dm = DistributionMapping::makeSFC(*fluid_cost[lev], false);
          }

          SetDistributionMap(lev, new_fluid_dm);

          macproj = std::make_unique<MacProjector>(Geom(0,finest_level),
                                         MLMG::Location::FaceCentroid,  // Location of mac_vec
                                         MLMG::Location::FaceCentroid,  // Location of beta
                                         MLMG::Location::CellCenter,    // Location of solution variable phi
                                         MLMG::Location::CellCentroid);// Location of MAC RHS

          RegridArrays(lev);

          // reset the ranks of fluid grids
          if (m_leveldata[lev]->ba_proc != nullptr)
            delete m_leveldata[lev]->ba_proc;

          m_leveldata[lev]->ba_proc = new MultiFab(grids[lev], new_fluid_dm, 1, 0,
                                                   MFInfo(), *ebfactory[lev]);
          const Real proc = Real(ParallelDescriptor::MyProc());
          for (MFIter mfi(*(m_leveldata[lev]->ba_proc), false); mfi.isValid(); ++mfi)
          {
            Array4<Real> const& bx_proc = m_leveldata[lev]->ba_proc->array(mfi);
            ParallelFor(mfi.validbox(), [bx_proc, proc] 
                        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        { bx_proc(i,j,k) = proc; });
          }

          if (fluid_cost[lev] != nullptr)
            delete fluid_cost[lev];

          fluid_cost[lev] = new MultiFab(grids[lev], new_fluid_dm, 1, 0);
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
          new_particle_dm = DistributionMapping::makeKnapSack(*particle_cost[lev],
                                                              knapsack_nmax);
        }
        else
        {
          new_particle_dm = DistributionMapping::makeSFC(*particle_cost[lev],
                                                         false);
        }

        pc->Regrid(new_particle_dm, pc->ParticleBoxArray(lev), lev);

        if (particle_cost[lev] != nullptr)
          delete particle_cost[lev];

        particle_cost[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                          new_particle_dm, 1, 0);
        particle_cost[lev]->setVal(0.0);

        // reset rank of particle grids
        if (particle_ba_proc[lev] != nullptr)
          delete particle_ba_proc[lev];
        particle_ba_proc[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                             new_particle_dm, 1, 0);
        const Real proc = Real(ParallelDescriptor::MyProc());
        for (MFIter mfi(*(particle_ba_proc[lev]), false); mfi.isValid(); ++mfi)
        {
          amrex::Array4<Real> const& par_bx_proc = particle_ba_proc[lev]->array(mfi);
          ParallelFor(mfi.validbox(), [par_bx_proc, proc]
                      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                      { par_bx_proc(i,j,k) = proc; });
        }

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
      if (FLUID::solve) {
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

      if (FLUID::solve)
        RegridArrays(base_lev);

      if (FLUID::solve)
      {
        if (fluid_cost[base_lev] != nullptr)
          delete fluid_cost[base_lev];

        fluid_cost[base_lev] = new MultiFab(grids[base_lev], newdm, 1, 0);
        fluid_cost[base_lev]->setVal(0.0);
      }

      if (DEM::solve)
      {
        if (particle_cost[base_lev] != nullptr)
          delete particle_cost[base_lev];

        particle_cost[base_lev] = new MultiFab(grids[base_lev], newdm, 1, 0);
        particle_cost[base_lev]->setVal(0.0);
      }

      if (DEM::solve || PIC::solve){
        pc->Regrid(dmap[base_lev], grids[base_lev], base_lev);
      }

      if (FLUID::solve) mfix_set_bc0();

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
  if (FLUID::solve)
    mfix_setup_solvers();

  BL_PROFILE_REGION_STOP("mfix::Regrid()");
}
