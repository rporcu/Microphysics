#include <AMReX_ParmParse.H>

#include <mfix.H>
#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <bc_mod_F.H>
#include <constant_mod_F.H>
#include <mfix_init_fluid.hpp>

#include <AMReX_EBAmrUtil.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_EBFabFactory.H>

#include <MFIX_DEM_Parms.H>

void
mfix::InitParams(int solve_fluid_in, int solve_dem_in, int call_udf_in)
{
    if (ooo_debug) amrex::Print() << "InitParams" << std::endl;
    // set n_error_buf (used in AmrMesh) to default (can overwrite later)
    for (int i = 0; i < n_error_buf.size(); i++)
        n_error_buf[i] = {8,8,8};

    {
        ParmParse pp("mfix");

        // Options to control time stepping
        pp.query("cfl", cfl );

        fixed_dt = -1.;
        pp.query("fixed_dt", fixed_dt );
        pp.query("dt_min", dt_min );
        pp.query("dt_max", dt_max );

        // Options to control verbosity level
        pp.query("verbose", m_verbose);

        // Options to control MLMG behavior
        pp.query( "mg_verbose"             , nodal_mg_verbose );
        pp.query( "mg_cg_verbose"          , nodal_mg_cg_verbose );
        pp.query( "mg_maxiter"             , nodal_mg_maxiter );
        pp.query( "mg_cg_maxiter"          , nodal_mg_cg_maxiter );
        pp.query( "mg_rtol"                , nodal_mg_rtol );
        pp.query( "mg_atol"                , nodal_mg_atol );
        pp.query( "mg_max_coarsening_level", nodal_mg_max_coarsening_level );

        // Default bottom solver here is bicgcg, but alternatives are
        // "smoother", "hypre", "cg", "cgbicg" or "bicgstab"
        nodal_bottom_solver_type = "bicgcg";
        pp.query( "bottom_solver_type",  nodal_bottom_solver_type );

        // Is this a steady-state calculation
        steady_state = 0;
        pp.query( "steady_state", steady_state );

        // Tolerance to check for steady state
        steady_state_tol = -1.;
        pp.query( "steady_state_tol", steady_state_tol );

        // Maximum number of iterations allowed to reach steady state
        pp.query( "steady_state_maxiter", steady_state_maxiter );

        if (steady_state > 0)
        {
           if (steady_state_tol < 0)
              amrex::Abort("Must set steady_state_tol if running to steady state!");

              amrex::Print() << "Running to steady state with maxiter = "
                             << steady_state_maxiter << " and tolerance "
                             << steady_state_tol << std::endl;

        } else {

           if (steady_state_tol > 0)
              amrex::Abort("steady_state_tol set but not steady_state!");
        }

        pp.query("ooo_debug", ooo_debug);

        // The default type is "AsciiFile" but we can over-write that in the inputs file
        //  with "Random"
        pp.query("particle_init_type", particle_init_type);

        // Options to control initial projections (mostly we use these for debugging)
        pp.query("initial_iterations", initial_iterations);
        pp.query("do_initial_proj"   , do_initial_proj);

        pp.query( "advect_density", advect_density );
        pp.query( "advect_tracer" , advect_tracer );
        pp.query( "test_tracer_conservation" , test_tracer_conservation );

        if (test_tracer_conservation && !advect_tracer)
           amrex::Abort("No point in testing tracer conservation with advect_tracer = false");

        if (advect_tracer && !advect_density)
           amrex::Abort("Cant advect tracer without advecting density");

        // The default type is "FixedSize" but we can over-write that in the inputs file
        //  with "KDTree" or "KnapSack"
        pp.query("load_balance_type", load_balance_type);
        pp.query("knapsack_weight_type", knapsack_weight_type);
        pp.query("load_balance_fluid", load_balance_fluid);

        // The default value for the rescaling ratio of the collision time is 50
        pp.query("tcoll_ratio", tcoll_ratio);

        ParmParse pp_mac("mac");
        pp_mac.query( "mg_verbose"   , mac_mg_verbose );
        pp_mac.query( "mg_cg_verbose", mac_mg_cg_verbose );
        pp_mac.query( "mg_rtol"      , mac_mg_rtol );
        pp_mac.query( "mg_atol"      , mac_mg_atol );
        pp_mac.query( "mg_maxiter"   , mac_mg_maxiter );
        pp_mac.query( "mg_cg_maxiter", mac_mg_cg_maxiter );
        pp_mac.query( "mg_max_coarsening_level", mac_mg_max_coarsening_level );

        // Default bottom solver here is bicgcg, but alternatives are
        // "smoother", "hypre", "cg", "cgbicg" or "bicgstab"
        mac_bottom_solver_type = "bicgcg";
        pp_mac.query( "bottom_solver_type", mac_bottom_solver_type );

        AMREX_ALWAYS_ASSERT(load_balance_type == "FixedSize" ||
                            load_balance_type == "KDTree"    ||
                            load_balance_type == "KnapSack"  ||
                            load_balance_type == "SFC"       );

        AMREX_ALWAYS_ASSERT(knapsack_weight_type == "RunTimeCosts" ||
                            knapsack_weight_type == "NumParticles"  );

        // We ensure that we always use the dual grid formulation when using KDTree, and
        //           that we  never use the dual grid formulation otherwise
        if (load_balance_type == "KDTree") {
           dual_grid = true;
        } else {
          ParmParse amr_pp("amr");
          amr_pp.query("dual_grid", dual_grid);
        }

        if (load_balance_type == "KnapSack")
            pp.query("knapsack_nmax", knapsack_nmax);

        // Parameters used be the level-set algorithm. Refer to LSFactory (or
        // mfix.H) for more details:
        //   -> refinement: how well resolved (fine) the (level-set/EB-facet)
        //                  grid needs to be (note: a fine level-set grid means
        //                  that distances and normals are computed accurately)
        //   -> pad:        how many (refined) grid points _outside_ the
        //                  problem domain the grid extends (avoids edge cases
        //                  in physical domain)
        pp.query("levelset__refinement", levelset__refinement);

        // Not needed here... the role of refining EB is filled with AMR level-set
        levelset__eb_refinement = 1;
        // Make sure that a coarsened level-set has a level-set pad of _at least_ 2;
        levelset__pad = 2*levelset__refinement;
        // Ensure that velocity_reconstruction has enough level-set to work off:
        // (2 => EB lives on the same grid resolution as fluid)
        levelset__eb_pad = std::max(2, levelset__pad);

        amrex::Print() << "Auto-generating level-set parameters:" << std::endl
                       << "eb_refinement = " << levelset__eb_refinement << std::endl
                       << "levelset_pad  = " << levelset__pad << std::endl
                       << "eb_pad        = " << levelset__eb_pad << std::endl;

    }

    {
        ParmParse pp_diff("diff");
        pp_diff.query( "mg_verbose"   , diff_mg_verbose );
        pp_diff.query( "mg_cg_verbose", diff_mg_cg_verbose );
        pp_diff.query( "mg_rtol"      , diff_mg_rtol );
        pp_diff.query( "mg_atol"      , diff_mg_atol );
        pp_diff.query( "mg_maxiter"   , diff_mg_maxiter );
        pp_diff.query( "mg_cg_maxiter", diff_mg_cg_maxiter );
        pp_diff.query( "mg_max_coarsening_level", diff_mg_max_coarsening_level );

        // Default bottom solver here is bicgcg, but alternatives are
        // "smoother", "hypre", "cg", "cgbicg" or "bicgstab"
        diff_bottom_solver_type = "bicgcg";
        pp_diff.query( "bottom_solver_type", diff_bottom_solver_type );
    }

    solve_fluid  = solve_fluid_in;
    solve_dem    = solve_dem_in;
    call_udf     = call_udf_in;

    if (solve_dem)
    {
        ParmParse pp("particles");

        pp.query("max_grid_size_x", particle_max_grid_size_x);
        pp.query("max_grid_size_y", particle_max_grid_size_y);
        pp.query("max_grid_size_z", particle_max_grid_size_z);

        // Keep particles that are initially touching the wall. Used by DEM tests.
        pp.query("removeOutOfRange", removeOutOfRange );
    }

    if (solve_dem && !solve_fluid)
    {
        if (fixed_dt <= 0.0)
            amrex::Abort("If running particle-only must specify a positive fixed_dt in the inputs file");
    }

    if (solve_dem && solve_fluid)
    {
      ParmParse pp("mfix");

      std::string drag_type = "None";
      pp.query("drag_type", drag_type);

      if (drag_type == "WenYu")
      {
        m_drag_type = DragType::WenYu;
      }
      else if (drag_type == "Gidaspow")
      {
        m_drag_type = DragType::Gidaspow;
      }
      else if (drag_type == "BVK2")
      {
        m_drag_type = DragType::BVK2;
      }
      else if (drag_type == "UserDrag")
      {
        m_drag_type = DragType::UserDrag;
      }
      else
      {
        amrex::Abort("Don't know this drag_type!");
      }
    }

    {
      ParmParse amr_pp("amr");

      amr_pp.queryarr("avg_p_g",   avg_p_g);
      amr_pp.queryarr("avg_ep_g",  avg_ep_g);
      amr_pp.queryarr("avg_vel_g", avg_vel_g);

      amr_pp.queryarr("avg_vel_p", avg_vel_p);

      // Regions geometry
      amr_pp.queryarr("avg_region_x_e", avg_region_x_e );
      amr_pp.queryarr("avg_region_x_w", avg_region_x_w );
      amr_pp.queryarr("avg_region_y_n", avg_region_y_n );
      amr_pp.queryarr("avg_region_y_s", avg_region_y_s );
      amr_pp.queryarr("avg_region_z_t", avg_region_z_t );
      amr_pp.queryarr("avg_region_z_b", avg_region_z_b );

    }

    get_gravity(gravity);
}


//! Tag using each EB level's volfrac. This requires that the `eb_levels` have
//! already been build.
void mfix::ErrorEst (int lev, TagBoxArray & tags, Real time, int ngrow){
    if (ooo_debug) amrex::Print() << "ErrorEst" << std::endl;
    //___________________________________________________________________________
    // Tag all cells with volfrac \in (0, 1)
    MultiFab volfrac(grids[lev], dmap[lev], 1, 1);
    eb_levels[lev]->fillVolFrac(volfrac, geom[lev]);

    amrex::TagVolfrac(tags, volfrac);
}


void mfix::Init( Real time)
{
    if (ooo_debug) amrex::Print() << "Init" << std::endl;
    InitIOChkData();
    InitIOPltData();

    // Note that finest_level = last level
    finest_level = nlev-1;

    /****************************************************************************
     *                                                                          *
     * Generate levels using ErrorEst tagging. This approach has been based on  *
     * the LSCoreBase::Init().                                                  *
     *                                                                          *
     ***************************************************************************/

    // This tells the AmrMesh class not to iterate when creating the initial
    // grid hierarchy
    SetIterateToFalse();

    // This tells the Cluster routine to use the new chopping routine which
    // rejects cuts if they don't improve the efficiency
    SetUseNewChop();

    /****************************************************************************
     *                                                                          *
     * MFIX-specific grid creation                                              *
     *                                                                          *
     ***************************************************************************/

    // Define coarse level BoxArray and DistributionMap
    const BoxArray& ba = MakeBaseGrids();
    DistributionMapping dm(ba, ParallelDescriptor::NProcs());
    MakeNewLevelFromScratch(0, time, ba, dm);

    for (int lev = 1; lev <= finest_level; lev++)
    {
       if (m_verbose > 0)
            std::cout << "Setting refined region at level " << lev
                      << " to " << grids[lev] << std::endl;

       MakeNewLevelFromScratch(lev, time, grids[lev], dmap[lev]);
    }

    /****************************************************************************
     *                                                                          *
     * Create particle container using mfix::ParGDB                             *
     *                                                                          *
     ***************************************************************************/

    if (solve_dem)
        pc = std::unique_ptr<MFIXParticleContainer>(new MFIXParticleContainer(this));


    /****************************************************************************
     *                                                                          *
     * MFIX-Specific Initialization                                             *
     *                                                                          *
     ***************************************************************************/

    // ******************************************************
    // We only do these at level 0
    // ******************************************************
    Real dx = geom[0].CellSize(0);
    Real dy = geom[0].CellSize(1);
    Real dz = geom[0].CellSize(2);

    Real xlen = geom[0].ProbHi(0) - geom[0].ProbLo(0);
    Real ylen = geom[0].ProbHi(1) - geom[0].ProbLo(1);
    Real zlen = geom[0].ProbHi(2) - geom[0].ProbLo(2);

    Box domain(geom[0].Domain());


    int cyc_x=0, cyc_y=0, cyc_z=0;
    if (geom[0].isPeriodic(0)) cyc_x = 1;
    if (geom[0].isPeriodic(1)) cyc_y = 1;
    if (geom[0].isPeriodic(2)) cyc_z = 1;

    mfix_set_cyclic(&cyc_x, &cyc_y, &cyc_z);

    // Since these involving writing to output files we only do these on the IOProcessor
    if ( ParallelDescriptor::IOProcessor() )
    {
       // Write the initial part of the standard output file
       write_out0(&time, &dx, &dy, &dz, &xlen, &ylen, &zlen,
                  domain.loVect(), domain.hiVect());

       // Write the initial part of the special output file(s)
       write_usr0();
    }

    // Set point sources.
    {
        int err_ps = 0;
        int is_ioproc = 0;
        if ( ParallelDescriptor::IOProcessor() )
            is_ioproc = 1;

        set_ps(&dx,&dy,&dz,&err_ps,&is_ioproc);

        if (err_ps == 1)
            amrex::Abort("Bad data in set_ps");
    }

    for (int lev = 0; lev < nlev; lev++)
        mfix_set_bc_type(lev);
}


BoxArray mfix::MakeBaseGrids () const
{
    if (ooo_debug) amrex::Print() << "MakeBaseGrids" << std::endl;
    BoxArray ba(geom[0].Domain());

    ba.maxSize(max_grid_size[0]);

    // We only call ChopGrids if dividing up the grid using max_grid_size didn't
    //    create enough grids to have at least one grid per processor.
    // This option is controlled by "refine_grid_layout" which defaults to true.

    if ( refine_grid_layout &&
         ba.size() < ParallelDescriptor::NProcs() &&
         (load_balance_type == "FixedSize" || load_balance_type == "KnapSack" || load_balance_type == "SFC") ) {
        ChopGrids(geom[0].Domain(), ba, ParallelDescriptor::NProcs());
    }

    if (ba == grids[0]) {
        ba = grids[0];  // to avoid dupliates
    }
    amrex::Print() << "In MakeBaseGrids: BA HAS " << ba.size() << " GRIDS " << std::endl;
    return ba;
}


void mfix::ChopGrids (const Box& domain, BoxArray& ba, int target_size) const
{
    if (ooo_debug) amrex::Print() << "ChopGrids" << std::endl;
    if ( ParallelDescriptor::IOProcessor() )
       amrex::Warning("Using max_grid_size only did not make enough grids for the number of processors");

    // Here we hard-wire the maximum number of times we divide the boxes.
    int n = 10;

    // Here we hard-wire the minimum size in any one direction the boxes can be
    int min_grid_size = 4;

    IntVect chunk(domain.length(0),domain.length(1),domain.length(2));

    int j;
    for (int cnt = 1; cnt <= n; ++cnt)
    {
        if (chunk[0] >= chunk[1] && chunk[0] >= chunk[2])
        {
            j = 0;
        }
        else if (chunk[1] >= chunk[0] && chunk[1] >= chunk[2])
        {
            j = 1;
        }
        else if (chunk[2] >= chunk[0] && chunk[2] >= chunk[1])
        {
            j = 2;
        }
        chunk[j] /= 2;

        if (chunk[j] >= min_grid_size)
        {
            ba.maxSize(chunk);
        }
        else
        {
            // chunk[j] was the biggest chunk -- if this is too small then we're done
            if ( ParallelDescriptor::IOProcessor() )
               amrex::Warning("ChopGrids was unable to make enough grids for the number of processors");
            return;
        }

        // Test if we now have enough grids
        if (ba.size() >= target_size) return;
    }
}


void mfix::MakeNewLevelFromScratch (int lev, Real time,
                                    const BoxArray& new_grids, const DistributionMapping& new_dmap)
{
    if (ooo_debug) amrex::Print() << "MakeNewLevelFromScratch" << std::endl;
    if (m_verbose > 0)
    {
        std::cout << "MAKING NEW LEVEL " << lev << std::endl;
        std::cout << "WITH BOX ARRAY   " << new_grids << std::endl;
    }

    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    amrex::Print() << "SETTING NEW GRIDS IN MAKE NEW LEVEL " << new_grids << std::endl;
    amrex::Print() << "SETTING NEW DMAP IN MAKE NEW LEVEL " << new_dmap << std::endl;

    if (lev == 0)
    {
        // This is being done by mfix::make_eb_geometry, otherwise it would be
        // here
        MakeBCArrays();
        check_data();

        Real dx = geom[lev].CellSize(0);
        Real dy = geom[lev].CellSize(1);
        Real dz = geom[lev].CellSize(2);

        // This is separate from check_data because it is only called on
        // initialization, not on restart
        Box domain(geom[0].Domain());
        if ( ParallelDescriptor::IOProcessor() )
            check_initial_conditions(&dx,&dy,&dz,domain.loVect(),domain.hiVect());
    }
}


void mfix::ReMakeNewLevelFromScratch (int lev,
                                 const BoxArray & new_grids, const DistributionMapping & new_dmap)
{
    if (ooo_debug) amrex::Print() << "ReMakeNewLevelFromScratch" << std::endl;
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

    if (lev == 0)
    {
       MakeBCArrays();
       check_data();
    }

    // We need to re-fill these arrays for the larger domain (after replication).
    mfix_set_bc_type(lev);
}


void mfix::check_data ()
{
    if (ooo_debug) amrex::Print() << "check_data" << std::endl;
    Real dx = geom[0].CellSize(0);
    Real dy = geom[0].CellSize(1);
    Real dz = geom[0].CellSize(2);

    Box domain(geom[0].Domain());

    // Only call this check on one processor since it has a bunch of print statements
    if ( ParallelDescriptor::IOProcessor() )
    {
//     check_boundary_conditions(&dx,&dy,&dz,&xlen,&ylen,&zlen,domain.loVect(),domain.hiVect());
       check_point_sources(&dx,&dy,&dz);
       check_bc_flow();
    }
}


void mfix::InitLevelData(Real time)
{
    if (ooo_debug) amrex::Print() << "InitLevelData" << std::endl;
    // Allocate the fluid data, NOTE: this depends on the ebfactories.
    if (solve_fluid)
       for (int lev = 0; lev < nlev; lev++)
          AllocateArrays(lev);

    // Allocate the particle data
    if (solve_dem)
    {
      Real strt_init_part = ParallelDescriptor::second();

      pc->AllocData();

      if (particle_init_type == "AsciiFile")
      {
        amrex::Print() << "Reading particles from particle_input.dat ..." << std::endl;
        pc->InitParticlesAscii("particle_input.dat");

      } else if (particle_init_type == "Random")
      {
        int n_per_cell = 1;

        amrex::Print() << "Randomly initializing " << n_per_cell
                       << " particles per cell ..."
                       << std::endl;

        Real  radius = 1.0;
        Real  volume = 1.0;
        Real    mass = 1.0;
        Real density = 1.0;
        Real    omoi = 1.0;
        Real    velx = 0.0;
        Real    vely = 0.0;
        Real    velz = 0.0;
        Real   dragx = 0.0;
        Real   dragy = 0.0;
        Real   dragz = 0.0;
        Real  omegax = 0.0;
        Real  omegay = 0.0;
        Real  omegaz = 0.0;
        int    phase = 1;
        int    state = 0;

        MFIXParticleContainer::ParticleInitData pdata = {radius,volume,mass,density,omoi,
                                                         velx,vely,velz,omegax,omegay,omegaz,
                                                         dragx,dragy,dragz,phase,state};

        pc->InitNRandomPerCell(n_per_cell, pdata);
        pc->WriteAsciiFileForInit ("random_particles");
        exit(0);

      } else if (particle_init_type == "Auto") {

         amrex::Print() << "Auto generating particles ..." << std::endl;

         pc->InitParticlesAuto();

      } else {

         amrex::Abort("Bad particle_init_type");
      }

      pc->Redistribute();

      // used in load balancing
      if (load_balance_type == "KnapSack" || load_balance_type == "SFC")
      {
          for (int lev = 0; lev < nlev; lev++)
          {
             particle_cost[lev].reset(new MultiFab(pc->ParticleBoxArray(lev),
                                                   pc->ParticleDistributionMap(lev), 1, 0));
             particle_cost[lev]->setVal(0.0);
          }
      }

      Real end_init_part = ParallelDescriptor::second() - strt_init_part;
      ParallelDescriptor::ReduceRealMax(end_init_part, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "Time spent in initializing particles " << end_init_part << std::endl;
    }

    if (solve_fluid)
    {
       if (load_balance_type == "KnapSack" || load_balance_type == "SFC")
       {
          for (int lev = 0; lev < nlev; lev++)
          {
             fluid_cost[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 0));
             fluid_cost[lev]->setVal(0.0);
          }
       }
    }
}

void
mfix::PostInit(Real& dt, Real time, int restart_flag, Real stop_time)
{
    if (ooo_debug) amrex::Print() << "PostInit" << std::endl;
    if (solve_dem)
    {
        // Auto generated particles may be out of the domain. This call will
        // remove them. Note that this has to occur after the EB geometry is
        // created. if (particle_init_type == "Auto" && !restart_flag &&
        // particle_ebfactory[finest_level])

      if (removeOutOfRange)
        {

          if ((nlev == 1) && (!restart_flag && particle_ebfactory[finest_level]))
          {
            //___________________________________________________________________
            // Only 1 refined level-set

            Print() << "Clean up auto-generated particles.\n" << std::endl;

            const MultiFab * ls_data = level_sets[1].get();
            iMultiFab ls_valid(ls_data->boxArray(), ls_data->DistributionMap(),
                               ls_data->nComp(), ls_data->nGrow());

            pc->RemoveOutOfRange(finest_level, particle_ebfactory[finest_level].get(),
                                 ls_data, levelset__refinement);
          }
          else if (!restart_flag && particle_ebfactory[finest_level])
          {
            //___________________________________________________________________
            // Multi-level everything

            Print() << "Clean up auto-generated particles.\n" << std::endl;

            for (int ilev = 0; ilev < nlev; ilev ++)
            {
              const MultiFab * ls_data = level_sets[ilev].get();
              iMultiFab ls_valid(ls_data->boxArray(), ls_data->DistributionMap(),
                                 ls_data->nComp(), ls_data->nGrow());

              pc->RemoveOutOfRange(ilev, particle_ebfactory[ilev].get(),
                                   ls_data, 1);
            }
          }

        }


        for (int lev = 0; lev < nlev; lev++)
        {
            // We need to do this *after* restart (hence putting this here not
            // in Init) because we may want to move from KDTree to Knapsack, or
            // change the particle_max_grid_size on restart.
            if ( (load_balance_type == "KnapSack" || load_balance_type == "SFC") &&
                 dual_grid && particle_max_grid_size_x > 0
                           && particle_max_grid_size_y > 0
                           && particle_max_grid_size_z > 0)
            {
                BoxArray particle_ba(geom[lev].Domain());
                IntVect particle_max_grid_size(particle_max_grid_size_x,
                                               particle_max_grid_size_y,
                                               particle_max_grid_size_z);
                particle_ba.maxSize(particle_max_grid_size);
                DistributionMapping particle_dm(particle_ba, ParallelDescriptor::NProcs());
                pc->Regrid(particle_dm, particle_ba);

                particle_cost[lev].reset(new MultiFab(pc->ParticleBoxArray(lev),
                                                      pc->ParticleDistributionMap(lev), 1, 0));
                particle_cost[lev]->setVal(0.0);

                // This calls re-creates a proper particle_ebfactories
                //  and regrids all the multifabs that depend on it
                if (solve_dem)
                    RegridLevelSetArray(lev);

            }
        }

        Real min_dp[10], min_ro[10];
        Real max_dp[10], max_ro[10];
        Real avg_dp[10], avg_ro[10];

        pc->GetParticleAvgProp( min_dp, min_ro,
                                max_dp, max_ro,
                                avg_dp, avg_ro );

        init_collision(min_dp, min_ro,
                       max_dp, max_ro,
                       avg_dp, avg_ro,
                       tcoll_ratio);

        DEMParams::Initialize();

        if (!solve_fluid)
            dt = fixed_dt;
    }

    if (solve_fluid)
        mfix_init_fluid(restart_flag,dt,stop_time);

    // Call user-defined subroutine to set constants, check data, etc.
    if (call_udf) mfix_usr0();
}

void
mfix::MakeBCArrays ()
{
    if (ooo_debug) amrex::Print() << "MakeBCArrays" << std::endl;
    bc_ilo.resize(nlev);
    bc_ihi.resize(nlev);
    bc_jlo.resize(nlev);
    bc_jhi.resize(nlev);
    bc_klo.resize(nlev);
    bc_khi.resize(nlev);

    for (int lev = 0; lev < nlev; lev++)
    {
       // Define and allocate the integer MultiFab that is the outside adjacent cells of the problem domain.
       Box domainx(geom[lev].Domain());
       domainx.grow(1,nghost);
       domainx.grow(2,nghost);
       Box box_ilo = amrex::adjCellLo(domainx,0,1);
       Box box_ihi = amrex::adjCellHi(domainx,0,1);

       Box domainy(geom[lev].Domain());
       domainy.grow(0,nghost);
       domainy.grow(2,nghost);
       Box box_jlo = amrex::adjCellLo(domainy,1,1);
       Box box_jhi = amrex::adjCellHi(domainy,1,1);

       Box domainz(geom[lev].Domain());
       domainz.grow(0,nghost);
       domainz.grow(1,nghost);
       Box box_klo = amrex::adjCellLo(domainz,2,1);
       Box box_khi = amrex::adjCellHi(domainz,2,1);

       // Note that each of these is a single IArrayBox so every process has a copy of them
       bc_ilo[lev].reset(new IArrayBox(box_ilo,2));
       bc_ihi[lev].reset(new IArrayBox(box_ihi,2));
       bc_jlo[lev].reset(new IArrayBox(box_jlo,2));
       bc_jhi[lev].reset(new IArrayBox(box_jhi,2));
       bc_klo[lev].reset(new IArrayBox(box_klo,2));
       bc_khi[lev].reset(new IArrayBox(box_khi,2));
   }
}

void
mfix::mfix_init_fluid( int is_restarting, Real dt, Real stop_time)
{
    if (ooo_debug) amrex::Print() << "mfix_init_fluid" << std::endl;

    Real xlen = geom[0].ProbHi(0) - geom[0].ProbLo(0);
    Real ylen = geom[0].ProbHi(1) - geom[0].ProbLo(1);
    Real zlen = geom[0].ProbHi(2) - geom[0].ProbLo(2);

    // Here we set bc values for p and u,v,w before the IC's are set
    mfix_set_bc0();

    for (int lev = 0; lev < nlev; lev++)
    {
       Box domain(geom[lev].Domain());

       Real dx = geom[lev].CellSize(0);
       Real dy = geom[lev].CellSize(1);
       Real dz = geom[lev].CellSize(2);

       // We deliberately don't tile this loop since we will be looping
       //    over bc's on faces and it makes more sense to do this one grid at a time
       for (MFIter mfi(*ep_g[lev],false); mfi.isValid(); ++mfi) {

          const Box& bx = mfi.validbox();
          const Box& sbx = (*ep_g[lev])[mfi].box();

        if ( is_restarting ) {

            init_fluid_restart(bx, (*mu_g[lev])[mfi]);

          } else {

            init_fluid(sbx, bx, domain,
                       (*ep_g[lev])[mfi], (*ro_g[lev])[mfi], 
                       (*trac[lev])[mfi], (*p_g[lev])[mfi],
                       (*vel_g[lev])[mfi], (*mu_g[lev])[mfi],
                       dx, dy, dz, xlen, ylen, zlen, test_tracer_conservation);
          }
       }

       // Make sure to fill the "old state" before we start ...
       MultiFab::Copy(  *ro_go[lev], *ro_g[lev], 0, 0, 1, 0);
       MultiFab::Copy( *trac_o[lev], *trac[lev], 0, 0, 1, 0);
    }

    mfix_set_p0();

  // Here we re-set the bc values for p and u,v,w just in case init_fluid
  //      over-wrote some of the bc values with ic values
  mfix_set_bc0();

  for (int lev = 0; lev < nlev; lev++)
  {
     ep_g[lev]->FillBoundary(geom[lev].periodicity());
     ro_g[lev]->FillBoundary(geom[lev].periodicity());
     mu_g[lev]->FillBoundary(geom[lev].periodicity());

    if (advect_tracer)
       trac[lev]->FillBoundary(geom[lev].periodicity());

     vel_g[lev]->FillBoundary(geom[lev].periodicity());
  }

  if (is_restarting == 0)
  {
     // Just for reference, we compute the volume inside the EB walls (as if there were no particles)
     ep_g[0]->setVal(1.0);

     sum_vol_orig = volWgtSum(0,*ep_g[0],0);

     Print() << "Enclosed domain volume is   " << sum_vol_orig << std::endl;

     Real domain_vol = sum_vol_orig;

     // Now initialize the volume fraction ep_g before the first projection
     mfix_calc_volume_fraction(sum_vol_orig);
     Print() << "Setting original sum_vol to " << sum_vol_orig << std::endl;

     Print() << "Difference is   " << (domain_vol - sum_vol_orig) << std::endl;

     // This sets bcs for ep_g and mu_g
     Real time = 0.0;
     mfix_set_scalar_bcs(time,ro_g ,trac  ,ep_g,mu_g);
     mfix_set_scalar_bcs(time,ro_go,trac_o,ep_g,mu_g);

     // Project the initial velocity field
     if (do_initial_proj)
        mfix_project_velocity();

     // Iterate to compute the initial pressure
     if (initial_iterations > 0)

        mfix_initial_iterations(dt,stop_time);

     } else {

        //Calculation of sum_vol_orig for a restarting point  
        sum_vol_orig = volWgtSum(0,*ep_g[0],0);

        Print() << "Setting original sum_vol to " << sum_vol_orig << std::endl;
     }
}

void
mfix::mfix_set_bc0()
{
    if (ooo_debug) amrex::Print() << "mfix_set_bc0" << std::endl;
    for (int lev = 0; lev < nlev; lev++)
    {

     Box domain(geom[lev].Domain());

     // Don't tile this -- at least for now
     for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
     {
       const Box& sbx = (*ep_g[lev])[mfi].box();

       set_bc0(sbx, &mfi, lev, domain);
     }

     ep_g[lev]->FillBoundary(geom[lev].periodicity());
     ro_g[lev]->FillBoundary(geom[lev].periodicity());
     if (advect_tracer)
        trac[lev]->FillBoundary(geom[lev].periodicity());
   }

   // Put velocity Dirichlet bc's on faces
   Real time = 0.0;
   int extrap_dir_bcs = 0;
   mfix_set_velocity_bcs(time,vel_g,extrap_dir_bcs);

  for (int lev = 0; lev < nlev; lev++)
     vel_g[lev]->FillBoundary(geom[lev].periodicity());
}

void
mfix::mfix_set_p0()
{
  if (ooo_debug) amrex::Print() << "mfix_set_p0" << std::endl;
  Real xlen = geom[0].ProbHi(0) - geom[0].ProbLo(0);
  Real ylen = geom[0].ProbHi(1) - geom[0].ProbLo(1);
  Real zlen = geom[0].ProbHi(2) - geom[0].ProbLo(2);

  int delp_dir;
  set_delp_dir(&delp_dir);

  IntVect press_per = IntVect(geom[0].isPeriodic(0),geom[0].isPeriodic(1),geom[0].isPeriodic(2));

  // Here we set a separate periodicity flag for p0_g because when we use
  // pressure drop (delp) boundary conditions we fill all variables *except* p0
  // periodically
  if (delp_dir > -1) press_per[delp_dir] = 0;
  p0_periodicity = Periodicity(press_per);

  // Initialize gp0 to 0
  gp0[0] = 0.0;
  gp0[1] = 0.0;
  gp0[2] = 0.0;

  for (int lev = 0; lev < nlev; lev++)
  {

     Real dx = geom[lev].CellSize(0);
     Real dy = geom[lev].CellSize(1);
     Real dz = geom[lev].CellSize(2);

     Box domain(geom[lev].Domain());

     // We put this outside the MFIter loop because we need gp0 even on ranks with no boxes
     // because we will use it in computing dt separately on every rank
     set_gp0(domain.loVect(), domain.hiVect(),
             gp0,
             &dx, &dy, &dz, &xlen, &ylen, &zlen, 
             bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
             bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
             bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
             &nghost, &delp_dir);

     // We deliberately don't tile this loop since we will be looping
     //    over bc's on faces and it makes more sense to do this one grid at a time
     for (MFIter mfi(*ep_g[lev],false); mfi.isValid(); ++mfi)
     {
       const Box& bx = mfi.validbox();

       set_p0(bx, &mfi, lev, domain, xlen, ylen, zlen, delp_dir);
     }

     p0_g[lev]->FillBoundary(p0_periodicity);
   }
}

void mfix::mfix_set_ls_near_inflow()
{
    if (ooo_debug) amrex::Print() << "mfix_ls_near_inflow" << std::endl;
    // This function is a bit Wonky... TODO: figure out why we need + nghost
    // (it's late at the moment, so I can't figure it it out right now...)
    const int levelset_nghost = levelset__eb_pad + nghost;

    if (nlev > 1)
    {
        // The level set is always at the same level of refinement as the
        // boundary condition arrays
        int n = 1;

        for (int lev = 0; lev < nlev; lev++)
        {
            Box domain(geom[lev].Domain());

            MultiFab* ls_phi = level_sets[lev].get();
            const Real* dx   = geom[lev].CellSize();

            // Don't tile this
            for (MFIter mfi(*ls_phi); mfi.isValid(); ++mfi)
            {
                FArrayBox & ls_fab = (* ls_phi)[mfi];

                set_ls_inflow(lev, ls_fab, domain, &levelset_nghost, n, dx);
            }
        }
    }
    else
    {
        // Here we assume that the particles and fluid only exist on one level
        int lev = 0;
        int lev_ref = 1;

        // ... but that the level set may be at a finer resolution.
        int n = levelset__refinement;
        {
            Box domain(geom[lev].Domain());
            const Real * dx   = geom[lev].CellSize();
            MultiFab * ls_phi = level_sets[lev_ref].get();

            // Don't tile this
            for (MFIter mfi(* ls_phi); mfi.isValid(); ++mfi)
            {
                FArrayBox & ls_fab = (* ls_phi)[mfi];

                set_ls_inflow( lev, ls_fab, domain, &levelset_nghost, n, dx);
            }
        }

        // ... now also "fix" the level-zero level-set (for consistency with
        // multi-level levelset)
        n = 1;
        {
            Box domain(geom[lev].Domain());
            const Real * dx   = geom[lev].CellSize();
            MultiFab * ls_phi = level_sets[lev].get();

            // Don't tile this
            for (MFIter mfi(* ls_phi); mfi.isValid(); ++mfi)
            {
                FArrayBox & ls_fab = (* ls_phi)[mfi];

                set_ls_inflow( lev, ls_fab, domain, &levelset_nghost, n, dx);
            }
        }
    }
}
