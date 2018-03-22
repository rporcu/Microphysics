#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

void
mfix_level::InitParams(int solve_fluid_in, int solve_dem_in,
                       int max_nit_in, int call_udf_in)
{
    {
        ParmParse pp("mfix");

        // Whether to use projection method
        pp.query("use_proj_method", use_proj_method );

        // Options to control time stepping
        pp.query("cfl", cfl );
        pp.query("fixed_dt", fixed_dt );

        // Option to control MGML behavior
        pp.query( "mg_verbose", mg_verbose );
        pp.query( "mg_cg_verbose", mg_cg_verbose );
        pp.query( "mg_max_iter", mg_max_iter );
        pp.query( "mg_cg_maxiter", mg_cg_maxiter );
        pp.query( "mg_max_fmg_iter", mg_max_fmg_iter );
        pp.query( "mg_rtol", mg_rtol );
        pp.query( "mg_atol", mg_atol );

        // Tolerance to check for steady state (projection only)
        pp.query( "steady_state_tol", steady_state_tol );

        // The default type is "AsciiFile" but we can over-write that in the inputs file
        //  with "Random"
        pp.query("particle_init_type", particle_init_type);

        // The default type is "FixedSize" but we can over-write that in the inputs file
        //  with "KDTree" or "KnapSack"
        pp.query("load_balance_type", load_balance_type);
        pp.query("knapsack_weight_type" , knapsack_weight_type);

        pp.query("dual_grid", dual_grid);

        AMREX_ALWAYS_ASSERT(load_balance_type == "FixedSize" ||
                            load_balance_type == "KDTree"    ||
                            load_balance_type == "KnapSack");

        AMREX_ALWAYS_ASSERT(knapsack_weight_type == "RunTimeCosts" ||
                            knapsack_weight_type == "NumParticles"  );

        // If subdt_io is true, des_time_loop calls output_manager
        subdt_io = false; // default to false (if not present in inputs file)
        pp.query("subdt_io", subdt_io);

        // If true, then compute particle/EB collisions directly using neighbouring eb-facets
        // WARNING: this mode can be slow, and EB-facets can sometimes not be "water-tight"
        pp.query("legacy__eb_collisions", legacy__eb_collisions);
        
        // Parameters used be the level-set algorithm. Refer to LSFactory (or mfix_level.H) for more details:
        //   -> refinement: how well resolved (fine) the (level-set/EB-facet) grid needs to be
        //                  (note: a fine level-set grid means that distances and normals are computed accurately)
        //   -> pad:        how many (refined) grid points _outside_ the problem domain the grid extends
        //                  (avoids edge cases in physical domain)
        pp.query("levelset__refinement", levelset__refinement);
        pp.query("levelset__eb_refinement", levelset__eb_refinement);
        pp.query("levelset__pad", levelset__pad);
        pp.query("levelset__eb_pad", levelset__eb_pad);
    }

    solve_fluid  = solve_fluid_in;
    solve_dem    = solve_dem_in;
    max_nit      = max_nit_in;
    call_udf     = call_udf_in;

    {
        ParmParse pp("amr");
        pp.query("dual_grid", dual_grid);
    }

    {
        ParmParse pp("particles");
        pp.query("max_grid_size_x", particle_max_grid_size_x);
        pp.query("max_grid_size_y", particle_max_grid_size_y);
        pp.query("max_grid_size_z", particle_max_grid_size_z);
    }
}

void mfix_level::Init(int lev, Real dt, Real time)
{
    BL_ASSERT(max_level == 0);

    InitIOData();

    // Define coarse level BoxArray and DistributionMap
    finest_level = 0;

    const BoxArray& ba = MakeBaseGrids();
    DistributionMapping dm(ba, ParallelDescriptor::NProcs());

    MakeNewLevelFromScratch(0, time, ba, dm);

    if (dual_grid                    &&
        particle_max_grid_size_x > 0 &&
        particle_max_grid_size_y > 0 &&
        particle_max_grid_size_z > 0) {
        BoxArray particle_ba(geom[lev].Domain());
        IntVect particle_max_grid_size(particle_max_grid_size_x,
                                       particle_max_grid_size_y,
                                       particle_max_grid_size_z);
        particle_ba.maxSize(particle_max_grid_size);
        DistributionMapping particle_dm(particle_ba, ParallelDescriptor::NProcs());
        pc->Regrid(particle_dm, particle_ba);
    }

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
    Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
    Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

    Box domain(geom[0].Domain());

    int cyc_x=0, cyc_y=0, cyc_z=0;
    if(geom[lev].isPeriodic(0)) cyc_x = 1;
    if(geom[lev].isPeriodic(1)) cyc_y = 1;
    if(geom[lev].isPeriodic(2)) cyc_z = 1;

    mfix_set_cyclic(&cyc_x, &cyc_y, &cyc_z);

    // Since these involving writing to output files we only do these on the IOProcessor
    if ( ParallelDescriptor::IOProcessor() )
    {
       // Write the initial part of the standard output file
       write_out0(&time, &dt, &dx, &dy, &dz, &xlen, &ylen, &zlen,
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

    mfix_set_bc_type(lev);

    // Allocate container for eb-normals
    eb_normals = unique_ptr<MultiFab>(new MultiFab);
    dummy = unique_ptr<MultiFab>(new MultiFab);

    // Level-Set: initialize container for level set
    // level-set MultiFab is defined here, and set to (fortran) huge(amrex_real)
    //            -> use min to intersect new eb boundaries (in update)
    level_set = std::unique_ptr<LSFactory>(new LSFactory(lev, levelset__refinement, levelset__eb_refinement,
                                                              levelset__pad, levelset__eb_pad, pc.get()));

}

BoxArray
mfix_level::MakeBaseGrids () const
{
    BoxArray ba(geom[0].Domain());

    ba.maxSize(max_grid_size[0]);

    // We only call ChopGrids if dividing up the grid using max_grid_size didn't
    //    create enough grids to have at least one grid per processor.
    // This option is controlled by "refine_grid_layout" which defaults to true.

    if ( refine_grid_layout &&
         ba.size() < ParallelDescriptor::NProcs() &&
         (load_balance_type == "FixedSize" || load_balance_type == "KnapSack") ) {
        ChopGrids(geom[0].Domain(), ba, ParallelDescriptor::NProcs());
    }

    if (ba == grids[0]) {
        ba = grids[0];  // to avoid dupliates
    }
    amrex::Print() << "In MakeBaseGrids: BA HAS " << ba.size() << " GRIDS " << std::endl;
    return ba;
}

void
mfix_level::ChopGrids (const Box& domain, BoxArray& ba, int target_size) const
{
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

void
mfix_level::MakeNewLevelFromScratch (int lev, Real time,
             const BoxArray& new_grids, const DistributionMapping& new_dmap)

{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

    MakeBCArrays();

    check_data(lev);

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    // This is separate from check_data because it is only called on initialization,
    // not on restart
    Box domain(geom[0].Domain());
    if ( ParallelDescriptor::IOProcessor() )
       check_initial_conditions(&dx,&dy,&dz,domain.loVect(),domain.hiVect());
}

void
mfix_level::ReMakeNewLevelFromScratch (int lev, const BoxArray& new_grids, const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

    MakeBCArrays();

    check_data(lev);

    // We need to re-fill these arrays for the larger domain (after replication).
    mfix_set_bc_type(lev);
}

void
mfix_level::check_data (int lev)
{

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
    Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
    Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

    Box domain(geom[0].Domain());

    // Convert (mass, volume) flows to velocities.
    set_bc_flow(&xlen, &ylen, &zlen, &dx,&dy,&dz);

    // Only call this check on one processor since it has a bunch of print statements
    if ( ParallelDescriptor::IOProcessor() )
    {
       check_boundary_conditions(&dx,&dy,&dz,&xlen,&ylen,&zlen,domain.loVect(),domain.hiVect());
       check_point_sources(&dx,&dy,&dz);
       check_bc_flow();
    }
}

void
mfix_level::AllocateArrays (int lev)
{
    int nghost = 2;

    // ********************************************************************************
    // Cell-based arrays
    // ********************************************************************************

    // Void fraction
    ep_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    ep_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    ep_g[lev]->setVal(1.0);
    ep_go[lev]->setVal(1.0);

    // Gas pressure fraction
    p_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    p_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    p_g[lev]->setVal(0.);
    p_go[lev]->setVal(0.);

    // Gas density
    ro_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    ro_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    ro_g[lev]->setVal(0.);
    ro_go[lev]->setVal(0.);

    // Gas bulk density
    rop_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    rop_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    rop_g[lev]->setVal(0.);
    rop_go[lev]->setVal(0.);

    // Base pressure that captures delp and/or p_in and p_out
    p0_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    p0_g[lev]->setVal(0.);

    // Pressure correction equation
    pp_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    pp_g[lev]->setVal(0.);

    // Molecular viscosity
    mu_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    mu_g[lev]->setVal(0.);

    //
    lambda_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    lambda_g[lev]->setVal(0.);
    //
    trD_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    trD_g[lev]->setVal(0.);

    // Vorticity
    vort[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    vort[lev]->setVal(0.);

    // Pressure increment
    phi[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    phi[lev]->setVal(0.);

    // diveu
    diveu[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    diveu[lev]->setVal(0.);

    // ********************************************************************************
    // X-face-based arrays
    // ********************************************************************************

    // Create a BoxArray on x-faces.
    BoxArray x_edge_ba = grids[lev];
    x_edge_ba.surroundingNodes(0);

    // X-axis gas velocity
    u_g[lev].reset(new MultiFab(x_edge_ba,dmap[lev],1,nghost));
    u_go[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    u_gt[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    fp_x[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    bcoeff[lev][0].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    u_g[lev]->setVal(0.);
    u_go[lev]->setVal(0.);
    u_gt[lev]->setVal(0.);
    fp_x[lev]->setVal(0.);
    bcoeff[lev][0]->setVal(0.);

    d_e[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    d_e[lev]->setVal(0.);

    tau_u_g[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    tau_u_g[lev]->setVal(0.);

    fluxX[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    fluxX[lev]->setVal(0.);

    ropX[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    ropX[lev]->setVal(0.);

    f_gds_u[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,1));
    f_gds_u[lev]->setVal(0.);

    drag_u[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,1));
    drag_u[lev]->setVal(0.);

    // u-velocity slopes. Note that the number of components is not 1, but 3!
    slopes_u[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],3,nghost));
    slopes_u[lev] -> setVal(0.);

    // u acceleration terms
    uacc[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    uacc[lev] -> setVal(0.);

    // ********************************************************************************
    // Y-face-based arrays
    // ********************************************************************************

    // Create a BoxArray on y-faces.
    BoxArray y_edge_ba = grids[lev];
    y_edge_ba.surroundingNodes(1);

    // Y-axis gas velocity
    v_g[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    v_go[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    v_gt[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    fp_y[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    bcoeff[lev][1].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    v_g[lev]->setVal(0.);
    v_go[lev]->setVal(0.);
    v_gt[lev]->setVal(0.);
    fp_y[lev]->setVal(0.);
    bcoeff[lev][1]->setVal(0.);

    d_n[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost));
    d_n[lev]->setVal(0.);

    tau_v_g[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost));
    tau_v_g[lev]->setVal(0.);

    fluxY[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost));
    fluxY[lev]->setVal(0.);

    ropY[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost));
    ropY[lev]->setVal(0.);

    f_gds_v[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,1));
    f_gds_v[lev]->setVal(0.);

    drag_v[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,1));
    drag_v[lev]->setVal(0.);

    // v-velocity slopes. Note that the number of components is not 1, but 3!
    slopes_v[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],3,nghost));
    slopes_v[lev] -> setVal(0.);

    // v acceleration terms
    vacc[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    vacc[lev] -> setVal(0.);



    // ********************************************************************************
    // Z-face-based arrays
    // ********************************************************************************

    // Create a BoxArray on y-faces.
    BoxArray z_edge_ba = grids[lev];
    z_edge_ba.surroundingNodes(2);

    // Z-axis gas velocity
    w_g[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    w_go[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    w_gt[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    fp_z[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    bcoeff[lev][2].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    w_g[lev]->setVal(0.);
    w_go[lev]->setVal(0.);
    w_gt[lev]->setVal(0.);
    fp_z[lev]->setVal(0.);
    bcoeff[lev][2]->setVal(0.);

    d_t[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    d_t[lev]->setVal(0.);

    tau_w_g[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    tau_w_g[lev]->setVal(0.);

    fluxZ[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    fluxZ[lev]->setVal(0.);

    ropZ[lev].reset(new MultiFab(z_edge_ba,dmap[lev],1,nghost));
    ropZ[lev]->setVal(0.);

    f_gds_w[lev].reset(new MultiFab(z_edge_ba,dmap[lev],1,1));
    f_gds_w[lev]->setVal(0.);

    drag_w[lev].reset(new MultiFab(z_edge_ba,dmap[lev],1,1));
    drag_w[lev]->setVal(0.);

    // w-velocity slopes. Note that the number of components is not 1, but 3!
    slopes_w[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],3,nghost));
    slopes_w[lev] -> setVal(0.);

    // w acceleration terms
    wacc[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    wacc[lev] -> setVal(0.);

}


void
mfix_level::InitLevelData(int lev, Real dt, Real time)
{

  // Allocate the fluid data
  if (solve_fluid)
     AllocateArrays(lev);

  // Allocate the particle data
  if (solve_dem)
  {
    //int lev = 0;
    pc -> AllocData();

    if (particle_init_type == "AsciiFile")
      {
        amrex::Print() << "Reading particles from particle_input.dat ..." << std::endl;
        pc -> InitParticlesAscii("particle_input.dat");

      } else if (particle_init_type == "Random")
      {
        int n_per_cell = 1;
        amrex::Print() << "Randomly initializing " << n_per_cell << " particles per cell ..." << std::endl;
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
        int phase = 1;
        int state = 0;
        MFIXParticleContainer::ParticleInitData pdata = {radius,volume,mass,density,omoi,
                velx,vely,velz,omegax,omegay,omegaz,dragx,dragy,dragz,phase,state};
        pc->InitNRandomPerCell(n_per_cell, pdata);
        pc->WriteAsciiFileForInit ("random_particles");
        exit(0);

      } else if (particle_init_type == "Auto") {

        amrex::Print() << "Auto generating particles ..." << std::endl;

        pc -> InitParticlesAuto(lev);

      } else {

      amrex::Abort("Bad particle_init_type");
    }

    // used in load balancing
    if (load_balance_type == "KnapSack") {
        particle_cost[lev].reset(new MultiFab(pc->ParticleBoxArray(lev),
                                              pc->ParticleDistributionMap(lev), 1, 0));
        particle_cost[lev]->setVal(0.0);

        if (solve_fluid)
        {
           fluid_cost[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 0));
           fluid_cost[lev]->setVal(0.0);
        }
    }
  }
}

void mfix_level::PostInit(int lev, Real dt, Real time, int nstep, int restart_flag, Real stop_time,
                          int steady_state)
{
  if (solve_dem) {

    // Auto generated particles may be out of the domain. This call will remove them.
    // Note that this has to occur after the EB geometry is created.
    if (particle_init_type == "Auto" && !restart_flag && particle_ebfactory)
      pc -> RemoveOutOfRange(lev, particle_ebfactory.get());

    Real avg_dp[10], avg_ro[10];
    pc -> GetParticleAvgProp( lev, avg_dp, avg_ro );
    init_collision(avg_dp, avg_ro);
  }

  // Initial fluid arrays: pressure, velocity, density, viscosity
  if (solve_fluid)
     mfix_init_fluid(lev,restart_flag,stop_time,steady_state);

  // Call user-defined subroutine to set constants, check data, etc.
  if (call_udf) mfix_usr0();

  // Calculate the initial volume fraction
  if (solve_fluid)
  {
     mfix_calc_volume_fraction(lev,sum_vol_orig);
     Print() << "Setting original sum_vol to " << sum_vol_orig << std::endl;
  }
}

void
mfix_level::MakeBCArrays ()
{
    int nghost = 2;

    // Define and allocate the integer MultiFab that is the outside adjacent cells of the problem domain.
    Box domainx(geom[0].Domain());
    domainx.grow(1,nghost);
    domainx.grow(2,nghost);
    Box box_ilo = amrex::adjCellLo(domainx,0,1);
    Box box_ihi = amrex::adjCellHi(domainx,0,1);

    Box domainy(geom[0].Domain());
    domainy.grow(0,nghost);
    domainy.grow(2,nghost);
    Box box_jlo = amrex::adjCellLo(domainy,1,1);
    Box box_jhi = amrex::adjCellHi(domainy,1,1);

    Box domainz(geom[0].Domain());
    domainz.grow(0,nghost);
    domainz.grow(1,nghost);
    Box box_klo = amrex::adjCellLo(domainz,2,1);
    Box box_khi = amrex::adjCellHi(domainz,2,1);

    // Note that each of these is a single IArrayBox so every process has a copy of them
    // No matter how thick the layer of BC nodes is, always resize with a factor of 2
    // 2 is the number on component (cfr the Fortran subroutines which handle BCs)
    bc_ilo.resize(box_ilo,2);
    bc_ihi.resize(box_ihi,2);
    bc_jlo.resize(box_jlo,2);
    bc_jhi.resize(box_jhi,2);
    bc_klo.resize(box_klo,2);
    bc_khi.resize(box_khi,2);
}

void
mfix_level::mfix_init_fluid(int lev, int is_restarting, Real stop_time, int steady_state)
{
  Box domain(geom[lev].Domain());

  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
  Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
  Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

  // Here we set bc values for p and u,v,w before the IC's are set
  mfix_set_bc0(lev);

  int delp_dir;

  // We deliberately don't tile this loop since we will be looping
  //    over bc's on faces and it makes more sense to do this one grid at a time
  for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.validbox();
    const Box& sbx = (*ep_g[lev])[mfi].box();

    if ( is_restarting ) {

      init_fluid_restart(sbx.loVect(), sbx.hiVect(), bx.loVect(),  bx.hiVect(),
           (*mu_g[lev])[mfi].dataPtr(), (*lambda_g[lev])[mfi].dataPtr());

    } else {
      const Box& ubx = (*u_g[lev])[mfi].box();
      const Box& vbx = (*v_g[lev])[mfi].box();
      const Box& wbx = (*w_g[lev])[mfi].box();

      init_fluid(sbx.loVect(), sbx.hiVect(),
           ubx.loVect(), ubx.hiVect(),
           vbx.loVect(), vbx.hiVect(),
           wbx.loVect(), wbx.hiVect(),
           bx.loVect(),  bx.hiVect(),
           domain.loVect(), domain.hiVect(),
           (*ep_g[lev])[mfi].dataPtr(),     (*ro_g[lev])[mfi].dataPtr(),
           (*rop_g[lev])[mfi].dataPtr(),    (*p_g[lev])[mfi].dataPtr(), (*p0_g[lev])[mfi].dataPtr(),
           (*u_g[lev])[mfi].dataPtr(),     (*v_g[lev])[mfi].dataPtr(),
           (*w_g[lev])[mfi].dataPtr(),
           (*mu_g[lev])[mfi].dataPtr(),   (*lambda_g[lev])[mfi].dataPtr(),
           &dx, &dy, &dz, &xlen, &ylen, &zlen);
    }

      set_p0(sbx.loVect(), sbx.hiVect(), bx.loVect(),  bx.hiVect(),
             domain.loVect(), domain.hiVect(),
             (*p0_g[lev])[mfi].dataPtr(),
             &dx, &dy, &dz, &xlen, &ylen, &zlen, &delp_dir);
  }

 // Here we set a separate periodicity flag for p0_g because when we use
 // pressure drop (delp) boundary conditions we fill all variables *except* p0
 // periodically
  IntVect press_per = IntVect(geom[lev].isPeriodic(0),geom[lev].isPeriodic(1),geom[lev].isPeriodic(2));
  if (delp_dir > -1) press_per[delp_dir] = 0;
  p0_periodicity = Periodicity(press_per);

  p0_g[lev]->FillBoundary(p0_periodicity);

  // Here we re-set the bc values for p and u,v,w just in case init_fluid
  //      over-wrote some of the bc values with ic values
  mfix_set_bc0(lev);

  // We deliberately don't tile this loop since we will be looping
  //    over bc's on faces and it makes more sense to do this one grid at a time
  if ( !is_restarting ) {

    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi) {

      const Box& sbx = (*ep_g[lev])[mfi].box();
      const Box& ubx = (*u_g[lev])[mfi].box();
      const Box& vbx = (*v_g[lev])[mfi].box();
      const Box& wbx = (*w_g[lev])[mfi].box();

      zero_wall_norm_vel(sbx.loVect(), sbx.hiVect(),
              ubx.loVect(), ubx.hiVect(),
              vbx.loVect(), vbx.hiVect(),
              wbx.loVect(), wbx.hiVect(),
              (*u_g[lev])[mfi].dataPtr(),
              (*v_g[lev])[mfi].dataPtr(),
              (*w_g[lev])[mfi].dataPtr(),
              bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
              bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect());
    }
  }

  if ( use_proj_method ) mfix_extrap_pressure(lev,p0_g[lev]);

  fill_mf_bc(lev,*ep_g[lev]);
  fill_mf_bc(lev,*ro_g[lev]);
  fill_mf_bc(lev,*rop_g[lev]);

  u_g[lev]->FillBoundary(geom[lev].periodicity());
  v_g[lev]->FillBoundary(geom[lev].periodicity());
  w_g[lev]->FillBoundary(geom[lev].periodicity());

  fill_mf_bc(lev,*mu_g[lev]);
  fill_mf_bc(lev,*lambda_g[lev]);

  if ( use_proj_method )
  {
     // We need to initialize the volume fraction ep_g before the first projection
     mfix_calc_volume_fraction(lev,sum_vol_orig);
     mfix_project_velocity(lev);
     mfix_initial_iterations(lev,stop_time,steady_state);
  }
}

void
mfix_level::mfix_set_bc0(int lev)
{
  Box domain(geom[lev].Domain());

  // Don't tile this -- at least for now
  for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
      const Box& sbx = (*ep_g[lev])[mfi].box();

      Box ubx((*u_g[lev])[mfi].box());
      Box vbx((*v_g[lev])[mfi].box());
      Box wbx((*w_g[lev])[mfi].box());

      set_bc0(sbx.loVect(), sbx.hiVect(),
              ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
              (*u_g[lev])[mfi].dataPtr(),     (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
              (*p0_g[lev])[mfi].dataPtr(),     (*ep_g[lev])[mfi].dataPtr(),
              (*ro_g[lev])[mfi].dataPtr(), (*rop_g[lev])[mfi].dataPtr(),
              (*mu_g[lev])[mfi].dataPtr(), (*lambda_g[lev])[mfi].dataPtr(),
              bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
              bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect());
    }

  fill_mf_bc(lev,*p_g[lev]);
  fill_mf_bc(lev,*ep_g[lev]);
  fill_mf_bc(lev,*ro_g[lev]);
  fill_mf_bc(lev,*rop_g[lev]);

  u_g[lev]->FillBoundary(geom[lev].periodicity());
  v_g[lev]->FillBoundary(geom[lev].periodicity());
  w_g[lev]->FillBoundary(geom[lev].periodicity());
}

void mfix_level::WriteEBSurface(int lev) {
  if (Geom(0).isAllPeriodic()) return;

  const Real* dx = Geom(lev).CellSize();

  BoxArray ba = grids[lev];

  // This creates the associated Distribution Mapping
  // DistributionMapping dm(ba, ParallelDescriptor::NProcs());

  MultiFab dummy(ba, dmap[lev], 1, 0, MFInfo(), *ebfactory);

  // // // Deliberately didn't time this loop.
  for (MFIter mfi(dummy); mfi.isValid(); ++mfi) {

    const auto& sfab = dynamic_cast<EBFArrayBox const&>((dummy)[mfi]);
    const auto& flag = sfab.getEBCellFlagFab();

    const Box& bx = mfi.validbox();

    if (flag.getType(bx) == FabType::covered or flag.getType(bx) == FabType::regular) continue;

    std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
    const MultiCutFab* bndrycent;

    areafrac  =  ebfactory->getAreaFrac();
    bndrycent = &(ebfactory->getBndryCent());

    mfix_eb_to_polygon(dx, bx.loVect(), bx.hiVect(),
         flag.dataPtr(), flag.loVect(), flag.hiVect(),
         (*bndrycent)[mfi].dataPtr(),
         (*bndrycent)[mfi].loVect(), (*bndrycent)[mfi].hiVect(),
         (*areafrac[0])[mfi].dataPtr(),
         (*areafrac[0])[mfi].loVect(), (*areafrac[0])[mfi].hiVect(),
         (*areafrac[1])[mfi].dataPtr(),
         (*areafrac[1])[mfi].loVect(), (*areafrac[1])[mfi].hiVect(),
         (*areafrac[2])[mfi].dataPtr(),
         (*areafrac[2])[mfi].loVect(), (*areafrac[2])[mfi].hiVect());
  }

  int cpu = ParallelDescriptor::MyProc();
  mfix_write_eb_vtp(&cpu);

  int nProcs = ParallelDescriptor::NProcs();
  if(ParallelDescriptor::IOProcessor())
    mfix_write_pvtp(&nProcs);

}
