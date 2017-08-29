#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

void
mfix_level::InitParams(int solve_fluid_in, int solve_dem_in,
                       int max_nit_in, int call_udf_in)
{
    solve_fluid  = solve_fluid_in;
    solve_dem    = solve_dem_in;
    max_nit      = max_nit_in;
    call_udf     = call_udf_in;
}

void mfix_level::Init(int lev, Real dt, Real time)
{
    BL_ASSERT(max_level == 0);

    // Define coarse level BoxArray and DistributionMap
    finest_level = 0;

    const BoxArray& ba = MakeBaseGrids();
    DistributionMapping dm(ba, ParallelDescriptor::NProcs());

    MakeNewLevelFromScratch(0, time, ba, dm);

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

    InitIOData ();

    mfix_set_bc_type(lev);
}

BoxArray
mfix_level::MakeBaseGrids () const
{
    BoxArray ba(geom[0].Domain());
    ba.maxSize(max_grid_size[0]);
    if (refine_grid_layout) {
        ChopGrids(0, ba, ParallelDescriptor::NProcs());
    }
    if (ba == grids[0]) {
        ba = grids[0];  // to avoid dupliates
    }
    if ( ParallelDescriptor::IOProcessor() )
        std::cout << "In MakeBaseGrids: BA HAS " << ba.size() << " GRIDS " << std::endl;
    return ba;
}

void
mfix_level::MakeNewLevelFromScratch (int lev, Real time,
             const BoxArray& new_grids, const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

#if 0
    t_new[lev] = time;
    t_old[lev] = time - 1.e200;
#endif

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
    bc_ilo.resize(box_ilo,nghost_bc);
    bc_ihi.resize(box_ihi,nghost_bc);
    bc_jlo.resize(box_jlo,nghost_bc);
    bc_jhi.resize(box_jhi,nghost_bc);
    bc_klo.resize(box_klo,nghost_bc);
    bc_khi.resize(box_khi,nghost_bc);

    check_data(lev);

    Box domain(geom[0].Domain());

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    // This is separate from check_data because it is only called on initialization,
    // not on restart
    if ( ParallelDescriptor::IOProcessor() )
       check_initial_conditions(&dx,&dy,&dz,domain.loVect(),domain.hiVect());
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
       check_initial_conditions(&dx,&dy,&dz,domain.loVect(),domain.hiVect());
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
    u_g[lev]->setVal(0.);
    u_go[lev]->setVal(0.);
    u_gt[lev]->setVal(0.);

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
    v_g[lev]->setVal(0.);
    v_go[lev]->setVal(0.);
    v_gt[lev]->setVal(0.);

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
    w_g[lev]->setVal(0.);
    w_go[lev]->setVal(0.);
    w_gt[lev]->setVal(0.);

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
}

void
mfix_level::InitLevelData(int lev, Real dt, Real time)
{
  AllocateArrays(lev);

  mfix_set_bc0(lev);

  // Initial fluid arrays: pressure, velocity, density, viscosity
  mfix_init_fluid(lev);

  // Allocate the particle arrays
  if (solve_dem)
    {
      pc -> AllocData();
      pc -> InitParticlesAscii("particle_input.dat");

      Real avg_dp[10], avg_ro[10];
      pc -> GetParticleAvgProp( lev, avg_dp, avg_ro );

      init_collision(avg_dp, avg_ro);

      //  Create mask for particle ghost cells
      pc -> InitLevelMask( lev, geom[lev], dmap[lev], grids[lev] );

      pc -> BalanceParticleLoad_KDTree ();
    }

  if (solve_dem)
  {
     mfix_calc_volume_fraction(lev,sum_vol_orig);
     Print() << "Setting original sum_vol to " << sum_vol_orig << std::endl;
  }

  // Call user-defined subroutine to set constants, check data, etc.
  if (call_udf) mfix_usr0();

  // Calculate all the coefficients once before entering the time loop
  int calc_flag = 2;
  mfix_calc_coeffs(lev,calc_flag);
}

void
mfix_level::InitLevelDataFromRestart(int lev, Real dt, Real time)
{

  if (solve_dem) {
      Real avg_dp[10], avg_ro[10];
      pc -> GetParticleAvgProp( lev, avg_dp, avg_ro );
      init_collision(avg_dp, avg_ro);
      pc -> BalanceParticleLoad_KDTree ();
  }

  mfix_set_bc0(lev);

  // Initial fluid arrays: pressure, velocity, density, viscosity
  mfix_init_fluid(lev,1);

  // Call user-defined subroutine to set constants, check data, etc.
  // if (call_udf)
  //  usr0();

  // Calculate all the coefficients once before entering the time loop
  int calc_flag = 2;
  mfix_calc_coeffs(lev,calc_flag);

  if (solve_dem)
  {
     mfix_calc_volume_fraction(lev,sum_vol_orig);
     Print() << "Setting original sum_vol to " << sum_vol_orig << std::endl;
  }
}

void
mfix_level::mfix_init_fluid(int lev, int is_restarting)
{
  Box domain(geom[lev].Domain());

  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
  Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
  Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

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
           (*rop_g[lev])[mfi].dataPtr(),     (*p_g[lev])[mfi].dataPtr(),
           (*u_g[lev])[mfi].dataPtr(),     (*v_g[lev])[mfi].dataPtr(),
           (*w_g[lev])[mfi].dataPtr(),
           (*mu_g[lev])[mfi].dataPtr(),   (*lambda_g[lev])[mfi].dataPtr(),
           &dx, &dy, &dz, &xlen, &ylen, &zlen );
    }
  }

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

  fill_mf_bc(lev,*p_g[lev]);
  fill_mf_bc(lev,*ep_g[lev]);
  fill_mf_bc(lev,*ro_g[lev]);
  fill_mf_bc(lev,*rop_g[lev]);

  u_g[lev]->FillBoundary(geom[lev].periodicity());
  v_g[lev]->FillBoundary(geom[lev].periodicity());
  w_g[lev]->FillBoundary(geom[lev].periodicity());

  fill_mf_bc(lev,*mu_g[lev]);
  fill_mf_bc(lev,*lambda_g[lev]);

}


void
mfix_level::mfix_set_bc0(int lev)
{

  // Array already allocated when restart is called
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
              (*p_g[lev])[mfi].dataPtr(),     (*ep_g[lev])[mfi].dataPtr(),
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
