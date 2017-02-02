#include <ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <BC_TYPES.H>

mfix_level::~mfix_level ()
{};

mfix_level::mfix_level (const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord)
 : AmrCore(rb,max_level_in,n_cell_in,coord)
{
    ReadParameters();

    // Geometry on all levels has just been defined in the AmrCore constructor

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;
#if 0
    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= maxLevel(); ++lev) {
       nsubsteps[lev] = MaxRefRatio(lev-1);
    }
#endif

    // Particle Container
    mypc = std::unique_ptr<MyParticleContainer> (new MyParticleContainer(this));

    flag.resize(nlevs_max);

    A_m.resize(nlevs_max);
    b_m.resize(nlevs_max);

    ep_g.resize(nlevs_max);
    ep_go.resize(nlevs_max);

    p_g.resize(nlevs_max);
    p_go.resize(nlevs_max);

    pp_g.resize(nlevs_max);

    ro_g.resize(nlevs_max);
    ro_go.resize(nlevs_max);

    rop_g.resize(nlevs_max);
    rop_go.resize(nlevs_max);

    u_g.resize(nlevs_max);
    u_go.resize(nlevs_max);
    u_gt.resize(nlevs_max);

    v_g.resize(nlevs_max);
    v_go.resize(nlevs_max);
    v_gt.resize(nlevs_max);

    w_g.resize(nlevs_max);
    w_go.resize(nlevs_max);
    w_gt.resize(nlevs_max);

    d_e.resize(nlevs_max);
    d_n.resize(nlevs_max);
    d_t.resize(nlevs_max);

    mu_g.resize(nlevs_max);
    lambda_g.resize(nlevs_max);
    trD_g.resize(nlevs_max);

    flux_gE.resize(nlevs_max);
    flux_gN.resize(nlevs_max);
    flux_gT.resize(nlevs_max);

    rop_gE.resize(nlevs_max);
    rop_gN.resize(nlevs_max);
    rop_gT.resize(nlevs_max);

    tau_u_g.resize(nlevs_max);
    tau_v_g.resize(nlevs_max);
    tau_w_g.resize(nlevs_max);

    f_gds.resize(nlevs_max);
    drag_bm.resize(nlevs_max);

    int nparticles = 100;

    particle_state.resize(  nparticles);
    particle_phase.resize(  nparticles);

    des_radius.resize    (  nparticles);
    ro_sol.resize        (  nparticles);
    pvol.resize          (  nparticles);
    pmass.resize         (  nparticles);
    omoi.resize          (  nparticles);
    des_pos_new.resize   (3*nparticles);
    des_vel_new.resize   (3*nparticles);
    des_usr_var.resize   (  nparticles);
    omega_new.resize     (3*nparticles);
    des_acc_old.resize   (3*nparticles);
    rot_acc_old.resize   (3*nparticles);
    drag_fc.resize       (3*nparticles);
    fc.resize            (3*nparticles);
    tow.resize           (3*nparticles);
    pairs.resize         (12*nparticles);
}

void
mfix_level::ReadParameters ()
{
    {
    ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
    pp.query("max_step", max_step);
    pp.query("stop_time", stop_time);
    }

    {
    ParmParse pp("amr"); // Traditionally, these have prefix, amr.
    pp.query("check_file", check_file);
    pp.query("check_int", check_int);
    pp.query("plot_file", plot_file);
    pp.query("plot_int", plot_int);
    pp.query("restart", restart_chkfile);
    }

    {
    ParmParse pp("mfix");
    pp.query("verbose", verbose);
    }
}

void
mfix_level::InitParams(int solve_fluid_in, int solve_dem_in, int cyclic_mf_in,
                       int max_nit_in, int call_udf_in)
{
   solve_fluid  = solve_fluid_in;
   solve_dem    = solve_dem_in;
   cyclic_mf    = cyclic_mf_in;
   max_nit      = max_nit_in;
   call_udf     = call_udf_in;
}

void
mfix_level::Init(int lev, Real dt, Real time)
{
    BL_ASSERT(max_level == 0);

    // define coarse level BoxArray and DistributionMap
    {
      finest_level = 0;

      const BoxArray& ba = MakeBaseGrids();
      DistributionMapping dm(ba, ParallelDescriptor::NProcs());

      MakeNewLevel(0, time, ba, dm);

      InitLevelData(lev,dt,time);
    }

    // if max_level > 0, define fine levels

    mypc->InitData();
}

void
mfix_level::Restart()
{
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
    std::cout << "BA " << ba << std::endl;
    return ba;
}

void
mfix_level::MakeNewLevel (int lev, Real time,
              const BoxArray& new_grids, const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

#if 0
    t_new[lev] = time;
    t_old[lev] = time - 1.e200;
#endif

    int nghost;
    if (ParallelDescriptor::NProcs() == 1) {
       nghost = 1;
    } else {
       nghost = 2;
    }

    // Define and allocate the integer MultiFab on BoxArray ba with
    // 4 components and nghost ghost cells.
    flag[lev].reset(new iMultiFab(grids[lev],4,nghost,dmap[lev],Fab_allocate));
    flag[lev]->setVal(0);

    // Call set_domain for each subdomain
    // Read input data, check data, do computations for IC and BC locations
    // and flows, and set geometry parameters such as X, X_E, DToDX, etc.

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
    {
       // const Box&  bx = mfi.validbox();
       const Box& sbx = (*flag[lev])[mfi].box();

       set_domain(sbx.loVect(),sbx.hiVect(),
                  (*flag[lev])[mfi].dataPtr(),&dx,&dy,&dz);
    }
    fill_mf_bc(lev,*flag[lev]);

    // Matrix and rhs vector
    A_m[lev].reset(new MultiFab(grids[lev],7,nghost,dmap[lev],Fab_allocate));
    b_m[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));

    // ********************************************************************************
    // Cell-based arrays
    // ********************************************************************************

    // Void fraction
    ep_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    ep_go[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    ep_g[lev]->setVal(1.0);
    ep_go[lev]->setVal(1.0);

    // Gas pressure fraction
    p_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    p_go[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    p_g[lev]->setVal(0.);
    p_go[lev]->setVal(0.);

    // Gas density
    ro_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    ro_go[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    ro_g[lev]->setVal(0.);
    ro_go[lev]->setVal(0.);

    // Gas bulk density
    rop_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    rop_go[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    rop_g[lev]->setVal(0.);
    rop_go[lev]->setVal(0.);

    // Pressure correction equation
    pp_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    pp_g[lev]->setVal(0.);

    // Molecular viscosity
    mu_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    mu_g[lev]->setVal(0.);

    //
    lambda_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    lambda_g[lev]->setVal(0.);

    //
    trD_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    trD_g[lev]->setVal(0.);

    //
    f_gds[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    f_gds[lev]->setVal(0.);

    //
    drag_bm[lev].reset(new MultiFab(grids[lev],3,nghost,dmap[lev],Fab_allocate));
    drag_bm[lev]->setVal(0.);

    // ********************************************************************************
    // X-face-based arrays
    // ********************************************************************************

    // X-axis gas velocity
    u_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    u_go[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    u_gt[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    u_g[lev]->setVal(0.);
    u_go[lev]->setVal(0.);
    u_gt[lev]->setVal(0.);

    d_e[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    d_e[lev]->setVal(0.);

    tau_u_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    tau_u_g[lev]->setVal(0.);

    flux_gE[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    flux_gE[lev]->setVal(0.);

    rop_gE[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    rop_gE[lev]->setVal(0.);

    // ********************************************************************************
    // Y-face-based arrays
    // ********************************************************************************

    // Y-axis gas velocity
    v_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    v_go[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    v_gt[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    v_g[lev]->setVal(0.);
    v_go[lev]->setVal(0.);
    v_gt[lev]->setVal(0.);

    d_n[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    d_n[lev]->setVal(0.);

    tau_v_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    tau_v_g[lev]->setVal(0.);

    flux_gN[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    flux_gN[lev]->setVal(0.);

    rop_gN[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    rop_gN[lev]->setVal(0.);

    // ********************************************************************************
    // Z-face-based arrays
    // ********************************************************************************

    // Z-axis gas velocity
    w_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    w_go[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    w_gt[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    w_g[lev]->setVal(0.);
    w_go[lev]->setVal(0.);
    w_gt[lev]->setVal(0.);

    d_t[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    d_t[lev]->setVal(0.);

    tau_w_g[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    tau_w_g[lev]->setVal(0.);

    flux_gT[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    flux_gT[lev]->setVal(0.);

    rop_gT[lev].reset(new MultiFab(grids[lev],1,nghost,dmap[lev],Fab_allocate));
    rop_gT[lev]->setVal(0.);

    // ********************************************************************************
}

void
mfix_level::evolve_fluid(int lev, int nstep, int set_normg,
                         Real dt, Real& prev_dt, Real time, Real normg)
{
      Real dx = geom[lev].CellSize(0);
      Real dy = geom[lev].CellSize(1);
      Real dz = geom[lev].CellSize(2);

      // Calculate transport coefficients
      mfix_calc_all_coeffs(lev);

      // Calculate the stress tensor trace and cross terms for all phases.
      mfix_calc_trd_and_tau(lev);

      // Backup field variable to old
      int nghost = ep_go[lev]->nGrow();
      MultiFab::Copy(*ep_go[lev],  *ep_g[lev],  0, 0, 1, nghost);
      MultiFab::Copy(*p_go[lev],   *p_g[lev],   0, 0, 1, nghost);
      MultiFab::Copy(*ro_go[lev],  *ro_g[lev],  0, 0, 1, nghost);
      MultiFab::Copy(*rop_go[lev], *rop_g[lev], 0, 0, 1, nghost);
      MultiFab::Copy(*u_go[lev],   *u_g[lev],   0, 0, 1, nghost);
      MultiFab::Copy(*v_go[lev],   *v_g[lev],   0, 0, 1, nghost);
      MultiFab::Copy(*w_go[lev],   *w_g[lev],   0, 0, 1, nghost);

      // Loop over iterate for auto time-step size adjustment
      int reiterate;
      do {
        prev_dt = dt;

        // Calculate bulk density (epg*ro_g) at cell faces
        mfix_conv_rop(lev,dt);

        // Calculate face mass fluxes
        mfix_calc_mflux(lev);

        int converged=0;
        int nit=0;          // number of iterations
        int gsmf=0;         // number of outer iterations for goal seek mass flux (GSMF)
        Real delP_MF=0.0L;  // actual GSMF pressure drop
        Real lMFlux=0.0L;   // actual GSMF mass flux
        Real resg=0.0L;     // fluid pressure residual

        // int lset_normg=1-set_normg;
        Real lnormg=normg;

        ///////////////// ---- call to iterate -------- /////////////////
        do {
          nit++;

          // User hooks
          for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
            mfix_usr2();

          // Calculate transport coefficients
          int calc_flag = 1;
          mfix_calc_coeffs(lev,calc_flag);

          mfix_solve_for_vels(lev,dt);

          // Calculate transport coefficients
          mfix_physical_prop(lev,0);

          // Calculate bulk density (epg*ro_g) at cell faces
          mfix_conv_rop(lev,dt);

          // Solve the pressure correction equation
          mfix_solve_for_pp(lev,dt,lnormg,resg);

          mfix_correct0(lev);

          // Update fluid density
          mfix_physical_prop(lev,0);

          // Update wall velocities
          mfix_set_wall_bc(lev);

          // Calculate face mass fluxes
          mfix_calc_mflux(lev);

          // Display current iteration residuals
          display_resid(&nit);

          // Check for convergence
          converged = check_convergence(&nit);

          // Iterate over cyclic mass flux bc
          if(cyclic_mf==1 && (converged==1 || nit >= max_nit))
            for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
            {
              const Box& sbx = (*flag[lev])[mfi].box();

              converged = goal_seek_mFlux(sbx.loVect(), sbx.hiVect(), &nit, &gsmf, &delP_MF, &lMFlux,
                (*flux_gE[lev])[mfi].dataPtr(),  (*flux_gN[lev])[mfi].dataPtr(),  (*flux_gT[lev])[mfi].dataPtr(),
                (*flag[lev])[mfi].dataPtr(), &dx, &dy, &dz);
            }

        } while(converged==0 && nit<max_nit);

        // Adjust time step if iteration failed.
        reiterate = mfix_adjustdt(&converged, &nit, &dt);
        if(reiterate == 1) {

          // Reset the field variables
          MultiFab::Copy(*ep_g[lev],  *ep_go[lev],  0, 0, 1, nghost);
          MultiFab::Copy(*p_g[lev],   *p_go[lev],   0, 0, 1, nghost);
          MultiFab::Copy(*ro_g[lev],  *ro_go[lev],  0, 0, 1, nghost);
          MultiFab::Copy(*rop_g[lev], *rop_go[lev], 0, 0, 1, nghost);
          MultiFab::Copy(*u_g[lev],   *u_go[lev],   0, 0, 1, nghost);
          MultiFab::Copy(*v_g[lev],   *v_go[lev],   0, 0, 1, nghost);
          MultiFab::Copy(*w_g[lev],   *w_go[lev],   0, 0, 1, nghost);

          // Recalculate all coefficients (JM: not sure why)
          mfix_calc_all_coeffs(lev);
        }
      } while (reiterate==1);
}

void
mfix_level::evolve_dem(int lev, int nstep, Real dt, Real time)
{
    int pair_count = 0;

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    const int max_pip = particle_state.size();

    for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
    {
       // const Box& bx = mfi.validbox();
       const Box& sbx = (*flag[lev])[mfi].box();

       mfix_des_time_march(sbx.loVect(), sbx.hiVect(), &max_pip,
        (*ep_g[lev])[mfi].dataPtr(), (*p_g[lev])[mfi].dataPtr(),
        (*u_g[lev])[mfi].dataPtr(),  (*v_g[lev])[mfi].dataPtr(), (*w_g[lev])[mfi].dataPtr(),
        (*ro_g[lev])[mfi].dataPtr(),
        (*mu_g[lev])[mfi].dataPtr(),
        particle_state.dataPtr(), particle_phase.dataPtr(),
        des_radius.dataPtr(),
        pvol.dataPtr(),           pmass.dataPtr(),
        omoi.dataPtr(),           des_usr_var.dataPtr(),
        des_pos_new.dataPtr(),    des_vel_new.dataPtr(),   omega_new.dataPtr(),
        des_acc_old.dataPtr(),    rot_acc_old.dataPtr(),
        drag_fc.dataPtr(),        fc.dataPtr(),            tow.dataPtr(),
        pairs.dataPtr(),          &pair_count,
        (*flag[lev])[mfi].dataPtr(),
        &time, &dt, &dx, &dy, &dz, &nstep);
    }

    fill_mf_bc(lev,*ep_g[lev]);
    fill_mf_bc(lev,*rop_g[lev]);
}

void
mfix_level::output(int lev, int estatus, int finish, int nstep, Real dt, Real time)
{
  const int max_pip = particle_state.size();

  for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
  {
     // const Box& sbx = (*flag[lev])[mfi].box();
     mfix_output_manager(&max_pip,
      &time, &dt, &nstep,
      particle_state.dataPtr(), des_radius.dataPtr(),
      des_pos_new.dataPtr(),
      des_vel_new.dataPtr(), des_usr_var.dataPtr(),
      omega_new.dataPtr(), &finish);
  }
}

void
mfix_level::InitLevelData(int lev, Real dt, Real time)
{
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
  {
     const Box& bx = mfi.validbox();
     const Box& sbx = (*flag[lev])[mfi].box();

     mfix_main1(sbx.loVect(), sbx.hiVect(),
                 bx.loVect(),  bx.hiVect(),
               &time, &dt,
               (*u_g[lev])[mfi].dataPtr(),     (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
               (*p_g[lev])[mfi].dataPtr(),     (*ep_g[lev])[mfi].dataPtr(),
               (*flag[lev])[mfi].dataPtr(), &dx, &dy, &dz );
  }

  fill_mf_bc(lev,*p_g[lev]);
  fill_mf_bc(lev,*ep_g[lev]);
  fill_mf_bc(lev,*ro_g[lev]);
  fill_mf_bc(lev,*rop_g[lev]);

  fill_mf_bc(lev,*u_g[lev]);
  fill_mf_bc(lev,*v_g[lev]);
  fill_mf_bc(lev,*w_g[lev]);

  // Allocate the particle arrays
  if (solve_dem)
  {
    for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi) {
      const int max_pip = particle_state.size();

      mfix_make_arrays_des(&max_pip,
        particle_state.dataPtr(),
        particle_phase.dataPtr(), des_radius.dataPtr(), ro_sol.dataPtr(),
        pvol.dataPtr(), pmass.dataPtr(), omoi.dataPtr(),
        des_pos_new.dataPtr(), des_vel_new.dataPtr(),
        des_usr_var.dataPtr(), omega_new.dataPtr(),
        fc.dataPtr(), tow.dataPtr());
    }

    for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi){
      const int max_pip = particle_state.size();
      mfix_write_des_data(&max_pip, particle_state.dataPtr(), des_radius.dataPtr(),
        des_pos_new.dataPtr(), des_vel_new.dataPtr(), des_usr_var.dataPtr());
    }
  }

  // Calculate volume fraction, ep_g.
  if (solve_dem)
    mfix_comp_mean_fields(lev);

  // Initial fluid arrays: pressure, velocity, density, viscosity
  mfix_init_fluid(lev);

  // Call user-defined subroutine to set constants, check data, etc.
  if (call_udf)
      mfix_usr0();

  // Calculate all the coefficients once before entering the time loop
  int calc_flag = 2;
  mfix_calc_coeffs(lev,calc_flag);

  mfix_finl_err_msg();
}

void
mfix_level::mfix_calc_coeffs(int lev, int calc_flag)
{
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
  {
     const Box& bx = mfi.validbox();
     const Box& sbx = (*flag[lev])[mfi].box();
     const int max_pip = particle_state.size();

     calc_coeff(sbx.loVect(), sbx.hiVect(),
                 &max_pip,
      (*flag[lev])[mfi].dataPtr(),    &calc_flag,
      (*ro_g[lev])[mfi].dataPtr(),    (*p_g[lev])[mfi].dataPtr(),
      (*ep_g[lev])[mfi].dataPtr(),    (*rop_g[lev])[mfi].dataPtr(),
      (*u_g[lev])[mfi].dataPtr(),     (*v_g[lev])[mfi].dataPtr(),
      (*w_g[lev])[mfi].dataPtr(),     (*mu_g[lev])[mfi].dataPtr(),
      (*f_gds[lev])[mfi].dataPtr(),   (*drag_bm[lev])[mfi].dataPtr(),
      particle_phase.dataPtr(), particle_state.dataPtr(),
      pvol.dataPtr(), des_pos_new.dataPtr(), des_vel_new.dataPtr(),
      des_radius.dataPtr(), &dx, &dy, &dz );
  }
  fill_mf_bc(lev,*ro_g[lev]);
  fill_mf_bc(lev,*rop_g[lev]);

  if (solve_dem)
  {
    fill_mf_bc(lev,*f_gds[lev]);
    fill_mf_bc(lev,*drag_bm[lev]);
  }
}

void
mfix_level::mfix_calc_all_coeffs(int lev)
{
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  const int max_pip = particle_state.size();

  for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
  {
     // const Box& bx = mfi.validbox();
     const Box& sbx = (*flag[lev])[mfi].box();

     calc_coeff_all(sbx.loVect(), sbx.hiVect(), &max_pip,
       (*ro_g[lev])[mfi].dataPtr(), (*p_g[lev])[mfi].dataPtr(),
       (*ep_g[lev])[mfi].dataPtr(), (*rop_g[lev])[mfi].dataPtr(),
       (*u_g[lev])[mfi].dataPtr(),  (*v_g[lev])[mfi].dataPtr(),   (*w_g[lev])[mfi].dataPtr(),
       (*mu_g[lev])[mfi].dataPtr(), (*f_gds[lev])[mfi].dataPtr(), (*drag_bm[lev])[mfi].dataPtr(),
       particle_phase.dataPtr(),  particle_state.dataPtr(),
       pvol.dataPtr(), des_pos_new.dataPtr(),
       des_vel_new.dataPtr(), des_radius.dataPtr(),
       (*flag[lev])[mfi].dataPtr(), &dx, &dy, &dz);
  }
  fill_mf_bc(lev,*ro_g[lev]);
  fill_mf_bc(lev,*rop_g[lev]);

  if (solve_dem)
  {
    fill_mf_bc(lev,*ep_g[lev]);
    fill_mf_bc(lev,*f_gds[lev]);
    fill_mf_bc(lev,*drag_bm[lev]);
  }
}

void
mfix_level::mfix_calc_trd_and_tau(int lev)
{
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
  {
     const Box& bx = mfi.validbox();
     const Box& sbx = (*flag[lev])[mfi].box();

     calc_trd_and_tau(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),
       (*tau_u_g[lev])[mfi].dataPtr(),  (*tau_v_g[lev])[mfi].dataPtr(), (*tau_w_g[lev])[mfi].dataPtr(),
       (*trD_g[lev])[mfi].dataPtr(),    (*ep_g[lev])[mfi].dataPtr(),
       (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),     (*w_g[lev])[mfi].dataPtr(),
       (*lambda_g[lev])[mfi].dataPtr(), (*mu_g[lev])[mfi].dataPtr(),
       (*flag[lev])[mfi].dataPtr(), &dx, &dy, &dz);
  }

  fill_mf_bc(lev,*tau_u_g[lev]);
  fill_mf_bc(lev,*tau_v_g[lev]);
  fill_mf_bc(lev,*tau_w_g[lev]);
  fill_mf_bc(lev,*trD_g[lev]);
}

void
mfix_level::mfix_init_fluid(int lev)
{
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
  {
     const Box& bx = mfi.validbox();
     const Box& sbx = (*flag[lev])[mfi].box();

     init_fluid(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),
       (*ep_g[lev])[mfi].dataPtr(),     (*ro_g[lev])[mfi].dataPtr(),
       (*rop_g[lev])[mfi].dataPtr(),     (*p_g[lev])[mfi].dataPtr(),
       (*u_g[lev])[mfi].dataPtr(),     (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
       (*mu_g[lev])[mfi].dataPtr(),   (*lambda_g[lev])[mfi].dataPtr(),
       (*flag[lev])[mfi].dataPtr(), &dx, &dy, &dz );
  }

  fill_mf_bc(lev,*p_g[lev]);
  fill_mf_bc(lev,*ep_g[lev]);
  fill_mf_bc(lev,*ro_g[lev]);
  fill_mf_bc(lev,*rop_g[lev]);

  fill_mf_bc(lev,*u_g[lev]);
  fill_mf_bc(lev,*v_g[lev]);
  fill_mf_bc(lev,*w_g[lev]);

  fill_mf_bc(lev,*mu_g[lev]);
  fill_mf_bc(lev,*lambda_g[lev]);
}

void
mfix_level::mfix_comp_mean_fields(int lev)
{
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    const int max_pip = particle_state.size();

    for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.validbox();
       const Box& sbx = (*flag[lev])[mfi].box();

       comp_mean_fields(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),
            &max_pip, (*ep_g[lev])[mfi].dataPtr(),
            particle_state.dataPtr(), des_pos_new.dataPtr(), pvol.dataPtr(),
            (*flag[lev])[mfi].dataPtr(), &dx, &dy, &dz );
    }
    fill_mf_bc(lev,*ep_g[lev]);
}

void
mfix_level::mfix_calc_mflux(int lev)
{
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
  {
     const Box& bx = mfi.validbox();
     const Box& sbx = (*flag[lev])[mfi].box();

     calc_mflux(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),
       (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
       (*rop_gE[lev])[mfi].dataPtr(),   (*rop_gN[lev])[mfi].dataPtr(),   (*rop_gT[lev])[mfi].dataPtr(),
       (*flux_gE[lev])[mfi].dataPtr(),  (*flux_gN[lev])[mfi].dataPtr(),  (*flux_gT[lev])[mfi].dataPtr(),
       (*flag[lev])[mfi].dataPtr(),     &dx, &dy, &dz);
  }
  fill_mf_bc(lev,*flux_gE[lev]);
  fill_mf_bc(lev,*flux_gN[lev]);
  fill_mf_bc(lev,*flux_gT[lev]);
}


void
mfix_level::mfix_set_wall_bc(int lev)
{
  for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
  {
     const Box& bx = mfi.validbox();
     const Box& sbx = (*flag[lev])[mfi].box();

     set_wall_bc(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),
       (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
       (*flag[lev])[mfi].dataPtr());
  }
    fill_mf_bc(lev,*u_g[lev]);
    fill_mf_bc(lev,*v_g[lev]);
    fill_mf_bc(lev,*w_g[lev]);
}

void
mfix_level::mfix_conv_rop(int lev, Real dt)
{
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.validbox();
       const Box& sbx = (*flag[lev])[mfi].box();

       conv_rop(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),
         (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
         (*rop_g[lev])[mfi].dataPtr(),
         (*rop_gE[lev])[mfi].dataPtr(),   (*rop_gN[lev])[mfi].dataPtr(),   (*rop_gT[lev])[mfi].dataPtr(),
         (*flag[lev])[mfi].dataPtr(),     &dt, &dx, &dy, &dz);
    }

    fill_mf_bc(lev,*rop_gE[lev]);
    fill_mf_bc(lev,*rop_gN[lev]);
    fill_mf_bc(lev,*rop_gT[lev]);
}

void
mfix_level::mfix_solve_for_vels(int lev, Real dt)
{
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    // Solve U-Momentum equation
    MultiFab::Copy(*u_gt[lev], *u_g[lev], 0, 0, 1, u_g[lev]->nGrow());
    for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.validbox();
      const Box& sbx = (*flag[lev])[mfi].box();

      Box ubx((*u_g[lev])[mfi].box()); //ubx.shift(0,-1);
      Box vbx((*v_g[lev])[mfi].box()); //vbx.shift(1,-1);
      Box wbx((*w_g[lev])[mfi].box()); //wbx.shift(2,-1);

      solve_u_g_star(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),
          ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
          (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
          (*u_go[lev])[mfi].dataPtr(),     (*p_g[lev])[mfi].dataPtr(),      (*ro_g[lev])[mfi].dataPtr(),
          (*rop_g[lev])[mfi].dataPtr(),    (*rop_go[lev])[mfi].dataPtr(),   (*ep_g[lev])[mfi].dataPtr(),
          (*tau_u_g[lev])[mfi].dataPtr(),  (*d_e[lev])[mfi].dataPtr(),
          (*flux_gE[lev])[mfi].dataPtr(),  (*flux_gN[lev])[mfi].dataPtr(),  (*flux_gT[lev])[mfi].dataPtr(),
          (*mu_g[lev])[mfi].dataPtr(),     (*f_gds[lev])[mfi].dataPtr(),
                (*A_m[lev])[mfi].dataPtr(),      (*b_m[lev])[mfi].dataPtr(),      (*drag_bm[lev])[mfi].dataPtr(),
          (*flag[lev])[mfi].dataPtr(),     &dt, &dx, &dy, &dz);
    }

    fill_mf_bc(lev,*A_m[lev]);
    fill_mf_bc(lev,*b_m[lev]);
    fill_mf_bc(lev,*d_e[lev]);

    int eq_id=3;
    mfix_solve_linear_equation(eq_id,lev,(*u_gt[lev]),(*A_m[lev]),(*b_m[lev]));

    // Solve V-Momentum equation
    MultiFab::Copy(*v_gt[lev], *v_g[lev], 0, 0, 1, v_g[lev]->nGrow());
    for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.validbox();
      const Box& sbx = (*flag[lev])[mfi].box();

      Box ubx((*u_g[lev])[mfi].box()); //ubx.shift(0,-1);
      Box vbx((*v_g[lev])[mfi].box()); //vbx.shift(1,-1);
      Box wbx((*w_g[lev])[mfi].box()); //wbx.shift(2,-1);

      solve_v_g_star(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),
          ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
          (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
          (*v_go[lev])[mfi].dataPtr(),     (*p_g[lev])[mfi].dataPtr(),      (*ro_g[lev])[mfi].dataPtr(),
          (*rop_g[lev])[mfi].dataPtr(),    (*rop_go[lev])[mfi].dataPtr(),   (*ep_g[lev])[mfi].dataPtr(),
          (*tau_v_g[lev])[mfi].dataPtr(),  (*d_n[lev])[mfi].dataPtr(),
          (*flux_gE[lev])[mfi].dataPtr(),  (*flux_gN[lev])[mfi].dataPtr(),  (*flux_gT[lev])[mfi].dataPtr(),
          (*mu_g[lev])[mfi].dataPtr(),     (*f_gds[lev])[mfi].dataPtr(),
          (*A_m[lev])[mfi].dataPtr(),      (*b_m[lev])[mfi].dataPtr(),      (*drag_bm[lev])[mfi].dataPtr(),
          (*flag[lev])[mfi].dataPtr(),     &dt, &dx, &dy, &dz);
    }

    fill_mf_bc(lev,*A_m[lev]);
    fill_mf_bc(lev,*b_m[lev]);
    fill_mf_bc(lev,*d_n[lev]);

    eq_id=4;
    mfix_solve_linear_equation(eq_id,lev,(*v_gt[lev]),(*A_m[lev]),(*b_m[lev]));

    // Solve W-Momentum equation
    MultiFab::Copy(*w_gt[lev], *w_g[lev], 0, 0, 1, w_g[lev]->nGrow());
    for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.validbox();
      const Box& sbx = (*flag[lev])[mfi].box();

      Box ubx((*u_g[lev])[mfi].box()); //ubx.shift(0,-1);
      Box vbx((*v_g[lev])[mfi].box()); //vbx.shift(1,-1);
      Box wbx((*w_g[lev])[mfi].box()); //wbx.shift(2,-1);

      solve_w_g_star(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),
          ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
          (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
          (*w_go[lev])[mfi].dataPtr(),     (*p_g[lev])[mfi].dataPtr(),      (*ro_g[lev])[mfi].dataPtr(),
          (*rop_g[lev])[mfi].dataPtr(),    (*rop_go[lev])[mfi].dataPtr(),   (*ep_g[lev])[mfi].dataPtr(),
          (*tau_w_g[lev])[mfi].dataPtr(),  (*d_t[lev])[mfi].dataPtr(),
          (*flux_gE[lev])[mfi].dataPtr(),  (*flux_gN[lev])[mfi].dataPtr(),  (*flux_gT[lev])[mfi].dataPtr(),
          (*mu_g[lev])[mfi].dataPtr(),     (*f_gds[lev])[mfi].dataPtr(),
          (*A_m[lev])[mfi].dataPtr(),      (*b_m[lev])[mfi].dataPtr(),      (*drag_bm[lev])[mfi].dataPtr(),
          (*flag[lev])[mfi].dataPtr(),     &dt, &dx, &dy, &dz);
    }

    fill_mf_bc(lev,*A_m[lev]);
    fill_mf_bc(lev,*b_m[lev]);
    fill_mf_bc(lev,*d_t[lev]);

    eq_id=5;
    mfix_solve_linear_equation(eq_id,lev,(*w_gt[lev]),(*A_m[lev]),(*b_m[lev]));

    int nghost = u_g[lev]->nGrow();
    MultiFab::Copy(*u_g[lev], *u_gt[lev], 0, 0, 1, nghost);
    MultiFab::Copy(*v_g[lev], *v_gt[lev], 0, 0, 1, nghost);
    MultiFab::Copy(*w_g[lev], *w_gt[lev], 0, 0, 1, nghost);

    fill_mf_bc(lev,*u_g[lev]);
    fill_mf_bc(lev,*v_g[lev]);
    fill_mf_bc(lev,*w_g[lev]);
}

void
mfix_level::mfix_solve_for_pp(int lev, Real dt, Real& lnormg, Real& resg)
{
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    // Solve the pressure correction equation
    for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.validbox();
      const Box& sbx = (*flag[lev])[mfi].box();

      solve_pp_g(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),
        (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
        (*p_g[lev])[mfi].dataPtr(),      (*ep_g[lev])[mfi].dataPtr(),
        (*rop_g[lev])[mfi].dataPtr(),    (*rop_go[lev])[mfi].dataPtr(),
        (*ro_g[lev])[mfi].dataPtr(),
        (*rop_gE[lev])[mfi].dataPtr(),   (*rop_gN[lev])[mfi].dataPtr(),   (*rop_gT[lev])[mfi].dataPtr(),
        (*d_e[lev])[mfi].dataPtr(),      (*d_n[lev])[mfi].dataPtr(),      (*d_t[lev])[mfi].dataPtr(),
        (*A_m[lev])[mfi].dataPtr(),      (*b_m[lev])[mfi].dataPtr(),
        (*flag[lev])[mfi].dataPtr(),     &dt,
        &lnormg,                 &resg,     &dx, &dy, &dz);
    }
    pp_g[lev]->setVal(0.);
    fill_mf_bc(lev,*A_m[lev]);
    fill_mf_bc(lev,*b_m[lev]);

    int eq_id=1;
    mfix_solve_linear_equation(eq_id,lev,(*pp_g[lev]),(*A_m[lev]),(*b_m[lev]));

    fill_mf_bc(lev,*pp_g[lev]);
}

void
mfix_level::mfix_correct0(int lev)
{
  for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
  {
     const Box& bx = mfi.validbox();
     const Box& sbx = (*flag[lev])[mfi].box();

     correct0(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),
      (*p_g[lev])[mfi].dataPtr(),      (*pp_g[lev])[mfi].dataPtr(),
      (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
      (*d_e[lev])[mfi].dataPtr(),      (*d_n[lev])[mfi].dataPtr(),      (*d_t[lev])[mfi].dataPtr(),
      (*flag[lev])[mfi].dataPtr());
  }
  fill_mf_bc(lev,*p_g[lev]);
  fill_mf_bc(lev,*u_g[lev]);
  fill_mf_bc(lev,*v_g[lev]);
  fill_mf_bc(lev,*w_g[lev]);
}

void
mfix_level::mfix_physical_prop(int lev, int calc_flag)
{
  for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
  {
     // const Box& bx = mfi.validbox();
     const Box& sbx = (*flag[lev])[mfi].box();

     physical_prop(sbx.loVect(), sbx.hiVect(), &calc_flag,
        (*ro_g[lev])[mfi].dataPtr(), (*p_g[lev])[mfi].dataPtr(),
        (*ep_g[lev])[mfi].dataPtr(), (*rop_g[lev])[mfi].dataPtr(),
        (*flag[lev])[mfi].dataPtr());
  }
  fill_mf_bc(lev,*ro_g[lev]);
  fill_mf_bc(lev,*rop_g[lev]);
}

void
mfix_level::usr3(int lev)
{
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  for (MFIter mfi(*flag[lev]); mfi.isValid(); ++mfi)
  {
     const Box& sbx = (*flag[lev])[mfi].box();

     mfix_usr3( sbx.loVect(),               sbx.hiVect(),
               (*u_g[lev])[mfi].dataPtr(),  (*v_g[lev])[mfi].dataPtr(),
               (*w_g[lev])[mfi].dataPtr(),  (*p_g[lev])[mfi].dataPtr(),
               &dx, &dy, &dz);
   }
}

void
mfix_level::mfix_solve_linear_equation(int eq_id,int lev,MultiFab& sol, MultiFab& matrix, MultiFab& rhs)
{
    int sweep_type, precond_type, max_it;
    Real tol;

    get_solver_params (&eq_id,&sweep_type,&precond_type,&max_it,&tol);

#if(0)
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.validbox();
       const Box& sbx = (*flag[lev])[mfi].box();

       mfix_solve_lin_eq(&eq_id, sol[mfi].dataPtr(), matrix[mfi].dataPtr(),
                         rhs[mfi].dataPtr(),
                         &sweep_type, &tol, &precond_type, &max_it, &ier,
                         sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect());
    }
#else
    solve_bicgstab(sol, rhs, matrix, sweep_type, precond_type, max_it, tol);
#endif
}

void
mfix_level::fill_mf_bc(int lev, MultiFab& mf)
{

  // NOTE -- double check the order of these! (ASA)
  Array<int> bc_type(6);
  bc_type[0] = FOEXTRAP;  // xlo
  bc_type[1] = FOEXTRAP;  // ylo
  bc_type[2] = FOEXTRAP;  // zlo
  bc_type[3] = FOEXTRAP;  // xhi
  bc_type[4] = FOEXTRAP;  // yhi
  bc_type[5] = FOEXTRAP;  // zhi

  // Fill all cell-centered arrays with first-order extrapolation at domain boundaries
  for (MFIter mfi(mf); mfi.isValid(); ++mfi)
  {
      const Box& sbx = mf[mfi].box();
      fill_bc(mf[mfi].dataPtr(),sbx.loVect(),sbx.hiVect(),bc_type.dataPtr());
  }

  // Impose periodic bc's at domain boundaries and fine-fine copies in the interio
  mf.FillBoundary(geom[lev].periodicity());
}

void
mfix_level::fill_mf_bc(int lev, iMultiFab& mf)
{
  mf.FillBoundary(geom[lev].periodicity());
}
