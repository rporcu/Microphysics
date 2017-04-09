#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

mfix_level::~mfix_level ()
{};

mfix_level::mfix_level ()
//mfix_level::mfix_level (const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord)
// AmrCore(rb,max_level_in,n_cell_in,coord)
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
    pc = std::unique_ptr<MFIXParticleContainer> (new MFIXParticleContainer(this));

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

    fluxX.resize(nlevs_max);
    fluxY.resize(nlevs_max);
    fluxZ.resize(nlevs_max);

    ropX.resize(nlevs_max);
    ropY.resize(nlevs_max);
    ropZ.resize(nlevs_max);

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
  // Traditionally, max_step and stop_time do not have prefix.
  {
    ParmParse pp;
    pp.query("max_step", max_step);
    pp.query("stop_time", stop_time);
  }

  // Traditionally, these have prefix "amr", but we will
  // give them prefix mfix to make it clear that they affect the
  // behavior of the solver and not amr (even thought they are read
  // via BoxLib
  {
    ParmParse pp("amr");
    pp.query("check_file", check_file);
    pp.query("check_int", check_int);
    pp.query("plot_file", plot_file);
    pp.query("plot_int", plot_int);
    pp.query("restart_chkfile", restart_chkfile);
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

      MakeNewLevelFromScratch(0, time, ba, dm);

      init_output_vars(&time, &dt);

      // Parse residual strings
      parse_resid_string();

      Real dx = geom[lev].CellSize(0);
      Real dy = geom[lev].CellSize(1);
      Real dz = geom[lev].CellSize(2);

      Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
      Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
      Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

      Box domain(geom[0].Domain());

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

      InitLevelData(lev,dt,time);

      InitIOData ();
    }

    // if max_level > 0, define fine levels

    pc->InitData();
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
    if ( ParallelDescriptor::IOProcessor() )
       std::cout << "BA " << ba << std::endl;
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
    bc_ilo.resize(box_ilo,2);
    bc_ihi.resize(box_ihi,2);
    bc_jlo.resize(box_jlo,2);
    bc_jhi.resize(box_jhi,2);
    bc_klo.resize(box_klo,2);
    bc_khi.resize(box_khi,2);

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
    Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
    Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

    set_domain();

    Box domain(geom[0].Domain());

    // Only call this check on one processor since it has a bunch of print statements
    // if ( ParallelDescriptor::IOProcessor() )
       check_domain(&dx,&dy,&dz,&xlen,&ylen,&zlen,domain.loVect(),domain.hiVect());

    set_bc_area(&dx,&dy,&dz);

    // Convert (mass, volume) flows to velocities.
    set_bc_flow();

    // Only call this check on one processor since it has a bunch of print statements
    if ( ParallelDescriptor::IOProcessor() )
       check_bc_flow();

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

    //
    f_gds[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    f_gds[lev]->setVal(0.);

    //
    drag_bm[lev].reset(new MultiFab(grids[lev],dmap[lev],3,nghost));
    drag_bm[lev]->setVal(0.);

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

    d_n[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    d_n[lev]->setVal(0.);

    tau_v_g[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    tau_v_g[lev]->setVal(0.);

    fluxY[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    fluxY[lev]->setVal(0.);

    ropY[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    ropY[lev]->setVal(0.);

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

    ropZ[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    ropZ[lev]->setVal(0.);

    // ********************************************************************************

    mfix_set_bc_type(lev);
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

          Real residuals[2*8];
          for (int i=0; i<=2*8; ++i)
            residuals[i] = 0.0L;

          // User hooks
          for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
            mfix_usr2();

          // Calculate transport coefficients
          int calc_flag = 1;
          mfix_calc_coeffs(lev,calc_flag);

          // Solve momentum equations
          mfix_solve_for_vels(lev, dt, residuals);

          // Calculate transport coefficients
          mfix_physical_prop(lev,0);

          // Calculate bulk density (epg*ro_g) at cell faces
          mfix_conv_rop(lev,dt);

          // Solve the pressure correction equation
          mfix_solve_for_pp(lev,dt,lnormg,resg, residuals);

          // Apply pressure correction to all Pg, Ug, Vg, Wg
          mfix_correct_0(lev);

          // Update fluid density
          mfix_physical_prop(lev,0);

          // Calculate face mass fluxes
          mfix_calc_mflux(lev);

          // Check for convergence
          ParallelDescriptor::ReduceRealSum(residuals,16);
          converged = check_convergence(&nit, residuals);

          // Display current iteration residuals
          if ( ParallelDescriptor::IOProcessor() )
             display_resid(&nit, residuals);

          // Iterate over cyclic mass flux bc
          if(cyclic_mf==1 && (converged==1 || nit >= max_nit))
            for (MFIter mfi(*fluxX[lev]); mfi.isValid(); ++mfi)
            {
              const Box& sbx = (*ep_g[lev])[mfi].box();

              converged = goal_seek_mflux(sbx.loVect(), sbx.hiVect(), &nit, &gsmf, &delP_MF, &lMFlux,
                (*fluxX[lev])[mfi].dataPtr(),  (*fluxY[lev])[mfi].dataPtr(),  (*fluxZ[lev])[mfi].dataPtr(),
                &dx, &dy, &dz);
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

    Box domain(geom[lev].Domain());

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
    Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
    Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

    const int max_pip = particle_state.size();

    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
       // const Box& bx = mfi.validbox();
       const Box& sbx = (*ep_g[lev])[mfi].box();
       const Box& bx = mfi.validbox();

       Box ubx((*u_g[lev])[mfi].box());
       Box vbx((*v_g[lev])[mfi].box());
       Box wbx((*w_g[lev])[mfi].box());

       mfix_des_time_march(&max_pip,
        sbx.loVect(), sbx.hiVect(),
        ubx.loVect(), ubx.hiVect(),
        vbx.loVect(), vbx.hiVect(),
        wbx.loVect(), wbx.hiVect(),
        bx.loVect(), bx.hiVect(),
        domain.loVect(), domain.hiVect(),
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
        &time, &dt, &dx, &dy, &dz, &xlen, &ylen, &zlen, &nstep);
    }

    fill_mf_bc(lev,*ep_g[lev]);
    fill_mf_bc(lev,*rop_g[lev]);
}

void
mfix_level::output(int lev, int estatus, int finish, int nstep, Real dt, Real time)
{
  const int max_pip = particle_state.size();

  Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
  Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
  Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

  for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
  {
     mfix_output_manager(&max_pip,
      &time, &dt, &xlen, &ylen, &zlen, &nstep,
      particle_state.dataPtr(), des_radius.dataPtr(),
      des_pos_new.dataPtr(),
      des_vel_new.dataPtr(), des_usr_var.dataPtr(),
      omega_new.dataPtr(), &finish);
  }
}

void
mfix_level::InitLevelData(int lev, Real dt, Real time)
{
  Box domain(geom[lev].Domain());

  for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
  {
     const Box& sbx = (*ep_g[lev])[mfi].box();

     Box ubx((*u_g[lev])[mfi].box());
     Box vbx((*v_g[lev])[mfi].box());
     Box wbx((*w_g[lev])[mfi].box());

     zero_norm_vel(sbx.loVect(), sbx.hiVect(),
                   ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
                   (*u_g[lev])[mfi].dataPtr(),     (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
                   bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                   bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect());

     set_bc0(sbx.loVect(), sbx.hiVect(),
             ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
             (*u_g[lev])[mfi].dataPtr(),     (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
             (*p_g[lev])[mfi].dataPtr(),     (*ep_g[lev])[mfi].dataPtr(),
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

  // Allocate the particle arrays
  if (solve_dem)
  {
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi) {
      const int max_pip = particle_state.size();

      mfix_make_arrays_des(&max_pip,
        particle_state.dataPtr(),
        particle_phase.dataPtr(), des_radius.dataPtr(), ro_sol.dataPtr(),
        pvol.dataPtr(), pmass.dataPtr(), omoi.dataPtr(),
        des_pos_new.dataPtr(), des_vel_new.dataPtr(),
        des_usr_var.dataPtr(), omega_new.dataPtr(),
        fc.dataPtr(), tow.dataPtr());
    }

    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi){
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

  // mfix_finl_err_msg();
}

void
mfix_level::mfix_calc_coeffs(int lev, int calc_flag)
{
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
  {
     const Box& sbx = (*ep_g[lev])[mfi].box();
     const Box& bx = mfi.validbox();

     Box ubx((*u_g[lev])[mfi].box());
     Box vbx((*v_g[lev])[mfi].box());
     Box wbx((*w_g[lev])[mfi].box());

     const int max_pip = particle_state.size();

     calc_coeff(sbx.loVect(), sbx.hiVect(),
      ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
      bx.loVect(), bx.hiVect(), &max_pip,   &calc_flag,
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

  for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
  {
     const Box& bx = mfi.validbox();
     const Box& sbx = (*ep_g[lev])[mfi].box();

     Box ubx((*u_g[lev])[mfi].box());
     Box vbx((*v_g[lev])[mfi].box());
     Box wbx((*w_g[lev])[mfi].box());

     calc_coeff_all(sbx.loVect(), sbx.hiVect(),
       ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
       bx.loVect(), bx.hiVect(), &max_pip,
       (*ro_g[lev])[mfi].dataPtr(), (*p_g[lev])[mfi].dataPtr(),
       (*ep_g[lev])[mfi].dataPtr(), (*rop_g[lev])[mfi].dataPtr(),
       (*u_g[lev])[mfi].dataPtr(),  (*v_g[lev])[mfi].dataPtr(),   (*w_g[lev])[mfi].dataPtr(),
       (*mu_g[lev])[mfi].dataPtr(), (*f_gds[lev])[mfi].dataPtr(), (*drag_bm[lev])[mfi].dataPtr(),
       particle_phase.dataPtr(),  particle_state.dataPtr(),
       pvol.dataPtr(), des_pos_new.dataPtr(),
       des_vel_new.dataPtr(), des_radius.dataPtr(), &dx, &dy, &dz);
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

  for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
  {
     const Box& bx = mfi.validbox();
     const Box& sbx = (*ep_g[lev])[mfi].box();

     Box ubx((*u_g[lev])[mfi].box());
     Box vbx((*v_g[lev])[mfi].box());
     Box wbx((*w_g[lev])[mfi].box());

     calc_trd_and_tau(sbx.loVect(), sbx.hiVect(),
       ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
        bx.loVect(),  bx.hiVect(),
       (*tau_u_g[lev])[mfi].dataPtr(),  (*tau_v_g[lev])[mfi].dataPtr(), (*tau_w_g[lev])[mfi].dataPtr(),
       (*trD_g[lev])[mfi].dataPtr(),
       (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),     (*w_g[lev])[mfi].dataPtr(),
       (*lambda_g[lev])[mfi].dataPtr(), (*mu_g[lev])[mfi].dataPtr(),
       &dx, &dy, &dz);
  }

  tau_u_g[lev]->FillBoundary(geom[lev].periodicity());
  tau_v_g[lev]->FillBoundary(geom[lev].periodicity());
  tau_w_g[lev]->FillBoundary(geom[lev].periodicity());

  fill_mf_bc(lev,*trD_g[lev]);
}

void
mfix_level::mfix_init_fluid(int lev)
{
  Box domain(geom[lev].Domain());

  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
  Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
  Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

  for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
  {
     const Box& bx = mfi.validbox();
     const Box& sbx = (*ep_g[lev])[mfi].box();

     Box ubx((*u_g[lev])[mfi].box());
     Box vbx((*v_g[lev])[mfi].box());
     Box wbx((*w_g[lev])[mfi].box());

     init_fluid(sbx.loVect(), sbx.hiVect(),
       ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
        bx.loVect(),  bx.hiVect(), domain.loVect(), domain.hiVect(),
       (*ep_g[lev])[mfi].dataPtr(),     (*ro_g[lev])[mfi].dataPtr(),
       (*rop_g[lev])[mfi].dataPtr(),     (*p_g[lev])[mfi].dataPtr(),
       (*u_g[lev])[mfi].dataPtr(),     (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
       (*mu_g[lev])[mfi].dataPtr(),   (*lambda_g[lev])[mfi].dataPtr(),
       &dx, &dy, &dz, &xlen, &ylen, &zlen );
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
mfix_level::mfix_comp_mean_fields(int lev)
{
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    const int max_pip = particle_state.size();

    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
       const Box& sbx = (*ep_g[lev])[mfi].box();

       comp_mean_fields(sbx.loVect(), sbx.hiVect(),
            &max_pip, (*ep_g[lev])[mfi].dataPtr(),
            particle_state.dataPtr(), des_pos_new.dataPtr(), pvol.dataPtr(),
            &dx, &dy, &dz );
    }
    fill_mf_bc(lev,*ep_g[lev]);
}

void
mfix_level::mfix_calc_mflux(int lev)
{
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  for (MFIter mfi(*u_g[lev]); mfi.isValid(); ++mfi)
  {
    Box ubx((*u_g[lev])[mfi].box());
    Box vbx((*v_g[lev])[mfi].box());
    Box wbx((*w_g[lev])[mfi].box());

     calc_mflux(
       ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
       (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
       (*ropX[lev])[mfi].dataPtr(),   (*ropY[lev])[mfi].dataPtr(),   (*ropZ[lev])[mfi].dataPtr(),
       (*fluxX[lev])[mfi].dataPtr(),  (*fluxY[lev])[mfi].dataPtr(),  (*fluxZ[lev])[mfi].dataPtr(),
       &dx, &dy, &dz);
  }

  // Impose periodic bc's at domain boundaries and fine-fine copies in the interio
  fluxX[lev]->FillBoundary(geom[lev].periodicity());
  fluxY[lev]->FillBoundary(geom[lev].periodicity());
  fluxZ[lev]->FillBoundary(geom[lev].periodicity());
}

void
mfix_level::mfix_conv_rop(int lev, Real dt)
{
    Box domain(geom[lev].Domain());

    for (MFIter mfi(*rop_g[lev]); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.validbox();
       const Box& sbx = (*rop_g[lev])[mfi].box();

       Box ubx((*u_g[lev])[mfi].box());
       Box vbx((*v_g[lev])[mfi].box());
       Box wbx((*w_g[lev])[mfi].box());

       conv_rop(sbx.loVect(), sbx.hiVect(),
         ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
          bx.loVect(),  bx.hiVect(),
         (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
         (*rop_g[lev])[mfi].dataPtr(),
         (*ropX[lev])[mfi].dataPtr(),   (*ropY[lev])[mfi].dataPtr(),   (*ropZ[lev])[mfi].dataPtr());
    }

    ropX[lev]->FillBoundary(geom[lev].periodicity());
    ropY[lev]->FillBoundary(geom[lev].periodicity());
    ropZ[lev]->FillBoundary(geom[lev].periodicity());
}

void
mfix_level::mfix_solve_for_vels(int lev, Real dt, Real (&residuals)[16])
{
    Box domain(geom[lev].Domain());

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    // Matrix and rhs vector
    BoxArray x_edge_ba = grids[lev];
    x_edge_ba.surroundingNodes(0);
    A_m[lev].reset(new MultiFab(x_edge_ba,dmap[lev],7,0));
    b_m[lev].reset(new MultiFab(x_edge_ba,dmap[lev],1,0));

    // Solve U-Momentum equation
    MultiFab::Copy(*u_gt[lev], *u_g[lev], 0, 0, 1, u_g[lev]->nGrow());
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.validbox();
      const Box& sbx = (*ep_g[lev])[mfi].box();
      Box abx((*A_m[lev])[mfi].box());

      Box ubx((*u_g[lev])[mfi].box());
      Box vbx((*v_g[lev])[mfi].box());
      Box wbx((*w_g[lev])[mfi].box());

      solve_u_g_star(sbx.loVect(), sbx.hiVect(),
          ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(),
          wbx.loVect(), wbx.hiVect(), abx.loVect(), abx.hiVect(),
           bx.loVect(),  bx.hiVect(),
          (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
          (*u_go[lev])[mfi].dataPtr(),     (*p_g[lev])[mfi].dataPtr(),      (*ro_g[lev])[mfi].dataPtr(),
          (*rop_g[lev])[mfi].dataPtr(),    (*rop_go[lev])[mfi].dataPtr(),   (*ep_g[lev])[mfi].dataPtr(),
          (*tau_u_g[lev])[mfi].dataPtr(),  (*d_e[lev])[mfi].dataPtr(),
          (*fluxX[lev])[mfi].dataPtr(),  (*fluxY[lev])[mfi].dataPtr(),  (*fluxZ[lev])[mfi].dataPtr(),
          (*mu_g[lev])[mfi].dataPtr(),     (*f_gds[lev])[mfi].dataPtr(),
          (*A_m[lev])[mfi].dataPtr(),      (*b_m[lev])[mfi].dataPtr(),      (*drag_bm[lev])[mfi].dataPtr(),
          bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
          bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect(),
          &dt, &dx, &dy, &dz, residuals);
    }


    int eq_id=3;
    mfix_solve_linear_equation(eq_id,lev,(*u_gt[lev]),(*A_m[lev]),(*b_m[lev]));

    // Solve V-Momentum equation

    // Matrix and rhs vector
    BoxArray y_edge_ba = grids[lev];
    y_edge_ba.surroundingNodes(1);
    A_m[lev].reset(new MultiFab(y_edge_ba,dmap[lev],7,0));
    b_m[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,0));

    MultiFab::Copy(*v_gt[lev], *v_g[lev], 0, 0, 1, v_g[lev]->nGrow());
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.validbox();
      const Box& sbx = (*ep_g[lev])[mfi].box();
      Box abx((*A_m[lev])[mfi].box());

      Box ubx((*u_g[lev])[mfi].box());
      Box vbx((*v_g[lev])[mfi].box());
      Box wbx((*w_g[lev])[mfi].box());

      solve_v_g_star(sbx.loVect(), sbx.hiVect(),
          ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(),
          wbx.loVect(), wbx.hiVect(), abx.loVect(), abx.hiVect(),
           bx.loVect(),  bx.hiVect(),
          (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
          (*v_go[lev])[mfi].dataPtr(),     (*p_g[lev])[mfi].dataPtr(),      (*ro_g[lev])[mfi].dataPtr(),
          (*rop_g[lev])[mfi].dataPtr(),    (*rop_go[lev])[mfi].dataPtr(),   (*ep_g[lev])[mfi].dataPtr(),
          (*tau_v_g[lev])[mfi].dataPtr(),  (*d_n[lev])[mfi].dataPtr(),
          (*fluxX[lev])[mfi].dataPtr(),  (*fluxY[lev])[mfi].dataPtr(),  (*fluxZ[lev])[mfi].dataPtr(),
          (*mu_g[lev])[mfi].dataPtr(),     (*f_gds[lev])[mfi].dataPtr(),
          (*A_m[lev])[mfi].dataPtr(),      (*b_m[lev])[mfi].dataPtr(),      (*drag_bm[lev])[mfi].dataPtr(),
          bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
          bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect(),
          &dt, &dx, &dy, &dz, residuals);
    }

    eq_id=4;
    mfix_solve_linear_equation(eq_id,lev,(*v_gt[lev]),(*A_m[lev]),(*b_m[lev]));

    // Solve W-Momentum equation

    // Matrix and rhs vector
    BoxArray z_edge_ba = grids[lev];
    z_edge_ba.surroundingNodes(2);
    A_m[lev].reset(new MultiFab(z_edge_ba,dmap[lev],7,0));
    b_m[lev].reset(new MultiFab(z_edge_ba,dmap[lev],1,0));

    MultiFab::Copy(*w_gt[lev], *w_g[lev], 0, 0, 1, w_g[lev]->nGrow());
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.validbox();
      const Box& sbx = (*ep_g[lev])[mfi].box();
      Box abx((*A_m[lev])[mfi].box());

      Box ubx((*u_g[lev])[mfi].box());
      Box vbx((*v_g[lev])[mfi].box());
      Box wbx((*w_g[lev])[mfi].box());

      solve_w_g_star(sbx.loVect(), sbx.hiVect(),
          ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(),
          wbx.loVect(), wbx.hiVect(), abx.loVect(), abx.hiVect(),
           bx.loVect(),  bx.hiVect(),
          (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
          (*w_go[lev])[mfi].dataPtr(),     (*p_g[lev])[mfi].dataPtr(),      (*ro_g[lev])[mfi].dataPtr(),
          (*rop_g[lev])[mfi].dataPtr(),    (*rop_go[lev])[mfi].dataPtr(),   (*ep_g[lev])[mfi].dataPtr(),
          (*tau_w_g[lev])[mfi].dataPtr(),  (*d_t[lev])[mfi].dataPtr(),
          (*fluxX[lev])[mfi].dataPtr(),  (*fluxY[lev])[mfi].dataPtr(),  (*fluxZ[lev])[mfi].dataPtr(),
          (*mu_g[lev])[mfi].dataPtr(),     (*f_gds[lev])[mfi].dataPtr(),
          (*A_m[lev])[mfi].dataPtr(),      (*b_m[lev])[mfi].dataPtr(),      (*drag_bm[lev])[mfi].dataPtr(),
          bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
          bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect(),
          &dt, &dx, &dy, &dz, residuals);
    }

    eq_id=5;
    mfix_solve_linear_equation(eq_id,lev,(*w_gt[lev]),(*A_m[lev]),(*b_m[lev]));

    MultiFab::Copy(*u_g[lev], *u_gt[lev], 0, 0, 1, u_g[lev]->nGrow());
    MultiFab::Copy(*v_g[lev], *v_gt[lev], 0, 0, 1, v_g[lev]->nGrow());
    MultiFab::Copy(*w_g[lev], *w_gt[lev], 0, 0, 1, w_g[lev]->nGrow());

    u_g[lev]->FillBoundary(geom[lev].periodicity());
    v_g[lev]->FillBoundary(geom[lev].periodicity());
    w_g[lev]->FillBoundary(geom[lev].periodicity());
}

void
mfix_level::mfix_solve_for_pp(int lev, Real dt, Real& lnormg, Real& resg, Real (&residuals)[16])
{
    Box domain(geom[lev].Domain());

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    // Matrix and rhs vector
    A_m[lev].reset(new MultiFab(grids[lev],dmap[lev],7,0));
    b_m[lev].reset(new MultiFab(grids[lev],dmap[lev],1,0));

    // Solve the pressure correction equation
    MultiFab b_mmax(b_m[lev]->boxArray(),dmap[lev],1,b_m[lev]->nGrow());
    b_mmax.setVal(0.);

    // Solve the pressure correction equation
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.validbox();
      const Box& sbx = (*ep_g[lev])[mfi].box();
      Box abx((*A_m[lev])[mfi].box());

      Box ubx((*u_g[lev])[mfi].box());
      Box vbx((*v_g[lev])[mfi].box());
      Box wbx((*w_g[lev])[mfi].box());

      solve_pp_g(sbx.loVect(), sbx.hiVect(),
        ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
        abx.loVect(), abx.hiVect(), bx.loVect(),  bx.hiVect(),
        (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
        (*p_g[lev])[mfi].dataPtr(),      (*ep_g[lev])[mfi].dataPtr(),
        (*rop_g[lev])[mfi].dataPtr(),    (*rop_go[lev])[mfi].dataPtr(),
        (*ro_g[lev])[mfi].dataPtr(),
        (*ropX[lev])[mfi].dataPtr(),   (*ropY[lev])[mfi].dataPtr(),   (*ropZ[lev])[mfi].dataPtr(),
        (*d_e[lev])[mfi].dataPtr(),      (*d_n[lev])[mfi].dataPtr(),      (*d_t[lev])[mfi].dataPtr(),
        (*A_m[lev])[mfi].dataPtr(),      (*b_m[lev])[mfi].dataPtr(),           b_mmax[mfi].dataPtr(),
        &dt, &dx, &dy, &dz, domain.loVect(), domain.hiVect(), residuals);
    }
    pp_g[lev]->setVal(0.);

    int eq_id=1;
    mfix_solve_linear_equation(eq_id,lev,(*pp_g[lev]),(*A_m[lev]),(*b_m[lev]));
    fill_mf_bc(lev,*pp_g[lev]);
}

void
mfix_level::mfix_correct_0(int lev)
{
  Box domain(geom[lev].Domain());
  for (MFIter mfi(*p_g[lev]); mfi.isValid(); ++mfi)
  {
     const Box& bx = mfi.validbox();
     const Box& sbx = (*p_g[lev])[mfi].box();

     Box ubx((*u_g[lev])[mfi].box());
     Box vbx((*v_g[lev])[mfi].box());
     Box wbx((*w_g[lev])[mfi].box());

     correct_0(sbx.loVect(), sbx.hiVect(),
               ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
                bx.loVect(),  bx.hiVect(), domain.loVect(), domain.hiVect(),
      (*p_g[lev])[mfi].dataPtr(),      (*pp_g[lev])[mfi].dataPtr(),
      (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
      (*d_e[lev])[mfi].dataPtr(),      (*d_n[lev])[mfi].dataPtr(),      (*d_t[lev])[mfi].dataPtr());
  }

  fill_mf_bc(lev,*p_g[lev]);

  u_g[lev]->FillBoundary(geom[lev].periodicity());
  v_g[lev]->FillBoundary(geom[lev].periodicity());
  w_g[lev]->FillBoundary(geom[lev].periodicity());
}

void
mfix_level::mfix_physical_prop(int lev, int calc_flag)
{
  for (MFIter mfi(*p_g[lev]); mfi.isValid(); ++mfi)
  {
     const Box& bx = mfi.validbox();
     const Box& sbx = (*p_g[lev])[mfi].box();

     physical_prop(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),&calc_flag,
        (*ro_g[lev])[mfi].dataPtr(), (*p_g[lev])[mfi].dataPtr(),
        (*ep_g[lev])[mfi].dataPtr(), (*rop_g[lev])[mfi].dataPtr());
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

  for (MFIter mfi(*p_g[lev]); mfi.isValid(); ++mfi)
  {
     const Box& sbx = (*p_g[lev])[mfi].box();
     Box ubx((*u_g[lev])[mfi].box());
     Box vbx((*v_g[lev])[mfi].box());
     Box wbx((*w_g[lev])[mfi].box());

     mfix_usr3((*u_g[lev])[mfi].dataPtr(), ubx.loVect(), ubx.hiVect(),
               (*v_g[lev])[mfi].dataPtr(), vbx.loVect(), vbx.hiVect(),
               (*w_g[lev])[mfi].dataPtr(), wbx.loVect(), wbx.hiVect(),
               (*p_g[lev])[mfi].dataPtr(), sbx.loVect(), sbx.hiVect(),
               &dx, &dy, &dz);
   }
}

void
mfix_level::mfix_solve_linear_equation(int eq_id,int lev,MultiFab& sol, MultiFab& matrix, MultiFab& rhs)
{
    int sweep_type, precond_type, max_it;
    Real tol;

    get_solver_params (&eq_id,&sweep_type,&precond_type,&max_it,&tol);

    solve_bicgstab(sol, rhs, matrix, sweep_type, precond_type, max_it, tol, lev);
}

void
mfix_level::mfix_set_bc_type(int lev)
{
  Box domain(geom[lev].Domain());
  for (MFIter mfi((*ep_g[lev])); mfi.isValid(); ++mfi)
  {
      const Box& sbx = (*ep_g[lev])[mfi].box();
      set_bc_type(sbx.loVect(),sbx.hiVect(), bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                  bc_jlo.dataPtr(), bc_jhi.dataPtr(), bc_klo.dataPtr(), bc_khi.dataPtr(),
                  domain.loVect(),domain.hiVect());
  }

}

void
mfix_level::fill_mf_bc(int lev, MultiFab& mf)
{
  Box domain(geom[lev].Domain());

  // Impose periodic bc's at domain boundaries and fine-fine copies in the interio
  mf.FillBoundary(geom[lev].periodicity());
  // Fill all cell-centered arrays with first-order extrapolation at domain boundaries
  for (MFIter mfi(mf); mfi.isValid(); ++mfi)
  {
      const Box& sbx = mf[mfi].box();
      fill_bc0(mf[mfi].dataPtr(),sbx.loVect(),sbx.hiVect(),
               bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
               bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect());
  }
}
