#include <fstream>
#include <iomanip>

#include <ParmParse.H>
#include <Geometry.H>
#include <VisMF.H>
#include <iMultiFab.H>

#include <mfix_level.H>
#include <mfix_F.H>

int main (int argc, char* argv[])
{

  // Copy arguments into MFIX
  for(int i=1; i < argc; i++) {
    int nlen = strlen(argv[i]);
    mfix_add_argument(argv[i], &nlen);
  }

  BoxLib::Initialize(argc,argv);

  Real strt_time = ParallelDescriptor::second();

  // Size of the entire domain
  int imax,jmax,kmax;
  int max_nit;
  int solve_fluid;
  int solve_dem;
  int steady_state;
  Real dt, dt_min, dt_max, tstop, time;
  int nstep=0;  // Number of time steps
  Real normg;
  int set_normg;
  int cyclic_x, cyclic_y, cyclic_z, cyclic_mf;

  mfix_get_data(
    &imax, &jmax, &kmax,
    &solve_fluid, 
    &solve_dem,
    &steady_state,
    &dt, &dt_min, &dt_max, &tstop, &time, &max_nit,
    &normg, &set_normg,
    &cyclic_x, &cyclic_y, &cyclic_z, &cyclic_mf);

  IntVect dom_lo(IntVect(D_DECL(0,0,0)));
  IntVect dom_hi(IntVect(D_DECL(imax-1, jmax-1, kmax-1)));

  Box domain(dom_lo, dom_hi);

  // Initialize the boxarray "ba" from the single box "bx"
  BoxArray ba;
  ba.define(domain);

  int max_grid_size = 1024;

  // Break up boxarray "ba" into chunks no larger than
  // "max_grid_size" along a direction.
  ba.maxSize(max_grid_size);

  // This defines the physical size of the box.
  // Right now the box is [0,1] in each direction.
  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    real_box.setLo(n, 0.0);
    real_box.setHi(n, 1.0);
  }

  // This says we are using Cartesian coordinates
  int coord = 0;

  // This sets the boundary conditions to be doubly or triply periodic
  int is_periodic[BL_SPACEDIM];
  for (int i = 0; i < BL_SPACEDIM; i++) {
    is_periodic[i] = 0;
  }

#if 0
  // This defines a Geometry object
  Geometry geom;
  geom.define(domain,&real_box,coord,is_periodic);
#endif

  int max_level = 0;
  int lev = 0;
  Array<int> n_cell(3);
  n_cell[0] = imax;
  n_cell[1] = jmax;
  n_cell[2] = kmax;

  // Note that the constructor constructs the Geometry object now.
  mfix_level my_mfix(&real_box,max_level,n_cell,coord);

  my_mfix.Init();

  // define dx[]
  // const Real* dx = mfix_level.geom.CellSize();

#if 0
  Array<int> particle_state (  nparticles);
  Array<int> particle_phase (  nparticles);

  Array<Real> des_radius    (  nparticles);
  Array<Real> ro_sol        (  nparticles);
  Array<Real> pvol          (  nparticles);
  Array<Real> pmass         (  nparticles);
  Array<Real> omoi          (  nparticles);
  Array<Real> des_pos_new   (3*nparticles);
  Array<Real> des_vel_new   (3*nparticles);
  Array<Real> des_usr_var   (  nparticles);
  Array<Real> omega_new     (3*nparticles);
  Array<Real> des_acc_old   (3*nparticles);
  Array<Real> rot_acc_old   (3*nparticles);
  Array<Real> drag_fc       (3*nparticles);
  Array<Real> fc            (3*nparticles);
  Array<Real> tow           (3*nparticles);
  Array<int> pairs          (12*nparticles);
#endif

  int pair_count=0;

  my_mfix.call_main(0);
#if 0
  for (MFIter mfi(flag); mfi.isValid(); ++mfi)
     mfix_MAIN(
               &time, &dt, &nstep,
               u_g[mfi].dataPtr(),     v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
               u_go[mfi].dataPtr(),    v_go[mfi].dataPtr(),     w_go[mfi].dataPtr(),
               p_g[mfi].dataPtr(),     p_go[mfi].dataPtr(),     pp_g[mfi].dataPtr(),
               ep_g[mfi].dataPtr(),    ep_go[mfi].dataPtr(),
               ro_g[mfi].dataPtr(),    ro_go[mfi].dataPtr(),
               rop_g[mfi].dataPtr(),   rop_go[mfi].dataPtr(),
               rop_gE[mfi].dataPtr(),  rop_gN[mfi].dataPtr(),   rop_gT[mfi].dataPtr(),
               d_e[mfi].dataPtr(),     d_n[mfi].dataPtr(),      d_t[mfi].dataPtr(),
               tau_u_g[mfi].dataPtr(), tau_v_g[mfi].dataPtr(),  tau_w_g[mfi].dataPtr(),
               flux_gE[mfi].dataPtr(), flux_gN[mfi].dataPtr(),  flux_gT[mfi].dataPtr(),
               trD_g[mfi].dataPtr(),   lambda_g[mfi].dataPtr(), mu_g[mfi].dataPtr(),
               f_gds[mfi].dataPtr(),   A_m[mfi].dataPtr(),      b_m[mfi].dataPtr(),
               drag_bm[mfi].dataPtr(),
               flag[mfi].dataPtr(),
               particle_state.dataPtr(),
               particle_phase.dataPtr(), des_radius.dataPtr(), ro_sol.dataPtr(),
               pvol.dataPtr(), pmass.dataPtr(), omoi.dataPtr(),
               des_pos_new.dataPtr(), des_vel_new.dataPtr(),
               des_usr_var.dataPtr(), omega_new.dataPtr(), des_acc_old.dataPtr(),
               rot_acc_old.dataPtr(), drag_fc.dataPtr(), fc.dataPtr(), tow.dataPtr());
#endif

  int finish, estatus;
  finish = 0;
  estatus = 0;

  Real prev_dt; // Actual dt used to solve fluid

  my_mfix.output(lev,estatus,finish);

#if 0
  // Call to output before entering time march loop
  for (MFIter mfi(flag); mfi.isValid(); ++mfi)
    mfix_output_manager(
      &time, &dt, &nstep,
      ep_g[mfi].dataPtr(),    p_g[mfi].dataPtr(),
      ro_g[mfi].dataPtr(),   rop_g[mfi].dataPtr(),
      u_g[mfi].dataPtr(),    v_g[mfi].dataPtr(),
      w_g[mfi].dataPtr(),
      particle_state.dataPtr(), des_radius.dataPtr(),
      ro_sol.dataPtr(), des_pos_new.dataPtr(),
      des_vel_new.dataPtr(), des_usr_var.dataPtr(),
      omega_new.dataPtr(), &estatus, &finish);
#endif

  my_mfix.evolve(0,estatus,finish);
#if 0
  do {
    for (MFIter mfi(flag); mfi.isValid(); ++mfi)
      mfix_usr1();

    if(solve_fluid) {

      // Update boundary conditions
      for (MFIter mfi(flag); mfi.isValid(); ++mfi)
        set_bc1(
          &time,                   &dt,
          p_g[mfi].dataPtr(),      ep_g[mfi].dataPtr(),
          ro_g[mfi].dataPtr(),     rop_g[mfi].dataPtr(),
          u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
          flux_gE[mfi].dataPtr(),  flux_gN[mfi].dataPtr(),  flux_gT[mfi].dataPtr(),
          flag[mfi].dataPtr());

      // Calculate transport coefficients
      for (MFIter mfi(flag); mfi.isValid(); ++mfi)
        calc_coeff_all(
          ro_g[mfi].dataPtr(), p_g[mfi].dataPtr(),
          ep_g[mfi].dataPtr(), rop_g[mfi].dataPtr(),
          u_g[mfi].dataPtr(),  v_g[mfi].dataPtr(),   w_g[mfi].dataPtr(),
          mu_g[mfi].dataPtr(), f_gds[mfi].dataPtr(), drag_bm[mfi].dataPtr(),
          particle_phase.dataPtr(),  particle_state.dataPtr(),
          pvol.dataPtr(), des_pos_new.dataPtr(),
          des_vel_new.dataPtr(), des_radius.dataPtr(),
          flag[mfi].dataPtr());

      // Calculate the stress tensor trace and cross terms for all phases.
      for (MFIter mfi(flag); mfi.isValid(); ++mfi)
        calc_trd_and_tau(
          tau_u_g[mfi].dataPtr(),  tau_v_g[mfi].dataPtr(), tau_w_g[mfi].dataPtr(),
          trD_g[mfi].dataPtr(),    ep_g[mfi].dataPtr(),
          u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),     w_g[mfi].dataPtr(),
          lambda_g[mfi].dataPtr(), mu_g[mfi].dataPtr(),
          flag[mfi].dataPtr());

      // Backup field variable to old
      MultiFab::Copy(ep_go,  ep_g,  0, 0, 1, nghost);
      MultiFab::Copy(p_go,   p_g,   0, 0, 1, nghost);
      MultiFab::Copy(ro_go,  ro_g,  0, 0, 1, nghost);
      MultiFab::Copy(rop_go, rop_g, 0, 0, 1, nghost);
      MultiFab::Copy(u_go,   u_g,   0, 0, 1, nghost);
      MultiFab::Copy(v_go,   v_g,   0, 0, 1, nghost);
      MultiFab::Copy(w_go,   w_g,   0, 0, 1, nghost);

      // Loop over iterate for auto time-step size adjustment
      int reiterate;
      do {
        prev_dt = dt;

        // Calculate bulk density (epg*ro_g) at cell faces
        for (MFIter mfi(flag); mfi.isValid(); ++mfi){
          const Box& bx=mfi.validbox();
          conv_rop(bx.loVect(), bx.hiVect(),
            u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
            rop_g[mfi].dataPtr(),
            rop_gE[mfi].dataPtr(),   rop_gN[mfi].dataPtr(),   rop_gT[mfi].dataPtr(),
            flag[mfi].dataPtr(),     &dt);
        }

        // Calculate face mass fluxes
        for (MFIter mfi(flag); mfi.isValid(); ++mfi)
          calc_mflux(
            u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
            rop_gE[mfi].dataPtr(),   rop_gN[mfi].dataPtr(),   rop_gT[mfi].dataPtr(),
            flux_gE[mfi].dataPtr(),  flux_gN[mfi].dataPtr(),  flux_gT[mfi].dataPtr(),
            flag[mfi].dataPtr());

        // Update boundary conditions
        for (MFIter mfi(flag); mfi.isValid(); ++mfi)
          set_bc1(
            &time,                   &dt,
            p_g[mfi].dataPtr(),      ep_g[mfi].dataPtr(),
            ro_g[mfi].dataPtr(),     rop_g[mfi].dataPtr(),
            u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
            flux_gE[mfi].dataPtr(),  flux_gN[mfi].dataPtr(),  flux_gT[mfi].dataPtr(),
            flag[mfi].dataPtr());


        int converged=0;
        int nit=0;          // number of iterations
        int gsmf=0;         // number of outer iterations for goal seek mass flux (GSMF)
        Real delP_MF=0.0L;  // actual GSMF pressure drop
        Real lMFlux=0.0L;   // actual GSMF mass flux
        Real resg=0.0L;     // fluid pressure residual

        int lset_normg=1-set_normg;
        Real lnormg=normg;

        ///////////////// ---- call to iterate -------- /////////////////
        do {
          nit++;

          // User hooks
          for (MFIter mfi(flag); mfi.isValid(); ++mfi)
            mfix_usr2();

          // Calculate transport coefficients
          int level=1;
          for (MFIter mfi(flag); mfi.isValid(); ++mfi)
            calc_coeff(
              flag[mfi].dataPtr(), &level,
              ro_g[mfi].dataPtr(), p_g[mfi].dataPtr(),
              ep_g[mfi].dataPtr(), rop_g[mfi].dataPtr(),
              u_g[mfi].dataPtr(),  v_g[mfi].dataPtr(),   w_g[mfi].dataPtr(),
              mu_g[mfi].dataPtr(), f_gds[mfi].dataPtr(), drag_bm[mfi].dataPtr(),
              particle_phase.dataPtr(),  particle_state.dataPtr(),
              pvol.dataPtr(), des_pos_new.dataPtr(),
              des_vel_new.dataPtr(), des_radius.dataPtr());

          // Solve U-Momentum equation
          {
            MultiFab::Copy(u_gt, u_g, 0, 0, 1, nghost);
            for (MFIter mfi(flag); mfi.isValid(); ++mfi)
              solve_u_g_star(
                u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
                u_go[mfi].dataPtr(),     p_g[mfi].dataPtr(),      ro_g[mfi].dataPtr(),
                rop_g[mfi].dataPtr(),    rop_go[mfi].dataPtr(),   ep_g[mfi].dataPtr(),
                tau_u_g[mfi].dataPtr(),  d_e[mfi].dataPtr(),
                flux_gE[mfi].dataPtr(),  flux_gN[mfi].dataPtr(),  flux_gT[mfi].dataPtr(),
                mu_g[mfi].dataPtr(),     f_gds[mfi].dataPtr(),
                A_m[mfi].dataPtr(),      b_m[mfi].dataPtr(),      drag_bm[mfi].dataPtr(),
                flag[mfi].dataPtr(),     &dt);

            int eq_id=3;
            for (MFIter mfi(flag); mfi.isValid(); ++mfi)
              mfix_solve_lin_eq(&eq_id, u_gt[mfi].dataPtr(),
                A_m[mfi].dataPtr(),      b_m[mfi].dataPtr());
          }

          // Solve V-Momentum equation
          {
            MultiFab::Copy(v_gt, v_g, 0, 0, 1, nghost);
            for (MFIter mfi(flag); mfi.isValid(); ++mfi)
              solve_v_g_star(
                u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
                v_go[mfi].dataPtr(),     p_g[mfi].dataPtr(),      ro_g[mfi].dataPtr(),
                rop_g[mfi].dataPtr(),    rop_go[mfi].dataPtr(),   ep_g[mfi].dataPtr(),
                tau_v_g[mfi].dataPtr(),  d_n[mfi].dataPtr(),
                flux_gE[mfi].dataPtr(),  flux_gN[mfi].dataPtr(),  flux_gT[mfi].dataPtr(),
                mu_g[mfi].dataPtr(),     f_gds[mfi].dataPtr(),
                A_m[mfi].dataPtr(),      b_m[mfi].dataPtr(),      drag_bm[mfi].dataPtr(),
                flag[mfi].dataPtr(),     &dt);

            int eq_id=4;
            for (MFIter mfi(flag); mfi.isValid(); ++mfi)
              mfix_solve_lin_eq(&eq_id, v_gt[mfi].dataPtr(),
                A_m[mfi].dataPtr(),      b_m[mfi].dataPtr());
          }

          // Solve W-Momentum equation
          {
            MultiFab::Copy(w_gt, w_g, 0, 0, 1, nghost);
            for (MFIter mfi(flag); mfi.isValid(); ++mfi)
              solve_w_g_star(
                u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
                w_go[mfi].dataPtr(),     p_g[mfi].dataPtr(),      ro_g[mfi].dataPtr(),
                rop_g[mfi].dataPtr(),    rop_go[mfi].dataPtr(),   ep_g[mfi].dataPtr(),
                tau_w_g[mfi].dataPtr(),  d_t[mfi].dataPtr(),
                flux_gE[mfi].dataPtr(),  flux_gN[mfi].dataPtr(),  flux_gT[mfi].dataPtr(),
                mu_g[mfi].dataPtr(),     f_gds[mfi].dataPtr(),
                A_m[mfi].dataPtr(),      b_m[mfi].dataPtr(),      drag_bm[mfi].dataPtr(),
                flag[mfi].dataPtr(),     &dt);

            int eq_id=5;
            for (MFIter mfi(flag); mfi.isValid(); ++mfi)
              mfix_solve_lin_eq(&eq_id, w_gt[mfi].dataPtr(),
                A_m[mfi].dataPtr(),      b_m[mfi].dataPtr());
          }

          MultiFab::Copy(u_g, u_gt, 0, 0, 1, nghost);
          MultiFab::Copy(v_g, v_gt, 0, 0, 1, nghost);
          MultiFab::Copy(w_g, w_gt, 0, 0, 1, nghost);

          // Calculate transport coefficients
          level=0;
          for (MFIter mfi(flag); mfi.isValid(); ++mfi)
            physical_prop(&level,
              ro_g[mfi].dataPtr(), p_g[mfi].dataPtr(),
              ep_g[mfi].dataPtr(), rop_g[mfi].dataPtr(),
              flag[mfi].dataPtr());

          // Calculate bulk density (epg*ro_g) at cell faces
          for (MFIter mfi(flag); mfi.isValid(); ++mfi){
            const Box& bx=mfi.validbox();
            conv_rop(bx.loVect(), bx.hiVect(),
              u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
              rop_g[mfi].dataPtr(),
              rop_gE[mfi].dataPtr(),   rop_gN[mfi].dataPtr(),   rop_gT[mfi].dataPtr(),
              flag[mfi].dataPtr(),     &dt);
          }

          // Solve the pressure correction equation
          for (MFIter mfi(flag); mfi.isValid(); ++mfi)
            solve_pp_g(
              u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
              p_g[mfi].dataPtr(),      ep_g[mfi].dataPtr(),
              rop_g[mfi].dataPtr(),    rop_go[mfi].dataPtr(),
              ro_g[mfi].dataPtr(),     pp_g[mfi].dataPtr(),
              rop_gE[mfi].dataPtr(),   rop_gN[mfi].dataPtr(),   rop_gT[mfi].dataPtr(),
              d_e[mfi].dataPtr(),      d_n[mfi].dataPtr(),      d_t[mfi].dataPtr(),
              A_m[mfi].dataPtr(),      b_m[mfi].dataPtr(),
              flag[mfi].dataPtr(),     &dt,
              &lnormg,                 &resg);

            int eq_id=1;
            for (MFIter mfi(flag); mfi.isValid(); ++mfi)
              mfix_solve_lin_eq(&eq_id,  pp_g[mfi].dataPtr(),
                A_m[mfi].dataPtr(),      b_m[mfi].dataPtr());

          // Correct fluid pressure and velocities
          for (MFIter mfi(flag); mfi.isValid(); ++mfi)
            correct0(
              p_g[mfi].dataPtr(),      pp_g[mfi].dataPtr(),
              u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
              d_e[mfi].dataPtr(),      d_n[mfi].dataPtr(),      d_t[mfi].dataPtr(),
              flag[mfi].dataPtr());

          // Update fluid density
          level=0;
          for (MFIter mfi(flag); mfi.isValid(); ++mfi)
            physical_prop(&level,
              ro_g[mfi].dataPtr(), p_g[mfi].dataPtr(),
              ep_g[mfi].dataPtr(), rop_g[mfi].dataPtr(),
              flag[mfi].dataPtr());

          // Update wall velocities
          for (MFIter mfi(flag); mfi.isValid(); ++mfi)
            set_wall_bc(
              u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
              flag[mfi].dataPtr());

          // Calculate face mass fluxes
          for (MFIter mfi(flag); mfi.isValid(); ++mfi)
            calc_mflux(
              u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
              rop_gE[mfi].dataPtr(),   rop_gN[mfi].dataPtr(),   rop_gT[mfi].dataPtr(),
              flux_gE[mfi].dataPtr(),  flux_gN[mfi].dataPtr(),  flux_gT[mfi].dataPtr(),
              flag[mfi].dataPtr());

          // Update boundary conditions
          for (MFIter mfi(flag); mfi.isValid(); ++mfi)
            set_bc1(
              &time,                   &dt,
              p_g[mfi].dataPtr(),      ep_g[mfi].dataPtr(),
              ro_g[mfi].dataPtr(),     rop_g[mfi].dataPtr(),
              u_g[mfi].dataPtr(),      v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
              flux_gE[mfi].dataPtr(),  flux_gN[mfi].dataPtr(),  flux_gT[mfi].dataPtr(),
              flag[mfi].dataPtr());

          // Display current iteration residuals
          display_resid(&nit);

          // Check for convergence
          converged = check_convergence(&nit);

          // Iterate over cyclic mass flux bc
          if(cyclic_mf==1 && (converged==1 || nit >= max_nit))
            for (MFIter mfi(flag); mfi.isValid(); ++mfi)
              converged = goal_seek_mFlux(&nit, &gsmf, &delP_MF, &lMFlux,
                flux_gE[mfi].dataPtr(),  flux_gN[mfi].dataPtr(),  flux_gT[mfi].dataPtr(),
                flag[mfi].dataPtr());

        } while(converged==0 && nit<max_nit);

        // Adjust time step if iteration failed.
        reiterate = mfix_adjustdt(&converged, &nit, &dt);
        if(reiterate == 1) {

          // Reset the field variables
          MultiFab::Copy(ep_g,  ep_go,  0, 0, 1, nghost);
          MultiFab::Copy(p_g,   p_go,   0, 0, 1, nghost);
          MultiFab::Copy(ro_g,  ro_go,  0, 0, 1, nghost);
          MultiFab::Copy(rop_g, rop_go, 0, 0, 1, nghost);
          MultiFab::Copy(u_g,   u_go,   0, 0, 1, nghost);
          MultiFab::Copy(v_g,   v_go,   0, 0, 1, nghost);
          MultiFab::Copy(w_g,   w_go,   0, 0, 1, nghost);

          // Recalculate all coefficients (JM: not sure why)
          for (MFIter mfi(flag); mfi.isValid(); ++mfi)
            calc_coeff_all(
              ro_g[mfi].dataPtr(), p_g[mfi].dataPtr(),
              ep_g[mfi].dataPtr(), rop_g[mfi].dataPtr(),
              u_g[mfi].dataPtr(),  v_g[mfi].dataPtr(),   w_g[mfi].dataPtr(),
              mu_g[mfi].dataPtr(), f_gds[mfi].dataPtr(), drag_bm[mfi].dataPtr(),
              particle_phase.dataPtr(),  particle_state.dataPtr(),
              pvol.dataPtr(), des_pos_new.dataPtr(),
              des_vel_new.dataPtr(), des_radius.dataPtr(),
              flag[mfi].dataPtr());
        }
      }while (reiterate==1);
    }

    if(solve_dem) {
      for (MFIter mfi(flag); mfi.isValid(); ++mfi)
        mfix_des_time_march(
          ep_g[mfi].dataPtr(),      p_g[mfi].dataPtr(),
          u_g[mfi].dataPtr(),       v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
          ro_g[mfi].dataPtr(),      rop_g[mfi].dataPtr(),    mu_g[mfi].dataPtr(),
          particle_state.dataPtr(), particle_phase.dataPtr(),
          des_radius.dataPtr(),     ro_sol.dataPtr(),
          pvol.dataPtr(),           pmass.dataPtr(),
          omoi.dataPtr(),           des_usr_var.dataPtr(),
          des_pos_new.dataPtr(),    des_vel_new.dataPtr(),   omega_new.dataPtr(),
          des_acc_old.dataPtr(),    rot_acc_old.dataPtr(),
          drag_fc.dataPtr(),        fc.dataPtr(),            tow.dataPtr(),
          pairs.dataPtr(),          &pair_count,
          flag[mfi].dataPtr(),
          &time, &dt, &nstep);
    }

    if(!steady_state) {
      time += prev_dt;
      nstep++;
    }

    for (MFIter mfi(flag); mfi.isValid(); ++mfi)
      mfix_output_manager(
        &time, &dt, &nstep,
        ep_g[mfi].dataPtr(),    p_g[mfi].dataPtr(),
        ro_g[mfi].dataPtr(),   rop_g[mfi].dataPtr(),
        u_g[mfi].dataPtr(),    v_g[mfi].dataPtr(),
        w_g[mfi].dataPtr(),
        particle_state.dataPtr(), des_radius.dataPtr(),
        ro_sol.dataPtr(), des_pos_new.dataPtr(),
        des_vel_new.dataPtr(), des_usr_var.dataPtr(),
        omega_new.dataPtr(), &estatus, &finish);

    // Mechanism to terminate MFIX normally.
    if(steady_state || time + 0.1L*dt >= tstop ||
       (solve_dem && !solve_fluid)) finish = 1;

  }while (finish==0);
#endif

  my_mfix.usr3(0);
#if 0
  for (MFIter mfi(flag); mfi.isValid(); ++mfi)
     mfix_usr3(u_g[mfi].dataPtr(),    v_g[mfi].dataPtr(),
               w_g[mfi].dataPtr(),    p_g[mfi].dataPtr());
#endif

  Real end_time = ParallelDescriptor::second() - strt_time;

  if (ParallelDescriptor::IOProcessor())
    std::cout << "Time spent in main " << end_time << std::endl;

  BoxLib::Finalize();
  return 0;
}
