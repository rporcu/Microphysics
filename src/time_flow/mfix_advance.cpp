#include <mfix_F.H>
#include <mfix.H>

#include <AMReX_VisMF.H>
#include <MFIX_MFHelpers.H>
#include <MFIX_DEM_Parms.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

void
mfix::EvolveFluid (int nstep, Real& dt,  Real& time, Real stop_time, Real coupling_timing)
{
    BL_PROFILE_REGION_START("mfix::EvolveFluid");
    BL_PROFILE("mfix::EvolveFluid");

#ifdef AMREX_MEM_PROFILING
    {
      std::ostringstream ss;
      ss << "EvolveFluid Start";
      MemProfiler::report(ss.str());
    }
#endif

    amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";

    // Extrapolate boundary values for ro_g, tracer, ep_g and mu_g
    // The subsequent call to mfix_set_scalar_bcs will only overwrite
    // ep_g ghost values for PINF and POUT
    for (int lev = 0; lev <= finest_level; lev++)
    {
      m_leveldata[lev]->ro_g->FillBoundary(geom[lev].periodicity());
      m_leveldata[lev]->trac->FillBoundary(geom[lev].periodicity());
      m_leveldata[lev]->ep_g->FillBoundary(geom[lev].periodicity());
      m_leveldata[lev]->mu_g->FillBoundary(geom[lev].periodicity());
    }

    // Fill ghost nodes and reimpose boundary conditions
    //mfix_set_velocity_bcs(time, vel_g, 0);

    mfix_set_density_bcs(time, get_ro_g());
    //mfix_set_scalar_bcs(time, trac, mu_g);

    //
    // Start loop: if we are not seeking a steady state solution,
    // the loop will execute only once
    //
    int keep_looping = 1;
    int iter = 1;

    // Create temporary multifabs to hold the old-time conv and divtau
    //    so we don't have to re-compute them in the corrector
    Vector< MultiFab* > conv_u_old;
    Vector< MultiFab* > conv_s_old;
    Vector< MultiFab* > divtau_old;
    Vector< MultiFab* >   laps_old;

    conv_u_old.resize(finest_level+1);
    conv_s_old.resize(finest_level+1);
    divtau_old.resize(finest_level+1);
      laps_old.resize(finest_level+1);

    for (int lev = 0; lev <= finest_level; lev++)
    {
       conv_u_old[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
       conv_s_old[lev] = new MultiFab(grids[lev], dmap[lev], 2, 0, MFInfo(), *ebfactory[lev]);
       divtau_old[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
         laps_old[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), *ebfactory[lev]);

       conv_u_old[lev]->setVal(0.0);
       conv_s_old[lev]->setVal(0.0);
       divtau_old[lev]->setVal(0.0);
         laps_old[lev]->setVal(0.0);
    }

    do
    {
        mfix_compute_dt(nstep, time, stop_time, dt);

        // Set new and old time to correctly use in fillpatching
        for (int lev = 0; lev <= finest_level; lev++)
        {
            t_old[lev] = time;
            t_new[lev] = time+dt;
        }

        if (steady_state)
        {
           amrex::Print() << "\n   Iteration " << iter << " with dt = " << dt << "\n" << std::endl;
        } else {
           amrex::Print() << "\n   Step " << nstep+1 << ": from old_time " \
                          << time << " to new time " << time+dt
                          << " with dt = " << dt << "\n" << std::endl;
        }

        for (int lev = 0; lev <= finest_level; lev++)
        {
          MultiFab& ep_g = *m_leveldata[lev]->ep_g;
          MultiFab& ep_go = *m_leveldata[lev]->ep_go;

          MultiFab& p_g = *m_leveldata[lev]->p_g;
          MultiFab& p_go = *m_leveldata[lev]->p_go;

          MultiFab& ro_g = *m_leveldata[lev]->ro_g;
          MultiFab& ro_go = *m_leveldata[lev]->ro_go;

          MultiFab& trac = *m_leveldata[lev]->trac;
          MultiFab& trac_o = *m_leveldata[lev]->trac_o;

          MultiFab& vel_g = *m_leveldata[lev]->vel_g;
          MultiFab& vel_go = *m_leveldata[lev]->vel_go;

          // Back up field variables to old
          MultiFab::Copy(ep_go, ep_g, 0, 0, ep_g.nComp(), ep_go.nGrow());
          MultiFab::Copy(p_go, p_g, 0, 0, p_g.nComp(), p_go.nGrow());
          MultiFab::Copy(ro_go, ro_g, 0, 0, ro_g.nComp(), ro_go.nGrow());
          MultiFab::Copy(trac_o, trac, 0, 0, trac.nComp(), trac_o.nGrow());
          MultiFab::Copy(vel_go, vel_g, 0, 0, vel_g.nComp(), vel_go.nGrow());

          // User hooks
          for (MFIter mfi(ep_g, false); mfi.isValid(); ++mfi)
             mfix_usr2();
        }

        //
        // Time integration step
        //
        Real new_time = time+dt;

        // Calculate drag coefficient
        if (DEM::solve) {
          Real start_drag = ParallelDescriptor::second();
          mfix_calc_drag_fluid(time);
          coupling_timing += ParallelDescriptor::second() - start_drag;
        }

        // Predictor step
        bool proj_2_pred = true;
        mfix_apply_predictor(conv_u_old, conv_s_old, divtau_old, laps_old, time, dt, proj_2_pred);

        // Calculate drag coefficient
        if (DEM::solve)
        {
          Real start_drag = ParallelDescriptor::second();
          amrex::Print() << "\nRecalculating drag ..." << std::endl;
          mfix_calc_drag_fluid(new_time);
          coupling_timing += ParallelDescriptor::second() - start_drag;
        }

        bool proj_2_corr = true;
        // Corrector step
        if (!steady_state)
           mfix_apply_corrector(conv_u_old, conv_s_old, divtau_old, laps_old, time, dt, proj_2_corr);

        //
        // Check whether to exit the loop or not
        //
        if (steady_state) {
          keep_looping = !steady_state_reached ( dt, iter);
        } else {
          keep_looping = 0;
        }


        // Update iteration count
        ++iter;
    }
    while ( keep_looping );

    if (test_tracer_conservation)
    {
       amrex::Print() << "Sum tracer volume wgt = "
                      << volWgtSum(0, *m_leveldata[0]->trac, 0)
                      << " " << volEpsWgtSum(0, *m_leveldata[0]->trac, 0)
                      << std::endl;
    }

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "EvolveFluid Stop";
            MemProfiler::report(ss.str());
        }
#endif

    for (int lev = 0; lev <= finest_level; lev++)
    {
       delete conv_u_old[lev];
       delete conv_s_old[lev];
       delete divtau_old[lev];
       delete   laps_old[lev];
    }

    BL_PROFILE_REGION_STOP("mfix::EvolveFluid");
}

void
mfix::mfix_project_velocity ()
{
    // Project velocity field to make sure initial velocity is divergence-free
    Real dummy_dt = 1.0;

    amrex::Print() << "Initial projection:\n";

    bool proj_2 = true;
    Real time = 0.0;

    // Apply projection -- depdt=0 for now
    Vector< MultiFab* > depdt(finest_level+1);
    for (int lev(0); lev <= finest_level; ++lev)
      depdt[lev] = MFHelpers::createFrom(*(m_leveldata[lev]->ep_g), 0.0, 1).release();

    mfix_apply_nodal_projection(depdt, time, dummy_dt, proj_2);

    for (int lev(0); lev <= finest_level; ++lev)
      delete depdt[lev];

    // We initialize p_g and gp back to zero (p0_g may still be still non-zero)
    for (int lev = 0; lev <= finest_level; lev++)
    {
      m_leveldata[lev]->p_g->setVal(0);
      m_leveldata[lev]->gp->setVal(0);
    }
}

void
mfix::mfix_initial_iterations (Real dt, Real stop_time)
{
  Real time = 0.0;
  int nstep = 0;

  mfix_compute_dt(nstep,time,stop_time,dt);

  amrex::Print() << "Doing initial pressure iterations with dt = " << dt << "\n";

  // Fill ghost nodes and reimpose boundary conditions
  mfix_set_velocity_bcs(time, get_vel_g(), 0);
  mfix_set_density_bcs(time, get_ro_g());
  mfix_set_scalar_bcs(time, get_trac(), get_mu_g());

  // Copy vel_g into vel_go
  for (int lev = 0; lev <= finest_level; lev++)
    MultiFab::Copy(*m_leveldata[lev]->vel_go, *m_leveldata[lev]->vel_g, 0, 0,
                   m_leveldata[lev]->vel_g->nComp(), m_leveldata[lev]->vel_go->nGrow());

  if (DEM::solve)
    mfix_calc_drag_fluid(time);

  // Create temporary multifabs to hold conv and divtau
  Vector< MultiFab* > conv_u(finest_level+1, nullptr);
  Vector< MultiFab* > conv_s(finest_level+1, nullptr);
  Vector< MultiFab* > divtau(finest_level+1, nullptr);
  Vector< MultiFab* >   laps(finest_level+1, nullptr);

  for (int lev = 0; lev <= finest_level; lev++)
  {
    conv_u[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
    conv_s[lev] = new MultiFab(grids[lev], dmap[lev], 2, 0, MFInfo(), *ebfactory[lev]);
    divtau[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
    laps[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), *ebfactory[lev]);

    conv_u[lev]->setVal(0.0);
    conv_s[lev]->setVal(0.0);
    divtau[lev]->setVal(0.0);
    laps[lev]->setVal(0.0);
  }

  for (int iter = 0; iter < initial_iterations; ++iter)
  {
    amrex::Print() << " " << std::endl;
    amrex::Print() << "In initial_iterations: iter = " << iter << "\n";

    bool proj_2 = false;

    mfix_apply_predictor(conv_u, conv_s, divtau, laps, time, dt, proj_2);

    // Reset any quantities which might have been updated
    for (int lev = 0; lev <= finest_level; lev++)
      MultiFab::Copy(*m_leveldata[lev]->vel_g, *m_leveldata[lev]->vel_go, 0, 0,
                     m_leveldata[lev]->vel_g->nComp(), m_leveldata[lev]->vel_g->nGrow());

    if (advect_density)
      for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Copy(*m_leveldata[lev]->ro_g, *m_leveldata[lev]->ro_go, 0, 0,
                        m_leveldata[lev]->ro_g->nComp(), m_leveldata[lev]->ro_g->nGrow());

    if (advect_tracer)
      for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Copy(*m_leveldata[lev]->trac, *m_leveldata[lev]->trac_o, 0, 0,
                       m_leveldata[lev]->trac->nComp(), m_leveldata[lev]->trac->nGrow());

    // Reset the boundary values (necessary if they are time-dependent)
    mfix_set_velocity_bcs(time, get_vel_g(), 0);
    mfix_set_density_bcs(time, get_ro_g());
    mfix_set_scalar_bcs(time, get_trac(), get_mu_g());
  }

  for (int lev = 0; lev <= finest_level; lev++)
  {
     delete conv_u[lev];
     delete conv_s[lev];
     delete divtau[lev];
     delete   laps[lev];
  }
}

//
// Explicit solve for the intermediate velocity.
// Currently this means accounting for the implicit part of the fluid/particle
// momentum exchange
//
void
mfix::mfix_add_drag_explicit (Real dt)
{
  /*
     This adds both components of the drag term
     So the drag term we add is beta * (particle_velocity - fluid_velocity)
                              = drag(0:2) - drag(3) * fluid_velocity
  */

  BL_PROFILE("mfix::mfix_add_drag_explicit");

  for (int lev = 0; lev <= finest_level; lev++)
  {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_leveldata[lev]->vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      // Tilebox
      Box bx = mfi.tilebox();

      const auto&  vel_fab = m_leveldata[lev]->vel_g->array(mfi);
      const auto& drag_fab = m_leveldata[lev]->drag->array(mfi);
      const auto&   ro_fab = m_leveldata[lev]->ro_g->array(mfi);
      const auto&   ep_fab = m_leveldata[lev]->ep_g->array(mfi);

      amrex::ParallelFor(bx,[dt,vel_fab,drag_fab,ro_fab,ep_fab]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
          Real orop  = dt / (ro_fab(i,j,k) * ep_fab(i,j,k));

          Real drag_0 = (drag_fab(i,j,k,0) - drag_fab(i,j,k,3)*vel_fab(i,j,k,0)) * orop;
          Real drag_1 = (drag_fab(i,j,k,1) - drag_fab(i,j,k,3)*vel_fab(i,j,k,1)) * orop;
          Real drag_2 = (drag_fab(i,j,k,2) - drag_fab(i,j,k,3)*vel_fab(i,j,k,2)) * orop;

          vel_fab(i,j,k,0) += drag_0;
          vel_fab(i,j,k,1) += drag_1;
          vel_fab(i,j,k,2) += drag_2;
      });
    }
  }
}

//
// Implicit solve for the intermediate velocity.
// Currently this means accounting for the implicit part of the fluid/particle
// momentum exchange
//
void
mfix::mfix_add_drag_implicit (Real dt)
{
  /*
     This adds both components of the drag term
     So the drag term we add is beta * (particle_velocity - fluid_velocity)
                              = drag(0:2) - drag(3) * fluid_velocity
  */

  BL_PROFILE("mfix::mfix_add_drag_implicit");

  for (int lev = 0; lev <= finest_level; lev++)
  {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_leveldata[lev]->vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      // Tilebox
      Box bx = mfi.tilebox();

      const auto&  vel_fab = m_leveldata[lev]->vel_g->array(mfi);
      const auto& drag_fab = m_leveldata[lev]->drag->array(mfi);
      const auto&   ro_fab = m_leveldata[lev]->ro_g->array(mfi);
      const auto&   ep_fab = m_leveldata[lev]->ep_g->array(mfi);

      amrex::ParallelFor(bx,[dt,vel_fab,drag_fab,ro_fab,ep_fab]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
          Real orop  = dt / (ro_fab(i,j,k) * ep_fab(i,j,k));
          Real denom = 1.0 / (1.0 + drag_fab(i,j,k,3) * orop);

          vel_fab(i,j,k,0) = (vel_fab(i,j,k,0) + drag_fab(i,j,k,0) * orop) * denom;
          vel_fab(i,j,k,1) = (vel_fab(i,j,k,1) + drag_fab(i,j,k,1) * orop) * denom;
          vel_fab(i,j,k,2) = (vel_fab(i,j,k,2) + drag_fab(i,j,k,2) * orop) * denom;
      });
    }
  }
}

//
// Check if steady state has been reached by verifying that
//
//      max(abs( u^(n+1) - u^(n) )) < tol * dt
//      max(abs( v^(n+1) - v^(n) )) < tol * dt
//      max(abs( w^(n+1) - w^(n) )) < tol * dt
//

int
mfix::steady_state_reached (Real dt, int iter)
{
    //
    // Count number of access
    //
    static int naccess = 0;

    int condition1[finest_level+1];
    int condition2[finest_level+1];

    Real time = 0.;

    mfix_set_velocity_bcs(time, get_vel_g(), 0);

    //
    // Make sure velocity is up to date
    //
    for (int lev = 0; lev <= finest_level; lev++)
    {

       //
       // Use temporaries to store the difference
       // between current and previous solution
       //
       MultiFab temp_vel(m_leveldata[lev]->vel_g->boxArray(),
                         m_leveldata[lev]->vel_g->DistributionMap(), 3, 0);
       MultiFab::LinComb(temp_vel, 1.0, *m_leveldata[lev]->vel_g, 0, -1.0,
                         *m_leveldata[lev]->vel_go, 0, 0, 3, 0);

       MultiFab tmp;

       const BoxArray & nd_grid = amrex::convert(grids[lev], IntVect{1,1,1});
       tmp.define(nd_grid, dmap[lev], 1, 0);

       MultiFab::LinComb(tmp, 1.0, *m_leveldata[lev]->p_g, 0, -1.0,
                         *m_leveldata[lev]->p_go, 0, 0, 1, 0);

       Real delta_u = temp_vel.norm0(0,0,false,true);
       Real delta_v = temp_vel.norm0(1,0,false,true);
       Real delta_w = temp_vel.norm0(2,0,false,true);
       Real delta_p = temp_vel.norm0(0,0,false,true);

       Real tol = steady_state_tol;

       condition1[lev] = (delta_u < tol*dt) and (delta_v < tol*dt ) and (delta_w < tol*dt);

       //
       // Second stop condition
       //
       Periodicity period = geom[lev].periodicity();

       Real du_n1 = temp_vel.norm1(0,period,true);
       Real dv_n1 = temp_vel.norm1(1,period,true);
       Real dw_n1 = temp_vel.norm1(2,period,true);
       Real dp_n1 =      tmp.norm1(0,period,true);
       Real uo_n1 = m_leveldata[lev]->vel_go->norm1(0,period,true);
       Real vo_n1 = m_leveldata[lev]->vel_go->norm1(1,period,true);
       Real wo_n1 = m_leveldata[lev]->vel_go->norm1(2,period,true);
       Real po_n1 = m_leveldata[lev]->p_go->norm1(0,period,true);

       Real tmp1, tmp2, tmp3, tmp4;

       Real local_tol = 1.0e-8;

       if ( uo_n1 < local_tol ) {
          tmp1 = 0.0;
       } else {
          tmp1 = du_n1 / uo_n1;
       };

       if ( vo_n1 < local_tol ) {
          tmp2 = 0.0;
       } else {
          tmp2 = dv_n1 / vo_n1;
       };

       if ( wo_n1 < local_tol ) {
          tmp3 = 0.0;
       } else {
          tmp3 = dw_n1 / wo_n1;
       };

       if ( po_n1 < local_tol ) {
          tmp4 = 0.0;
       } else {
          tmp4 = dp_n1 / po_n1;
       };

       condition2[lev] = (tmp1 < tol) and (tmp2 < tol) and (tmp3 < tol); // && (tmp4 < tol);

       //
       // Print out info on steady state checks
       //
       amrex::Print() << "\nSteady state check at level " << lev << ":\n";
       amrex::Print() << "||u-uo||/||uo|| , du/dt  = " << tmp1 <<" , "<< delta_u/dt << "\n";
       amrex::Print() << "||v-vo||/||vo|| , dv/dt  = " << tmp2 <<" , "<< delta_v/dt << "\n";
       amrex::Print() << "||w-wo||/||wo|| , dw/dt  = " << tmp3 <<" , "<< delta_w/dt << "\n";
       amrex::Print() << "||p-po||/||po|| , dp/dt  = " << tmp4 <<" , "<< delta_p/dt << "\n";
    }

    int reached = 1;
    for (int lev = 0; lev <= finest_level; lev++)
    {
       reached = reached and (condition1[lev] or condition2[lev]);
    }

    reached = reached or (iter >= steady_state_maxiter);

    // Count # access
    naccess++;

    //
    //  Always return negative to first access. This way
    //  initial zero velocity field do not test for false positive
    //
    if ( naccess == 1 ) {
       return 0;
    } else {
       return reached;
    };
}
