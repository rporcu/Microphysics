#include <mfix_F.H>
#include <mfix.H>

#include <AMReX_VisMF.H>
#include <MFIX_MFHelpers.H>
#include <MFIX_DEM_Parms.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

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
  mfix_set_temperature_bcs(time, get_T_g());
  mfix_set_scalar_bcs(time, get_trac(), get_mu_g());

  // Copy vel_g into vel_go
  for (int lev = 0; lev <= finest_level; lev++)
    MultiFab::Copy(*m_leveldata[lev]->vel_go, *m_leveldata[lev]->vel_g, 0, 0,
                   m_leveldata[lev]->vel_g->nComp(), m_leveldata[lev]->vel_go->nGrow());

  if (DEM::solve)
    mfix_calc_drag_fluid(time);

  // Create temporary multifabs to hold conv and divtau
  Vector< MultiFab* >  conv_u(finest_level+1, nullptr);
  Vector< MultiFab* >  conv_s(finest_level+1, nullptr);
  Vector< MultiFab* >  divtau(finest_level+1, nullptr);
  Vector< MultiFab* >    laps(finest_level+1, nullptr);
  Vector< MultiFab* > laptemp(finest_level+1, nullptr);

  for (int lev = 0; lev <= finest_level; lev++)
  {
    conv_u[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
    // TODO: check that 3 components is correct since now we have one for
    // density, one for tracer and one for temperature
    conv_s[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
    divtau[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
    laps[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), *ebfactory[lev]);
    laptemp[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);

    conv_u[lev]->setVal(0.0);
    conv_s[lev]->setVal(0.0);
    divtau[lev]->setVal(0.0);
    laps[lev]->setVal(0.0);
    laptemp[lev]->setVal(0.0);
  }

  for (int iter = 0; iter < initial_iterations; ++iter)
  {
    amrex::Print() << " " << std::endl;
    amrex::Print() << "In initial_iterations: iter = " << iter << "\n";

    bool proj_2 = false;

    mfix_apply_predictor(conv_u, conv_s, divtau, laps, laptemp, time, dt, proj_2);

    // Reset any quantities which might have been updated
    for (int lev = 0; lev <= finest_level; lev++)
      MultiFab::Copy(*m_leveldata[lev]->vel_g, *m_leveldata[lev]->vel_go, 0, 0,
                     m_leveldata[lev]->vel_g->nComp(), m_leveldata[lev]->vel_g->nGrow());

    if (advect_density)
      for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Copy(*m_leveldata[lev]->ro_g, *m_leveldata[lev]->ro_go, 0, 0,
                        m_leveldata[lev]->ro_g->nComp(), m_leveldata[lev]->ro_g->nGrow());

    if (advect_temperature)
      for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Copy(*m_leveldata[lev]->T_g, *m_leveldata[lev]->T_go, 0, 0,
                       m_leveldata[lev]->T_g->nComp(), m_leveldata[lev]->T_g->nGrow());

    if (advect_tracer)
      for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Copy(*m_leveldata[lev]->trac, *m_leveldata[lev]->trac_o, 0, 0,
                       m_leveldata[lev]->trac->nComp(), m_leveldata[lev]->trac->nGrow());

    // Reset the boundary values (necessary if they are time-dependent)
    mfix_set_velocity_bcs(time, get_vel_g(), 0);
    mfix_set_density_bcs(time, get_ro_g());
    mfix_set_temperature_bcs(time, get_T_g());
    mfix_set_scalar_bcs(time, get_trac(), get_mu_g());
  }

  for (int lev = 0; lev <= finest_level; lev++)
  {
     delete  conv_u[lev];
     delete  conv_s[lev];
     delete  divtau[lev];
     delete    laps[lev];
     delete laptemp[lev];
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
        const Real orop  = dt / (ro_fab(i,j,k) * ep_fab(i,j,k));

        const Real A = drag_fab(i,j,k,3);
        const Real vel_x = vel_fab(i,j,k,0);
        const Real vel_y = vel_fab(i,j,k,1);
        const Real vel_z = vel_fab(i,j,k,2);

        const Real drag_0 = (drag_fab(i,j,k,0) - A*vel_x) * orop;
        const Real drag_1 = (drag_fab(i,j,k,1) - A*vel_y) * orop;
        const Real drag_2 = (drag_fab(i,j,k,2) - A*vel_z) * orop;

        vel_fab(i,j,k,0) = vel_x + drag_0;
        vel_fab(i,j,k,1) = vel_y + drag_1;
        vel_fab(i,j,k,2) = vel_z + drag_2;
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
