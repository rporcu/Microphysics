#include <mfix.H>

#include <limits>

using namespace std;

void
mfix::mfix_compute_dt (int nstep, Real time, Real stop_time, Real& dt, Real& prev_dt)
{
    // dt is always computed even when fixed_dt is set,
    // so we can issue a warning if the value of fixed dt does not satisfy the CFL condition.

    Real dt_new;
    Real old_dt = dt;

    // Store the dt we've just used in the previous time step as prev_dt
    prev_dt = dt;

    /*
       Compute new dt by using the formula derived in
       "A Boundary Condition Capturing Method for Multiphase Incompressible Flow"
       by Kang et al. (JCP).

       dt/2 * ( C+V + sqrt( (C+V)**2 + 4Fx/dx + 4Fy/dy + 4Fz/dz )

      where

      C = max(|U|)/dx + max(|V|)/dy + max(|W|)/dz    --> Convection

      V = 2 * max(mu/ro) * (1/dx^2 + 1/dy^2 +1/dz^2) --> Diffusion

      Fx, Fy, Fz = net acceleration due to external forces

    */

    // Max CFL factor for all levels
    Real cfl_max(0.0);

#ifdef AMREX_USE_GPU
    {
      Gpu::DeviceScalar<Real> cfl_max_gpu(cfl_max);
      Real *cfl_max_ptr = cfl_max_gpu.dataPtr();
#endif

    for (int lev(0); lev <= finest_level; ++lev) {

        const Real* dx = geom[lev].CellSize();

        Real odx(1.0 / dx[0]);
        Real ody(1.0 / dx[1]);
        Real odz(1.0 / dx[2]);

#ifdef _OPENMP
#pragma omp parallel reduction(max:cfl_max) if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*m_leveldata[lev]->vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const auto& vel       = m_leveldata[lev]->vel_g->array(mfi);
            const auto& ep        = m_leveldata[lev]->ep_g->array(mfi);
            const auto& ro        = m_leveldata[lev]->ro_g->array(mfi);
            const auto& mu        = m_leveldata[lev]->mu_g->array(mfi);
            const auto& gradp     = m_leveldata[lev]->gp->array(mfi);
            const auto& txfr_fab  = m_leveldata[lev]->txfr->array(mfi);

            Box bx(mfi.tilebox());

            const auto& vel_fab   =
              static_cast<EBFArrayBox const&>((*m_leveldata[lev]->vel_g)[mfi]);

            const auto& flags     = vel_fab.getEBCellFlagFab();
            const auto& flags_fab = flags.array();

            // ew need this until we remove static attribute from mfix::gp0
            const RealVect gp0_dev(gp0);
            const RealVect gravity_dev(gravity);

            // Compute CFL on a per cell basis
            if (flags.getType(bx) != FabType::covered)
            {
              amrex::ParallelFor(bx,
                  [ro,ep,gp0_dev,gradp,txfr_fab,gravity_dev,vel,odx,ody,odz,flags_fab,mu,
#ifdef AMREX_USE_GPU
                  cfl_max_ptr]
#else
                  &cfl_max]
#endif
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                  {
                    if (!flags_fab(i,j,k).isCovered())
                    {
                        RealVect acc(0.);
                        Real qro  = 1.0/ro(i,j,k);
                        Real qep  = 1.0/ep(i,j,k);

                        // Compute the three components of the net acceleration
                        // Explicit particle forcing is given by
                        for (int n(0); n < 3; ++n) {
                            Real delp = gp0_dev[n] + gradp(i,j,k,n);
                            Real fp   = txfr_fab(i,j,k,n) - txfr_fab(i,j,k,3) * vel(i,j,k,n);

                            acc[n] = gravity_dev[n] + qro * ( - delp + fp*qep );
                        }

                        Real c_cfl   = amrex::Math::abs(vel(i,j,k,0))*odx +
                                       amrex::Math::abs(vel(i,j,k,1))*ody +
                                       amrex::Math::abs(vel(i,j,k,2))*odz;
                        Real v_cfl   = 2.0 * mu(i,j,k) * qro * (odx*odx + ody*ody + odz*odz);
                        Real cpv_cfl = c_cfl + v_cfl;

                        // MAX CFL factor on cell (i,j,k)
                        Real cfl_max_cell = cpv_cfl + std::sqrt( cpv_cfl*cpv_cfl +
                                                                 4.0*amrex::Math::abs(acc[0])*odx  +
                                                                 4.0*amrex::Math::abs(acc[1])*ody  +
                                                                 4.0*amrex::Math::abs(acc[2])*odz  );
#ifdef AMREX_USE_GPU
                        Gpu::Atomic::Max(cfl_max_ptr, cfl_max_cell);
#else
                        cfl_max = amrex::max(cfl_max, cfl_max_cell);
#endif
                    }
                  });
            }
        }
    }

#ifdef AMREX_USE_GPU
      cfl_max = cfl_max_gpu.dataValue();
    }
#endif

    // Do global max operation
    ParallelDescriptor::ReduceRealMax(cfl_max);

    // New dt
    dt_new = cfl * 2.0 / cfl_max;

    // Protect against cfl_max very small
    // This may happen, for example, when the initial velocity field
    // is zero for an inviscid flow with no external forcing
    Real eps = numeric_limits<Real>::epsilon();
    if ( nstep > 1 && cfl_max <= eps ) dt_new = 0.5 * old_dt;

    // Don't let the timestep grow by more than 1% per step.
    //       unless the previous time step was unduly shrunk to match plot_per_exact
    if ( nstep > 1 && !(plot_per_exact > 0 && last_plt == nstep && nstep > 0) )
        dt_new = amrex::min( dt_new, 1.01*old_dt );

    // Don't overshoot the final time if not running to steady state
    if (steady_state == 0 && stop_time > 0.)
       if (time+dt_new > stop_time)
           dt_new = stop_time - time;

    // Don't overshoot specified plot times
    if(plot_per_exact > 0.0 &&
            (trunc((time + dt_new + eps) / plot_per_exact) > trunc((time + eps) / plot_per_exact)))
    {
        dt_new = trunc((time + dt_new) / plot_per_exact) * plot_per_exact - time;
    }

    // dt_new is the step calculated with a cfl constraint; dt is the value set by fixed_dt
    // When the test was on dt > dt_new, there were cases where they were effectively equal
    //   but (dt > dt_new) was being set to true due to precision issues.
    Real ope(1.0 + 1.e-8);

    if ( fixed_dt > 0.)
    {
        if ( fixed_dt > dt_new*ope && cfl > 0)
        {
            amrex::Print() << "WARNING: fixed dt does not satisfy CFL condition: "
                           << " fixed dt = "  << fixed_dt
                           << " > dt based on cfl: " << dt_new
                           << std::endl;
            amrex::Abort ("Fixed dt is too large for fluid solve");
        } else {

            dt = fixed_dt;

        }
    }
    else
    {
        dt = amrex::min( dt_new, dt_max );

        if ( dt < dt_min )
            amrex::Abort ("Current dt is smaller than dt_min");

    }
}
