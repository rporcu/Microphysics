#include "mfix_proj_F.H"
#include <mfix.H>

#include <cmath>
#include <limits>

using namespace std;

void
mfix::mfix_compute_dt(int nstep, Real time, Real stop_time, Real& dt)
{
    // dt is always computed even when fixed_dt is set, 
    // so we can issue a warning if the value of fixed dt does not satisfy the CFL condition.

    Real dt_new;
    Real old_dt = dt;

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

#ifdef AMREX_USE_CUDA
    bool notInLaunchRegionStatus = Gpu::notInLaunchRegion();

    if(notInLaunchRegionStatus == true)
      Gpu::setLaunchRegion(true);

    {
      Gpu::DeviceScalar<Real> cfl_max_gpu(cfl_max);
      Real *cfl_max_ptr = cfl_max_gpu.dataPtr();
#endif
    
    for (int lev(0); lev < nlev; ++lev) {

        const Real* dx = geom[lev].CellSize();

        Real odx(1.0 / dx[0]);
        Real ody(1.0 / dx[1]);
        Real odz(1.0 / dx[2]);

#ifdef _OPENMP
#pragma omp parallel reduction(max:cfl_max) if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*vel_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const auto& vel       =   vel_g[lev] -> array(mfi);
            const auto& ep        =    ep_g[lev] -> array(mfi);
            const auto& ro        =    ro_g[lev] -> array(mfi);
            const auto& mu        =    mu_g[lev] -> array(mfi);
            const auto& gradp     =      gp[lev] -> array(mfi);
            const auto& drag_fab  =    drag[lev] -> array(mfi);
            
            Box bx(mfi.tilebox());

            const auto&  vel_fab   = static_cast<EBFArrayBox const&>((*vel_g[lev])[mfi]);
            const auto&  flags     = vel_fab.getEBCellFlagFab();
            const auto&  flags_fab = flags.array();

            const GpuArray<Real, 3> gp0_dev = {gp0[0], gp0[1], gp0[2]};
            const GpuArray<Real, 3> gravity_dev = {gravity[0], gravity[1], gravity[2]};

            // Compute CFL on a per cell basis
            if (flags.getType(bx) != FabType::covered)
            {
                AMREX_FOR_3D(bx, i, j, k,                                    
                {
                    if (!flags_fab(i,j,k).isCovered()) {
                        
                        Real acc[3];
                        Real qro  = 1.0/ro(i,j,k);
                        Real qep  = 1.0/ep(i,j,k);

                        // Compute the three components of the net acceleration
                        // Explicit particle forcing is given by 
                        for (int n(0); n < 3; ++n) {
                            Real delp = gp0_dev[n] + gradp(i,j,k,n);
                            Real fp   = drag_fab(i,j,k,n) - drag_fab(i,j,k,3) * vel(i,j,k,n);
                            
                            acc[n] = gravity_dev[n] + qro * ( - delp + fp*qep );
                        }
                        
                        Real c_cfl   = abs(vel(i,j,k,0))*odx + abs(vel(i,j,k,1))*ody + abs(vel(i,j,k,2))*odz;                        
                        Real v_cfl   = 2.0 * mu(i,j,k) * qro * (odx*odx + ody*ody + odz*odz);
                        Real cpv_cfl = c_cfl + v_cfl;

                        // MAX CFL factor on cell (i,j,k)
                        Real cfl_max_cell = cpv_cfl + std::sqrt( cpv_cfl*cpv_cfl +
                                                                 4.0*abs(acc[0])*odx  +
                                                                 4.0*abs(acc[1])*ody  +
                                                                 4.0*abs(acc[2])*odz  );
#ifdef AMREX_USE_CUDA
                        Gpu::Atomic::Max(cfl_max_ptr, cfl_max_cell);
#else
                        cfl_max = std::max(cfl_max, cfl_max_cell);
#endif
                    }
                });
            }
        }
    }

#ifdef AMREX_USE_CUDA
      cfl_max = cfl_max_gpu.dataValue();
    }

    if(notInLaunchRegionStatus == true)
      Gpu::setLaunchRegion(notInLaunchRegionStatus);
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
    if (nstep > 1) dt_new = std::min( dt_new, 1.01*old_dt );

    // Don't overshoot the final time if not running to steady state
    if (steady_state == 0 && stop_time > 0.) 
    {
       if (time+dt_new > stop_time) dt_new = stop_time - time;
    }

    // dt_new is the step calculated with a cfl contraint; dt is the value set by fixed_dt
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
