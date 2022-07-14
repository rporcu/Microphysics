#include <mfix.H>

using namespace amrex;

void mfix::compute_tra_forces (Vector<MultiFab*> const& tra_forces,
                               Vector<MultiFab const*> const& density)
{
    // NOTE: this routine must return the force term for the update of (rho s), NOT just s.
    if (fluid.solve_tracer()) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (MFIter mfi(*tra_forces[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real>       const& tra_f = tra_forces[lev]->array(mfi);
                Array4<Real const> const& rho   =    density[lev]->const_array(mfi);

                amrex::ParallelFor(bx, ntrac,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    // For now we don't have any external forces on the scalars
                    tra_f(i,j,k,n) = 0.0;

                    // Return the force term for the update of (rho s), NOT just s.
                    tra_f(i,j,k,n) *= rho(i,j,k);
                });
            }
        }
    }
}

void mfix::compute_vel_forces (Vector<MultiFab*      > const& vel_forces,
                               Vector<MultiFab const*> const& velocity,
                               Vector<MultiFab const*> const& density,
                               Vector<MultiFab const*> const& txfr_in,
                               bool include_pressure_gradient,
                               bool include_drag_force)
{
  if ( m_verbose > 1 ) {
    if ( include_pressure_gradient ) {
      amrex::Print() << "\nIncluding pressure gradient in vel forces\n";
    } else {
      amrex::Print() << "\nNOT including pressure gradient in vel forces\n";
    }

    if ( include_drag_force ) {
      amrex::Print() << "Including drag force in vel forces\n\n";
    } else {
      amrex::Print() << "NOT Including drag force in vel forces\n\n";
    }
  }
  for (int lev = 0; lev <= finest_level; ++lev) {
     compute_vel_forces_on_level (lev, *vel_forces[lev], *velocity[lev], *density[lev],
                                  *txfr_in[lev], include_pressure_gradient,
                                  include_drag_force);
  }
}

void mfix::compute_vel_forces_on_level (int lev,
                                        MultiFab& vel_forces,
                                        const MultiFab& velocity,
                                        const MultiFab& density,
                                        const MultiFab& txfr_in,
                                        bool include_pressure_gradient,
                                        bool include_drag_force)
{
    GpuArray<Real,3> l_gravity{gravity[0],gravity[1],gravity[2]};
    GpuArray<Real,3> l_gp0{gp0[0], gp0[1], gp0[2]};

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel_forces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.tilebox();
      Array4<Real>       const& vel_f =  vel_forces.array(mfi);
      Array4<Real const> const&   rho =     density.const_array(mfi);
      Array4<Real const> const& gradp = m_leveldata[lev]->gp->const_array(mfi);

      if(!include_pressure_gradient && !include_drag_force){
        amrex::ParallelFor(bx, [vel_f, rho, l_gravity, gradp, l_gp0]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const Real rhoinv = 1.0/rho(i,j,k);

          vel_f(i,j,k,0) = l_gravity[0]-(               l_gp0[0])*rhoinv;
          vel_f(i,j,k,1) = l_gravity[1]-(               l_gp0[1])*rhoinv;
          vel_f(i,j,k,2) = l_gravity[2]-(               l_gp0[2])*rhoinv;

        });

      } else if (include_pressure_gradient && !include_drag_force) {
        amrex::ParallelFor(bx, [vel_f, rho, l_gravity, gradp, l_gp0]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const Real rhoinv = 1.0/rho(i,j,k);

          vel_f(i,j,k,0) = l_gravity[0]-(gradp(i,j,k,0)+l_gp0[0])*rhoinv;
          vel_f(i,j,k,1) = l_gravity[1]-(gradp(i,j,k,1)+l_gp0[1])*rhoinv;
          vel_f(i,j,k,2) = l_gravity[2]-(gradp(i,j,k,2)+l_gp0[2])*rhoinv;
        });

      } else if (include_pressure_gradient && include_drag_force) {

        Array4<Real const> const& vel_g = velocity.const_array(mfi);
        Array4<Real const> const& txfr  = txfr_in.const_array(mfi);
        Array4<Real const> const& ep_g  = m_leveldata[lev]->ep_g->const_array(mfi);

        amrex::ParallelFor(bx, [vel_f, rho, l_gravity, gradp, l_gp0, vel_g, ep_g, txfr]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const Real rhoinv = 1.0/rho(i,j,k);
          const Real epginv = 1.0/ep_g(i,j,k);

          const Real beta = txfr(i,j,k,Transfer::drag_coeff);

          const Real drag_x = (txfr(i,j,k,Transfer::velx) - beta*vel_g(i,j,k,0))*epginv;
          const Real drag_y = (txfr(i,j,k,Transfer::vely) - beta*vel_g(i,j,k,1))*epginv;
          const Real drag_z = (txfr(i,j,k,Transfer::velz) - beta*vel_g(i,j,k,2))*epginv;

          vel_f(i,j,k,0) = l_gravity[0]-(gradp(i,j,k,0)+l_gp0[0]+drag_x)*rhoinv;
          vel_f(i,j,k,1) = l_gravity[1]-(gradp(i,j,k,1)+l_gp0[1]+drag_y)*rhoinv;
          vel_f(i,j,k,2) = l_gravity[2]-(gradp(i,j,k,2)+l_gp0[2]+drag_z)*rhoinv;
        });

      } else {
        amrex::Abort("Bad combo of inputs in compute vel forces");
      }

    }
}
