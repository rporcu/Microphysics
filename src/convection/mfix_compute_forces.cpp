#include <mfix.H>

using namespace amrex;

void mfix::compute_tra_forces (Vector<MultiFab*> const& tra_forces,
                               Vector<MultiFab const*> const& density)
{
    // NOTE: this routine must return the force term for the update of (rho s), NOT just s.
    if (advect_tracer) {
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

void mfix::compute_vel_forces (amrex::Vector<amrex::MultiFab*      > const& vel_forces,
                               amrex::Vector<amrex::MultiFab const*> const& velocity,
                               amrex::Vector<amrex::MultiFab const*> const& density,
                               bool include_pressure_gradient)
{
    for (int lev = 0; lev <= finest_level; ++lev)
       compute_vel_forces_on_level (lev, *vel_forces[lev], *velocity[lev], *density[lev],
                                    include_pressure_gradient);
}

void mfix::compute_vel_forces_on_level (int lev,
                                          MultiFab& vel_forces,
                                          const MultiFab& /*velocity*/,
                                          const MultiFab& density,
                                          bool include_pressure_gradient)
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

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real rhoinv = 1.0/rho(i,j,k);

        if (include_pressure_gradient) {
          AMREX_D_TERM(vel_f(i,j,k,0) = l_gravity[0]-(gradp(i,j,k,0)+l_gp0[0])*rhoinv;,
                       vel_f(i,j,k,1) = l_gravity[1]-(gradp(i,j,k,1)+l_gp0[1])*rhoinv;,
                       vel_f(i,j,k,2) = l_gravity[2]-(gradp(i,j,k,2)+l_gp0[2])*rhoinv;);
        } else {
          AMREX_D_TERM(vel_f(i,j,k,0) = l_gravity[0]-(               l_gp0[0])*rhoinv;,
                       vel_f(i,j,k,1) = l_gravity[1]-(               l_gp0[1])*rhoinv;,
                       vel_f(i,j,k,2) = l_gravity[2]-(               l_gp0[2])*rhoinv;);
        }
      });
    }
}
