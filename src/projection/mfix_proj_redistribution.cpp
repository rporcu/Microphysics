#include <mfix.H>
#include <hydro_redistribution.H>

#include <mfix_fluid_parms.H>
void
mfix::PostProjectionRedistribution (Real l_time, Real l_dt, 
                                    const Vector<MultiFab*>& sigma)
{
    if ( m_redistribution_type == "StateRedist")
    {
      // We must fill internal ghost values before calling redistribution
      // We also need any physical boundary conditions imposed if we are
      //    calling state redistribution (because that calls the slope routine)

      // We really only need to fillpatch velocity but what the heck -- for now 
      //    we'll do them all
      fillpatch_all(get_vel_g(), get_ro_g(), get_h_g(),
                    get_trac(), get_X_gk(), l_time);

      for (int lev = 0; lev <= finest_level; lev++)
      {
        auto& ld = *m_leveldata[lev];

        int ncomp = AMREX_SPACEDIM;

        // Make a temporary to the redistributed velocity into
        MultiFab new_vel(grids[lev], dmap[lev], ncomp, 0);
        new_vel.setVal(0.);

        auto const& bc_vel = get_hydro_velocity_bcrec_device_ptr();

        for (MFIter mfi(*ld.vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& fact = EBFactory(lev);

            EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
            Array4<EBCellFlag const> const& flag = flagfab.const_array();

            Array4<Real> const& vel_redist = new_vel.array(mfi);
            Array4<Real> const& vel_orig   = ld.vel_g->array(mfi);
            Array4<Real> const& gradp      = ld.gp->array(mfi);
            Array4<Real> const& ep_g       = ld.ep_g->array(mfi);
            Array4<Real> const& sig        = sigma[lev]->array(mfi);

            if ( (flagfab.getType(amrex::grow(bx,1)) != FabType::covered) &&
                 (flagfab.getType(amrex::grow(bx,1)) != FabType::regular) )
            {
                Array4<Real const> fcx, fcy, fcz, ccc, vfrac, apx, apy, apz;
                fcx = fact.getFaceCent()[0]->const_array(mfi);
                fcy = fact.getFaceCent()[1]->const_array(mfi);
                fcz = fact.getFaceCent()[2]->const_array(mfi);
                ccc   = fact.getCentroid().const_array(mfi);
                apx = fact.getAreaFrac()[0]->const_array(mfi);
                apy = fact.getAreaFrac()[1]->const_array(mfi);
                apz = fact.getAreaFrac()[2]->const_array(mfi);
                vfrac = fact.getVolFrac().const_array(mfi);

                Redistribution::ApplyToInitialData( bx, ncomp, vel_redist, vel_orig,
                                                    flag, apx, apy, apz, vfrac,
                                                    fcx, fcy, fcz, ccc, bc_vel,
                                                    geom[lev],
                                                    (m_redistribution_type == "StateRedist") ? 
                                                         "NewStateRedist" : m_redistribution_type);

                // We update gradp so that (vel_redist + dt gradp_redistnew/rho) == (vel_orig + dt gradp_orig/rho)
                // Note that we do not change rho in the redistribution
                amrex::ParallelFor(bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    Real delta_vel = vel_redist(i,j,k,n) - vel_orig(i,j,k,n);
                    gradp(i,j,k,n) -= delta_vel * ep_g(i,j,k) / (l_dt * sig(i,j,k));
                });

            } else {
                amrex::ParallelFor(bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    vel_redist(i,j,k,n) = vel_orig(i,j,k,n);
                });
            }
        }

        // Copy back into the original MultiFab
        MultiFab::Copy(*ld.vel_g, new_vel, 0,0, AMREX_SPACEDIM, 0);

        // We fill internal ghost values after calling redistribution
        ld.vel_g->FillBoundary();
        ld.gp->FillBoundary();
      }
    }
}

void
mfix::PreProjectionRedistribution (Real l_time)
{
    if ( m_redistribution_type == "StateRedist")
    {
      // We must fill internal ghost values before calling redistribution
      // We also need any physical boundary conditions imposed if we are
      //    calling state redistribution (because that calls the slope routine)

      // We really only need to fillpatch velocity but what the heck -- for now 
      //    we'll do them all
      fillpatch_all(get_vel_g(), get_ro_g(), get_h_g(),
                    get_trac(), get_X_gk(), l_time);

      for (int lev = 0; lev <= finest_level; lev++)
      {
        auto& ld = *m_leveldata[lev];

        int ncomp = AMREX_SPACEDIM;

        // Make a temporary to the redistributed velocity into
        MultiFab new_vel(grids[lev], dmap[lev], ncomp, 0);
        new_vel.setVal(0.);

        auto const& bc_vel = get_hydro_velocity_bcrec_device_ptr();

        for (MFIter mfi(*ld.vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& fact = EBFactory(lev);

            EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
            Array4<EBCellFlag const> const& flag = flagfab.const_array();

            Array4<Real> const& vel_redist = new_vel.array(mfi);
            Array4<Real> const& vel_orig   = ld.vel_g->array(mfi);

            if ( (flagfab.getType(amrex::grow(bx,1)) != FabType::covered) &&
                 (flagfab.getType(amrex::grow(bx,1)) != FabType::regular) )
            {
                Array4<Real const> fcx, fcy, fcz, ccc, vfrac, apx, apy, apz;
                fcx = fact.getFaceCent()[0]->const_array(mfi);
                fcy = fact.getFaceCent()[1]->const_array(mfi);
                fcz = fact.getFaceCent()[2]->const_array(mfi);
                ccc   = fact.getCentroid().const_array(mfi);
                apx = fact.getAreaFrac()[0]->const_array(mfi);
                apy = fact.getAreaFrac()[1]->const_array(mfi);
                apz = fact.getAreaFrac()[2]->const_array(mfi);
                vfrac = fact.getVolFrac().const_array(mfi);

                Redistribution::ApplyToInitialData( bx, ncomp, vel_redist, vel_orig,
                                                    flag, apx, apy, apz, vfrac,
                                                    fcx, fcy, fcz, ccc, bc_vel,
                                                    geom[lev],
                                                    (m_redistribution_type == "StateRedist") ? 
                                                         "NewStateRedist" : m_redistribution_type);

            } else {
                amrex::ParallelFor(bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    vel_redist(i,j,k,n) = vel_orig(i,j,k,n);
                });
            }
        }

        // Copy back into the original MultiFab
        MultiFab::Copy(*ld.vel_g, new_vel, 0,0, AMREX_SPACEDIM, 0);

        // We fill internal ghost values after calling redistribution
        ld.vel_g->FillBoundary();
      }
    }
}
