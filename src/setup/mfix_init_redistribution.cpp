#include <mfix.H>
#include <hydro_redistribution.H>

#include <mfix_fluid_parms.H>
void
mfix::InitialRedistribution (Real l_time)
{
    // Next we must redistribute the initial solution if we are going to use
    // StateRedist redistribution scheme
    if ( m_redistribution_type == "StateRedist")
    {
      // We use the "old" data as the input here
      // We must fill internal ghost values before calling redistribution
      // We also need any physical boundary conditions imposed if we are
      //    calling state redistribution (because that calls the slope routine)

      fillpatch_all(get_vel_g_old(), get_ro_g_old(), get_h_g_old(),
                    get_trac_old(), get_X_gk_old(), l_time);


      for (int lev = 0; lev <= finest_level; lev++)
      {
        auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*ld.ro_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const& fact = EBFactory(lev);

            EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
            Array4<EBCellFlag const> const& flag = flagfab.const_array();

            if ( (flagfab.getType(bx)                != FabType::covered) &&
                 (flagfab.getType(amrex::grow(bx,4)) != FabType::regular) )
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

                int ncomp = AMREX_SPACEDIM;

                auto const& bc_vel = get_hydro_velocity_bcrec_device_ptr();
                Redistribution::ApplyToInitialData( bx,ncomp,
                                          ld.vel_g->array(mfi), ld.vel_go->array(mfi),
                                          flag, apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                          bc_vel, geom[lev],
                                          (m_redistribution_type == "StateRedist") ? 
                                                "NewStateRedist" : m_redistribution_type);

                if (advect_density) {

                  ncomp = 1;

                  auto const& bc_den = get_density_bcrec_device_ptr();
                  Redistribution::ApplyToInitialData( bx,ncomp,
                                           ld.ro_g->array(mfi), ld.ro_go->array(mfi),
                                           flag, apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                           bc_den, geom[lev],
                                           (m_redistribution_type == "StateRedist") ? 
                                                "NewStateRedist" : m_redistribution_type);
                }
                if (advect_enthalpy) {

                  ncomp = 1;

                  auto const& bc_h = get_enthalpy_bcrec_device_ptr();
                  Redistribution::ApplyToInitialData( bx,ncomp,
                                           ld.h_g->array(mfi), ld.h_go->array(mfi),
                                           flag, apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                           bc_h, geom[lev],
                                           (m_redistribution_type == "StateRedist") ? 
                                                "NewStateRedist" : m_redistribution_type);
                }
                if (advect_tracer) {

                  ncomp = ntrac;

                  auto const& bc_t = get_tracer_bcrec_device_ptr();
                  Redistribution::ApplyToInitialData( bx,ncomp,
                                           ld.trac->array(mfi), ld.trac_o->array(mfi),
                                           flag, apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                           bc_t, geom[lev],
                                           (m_redistribution_type == "StateRedist") ? 
                                                "NewStateRedist" : m_redistribution_type);
                }
                if (solve_species) {

                  ncomp = fluid.nspecies;

                  auto const& bc_X = get_species_bcrec_device_ptr();
                  Redistribution::ApplyToInitialData( bx,ncomp,
                                           ld.X_gk->array(mfi), ld.X_gko->array(mfi),
                                           flag, apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                           bc_X, geom[lev],
                                           (m_redistribution_type == "StateRedist") ? 
                                                "NewStateRedist" : m_redistribution_type);
                }

            }
        }

        // We fill internal ghost values after calling redistribution
        ld.vel_g->FillBoundary();

        if(advect_density)
          ld.ro_g->FillBoundary();
        if(advect_enthalpy)
          ld.h_g->FillBoundary();

        if(ntrac > 0)
          ld.trac->FillBoundary();

        if(solve_species)
          ld.X_gk->FillBoundary();
    }
  }
}
