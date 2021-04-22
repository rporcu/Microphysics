#include <mfix.H>
#include <Redistribution.H>

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

        for (MFIter mfi(*ld.ro_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const& fact = EBFactory(lev);

            EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
            Array4<EBCellFlag const> const& flag = flagfab.const_array();

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

                int ncomp = AMREX_SPACEDIM;
                int icomp = 0;
                redistribution::redistribute_data( bx,ncomp, icomp,
                                          ld.vel_g->array(mfi), ld.vel_go->array(mfi),
                                          flag, apx, apy, apz, vfrac, fcx, fcy, fcz,
                                          ccc,geom[lev],m_redistribution_type);

                if (advect_density) {

                  ncomp = 1;
                  redistribution::redistribute_data( bx,ncomp, icomp,
                                           ld.ro_g->array(mfi), ld.ro_go->array(mfi),
                                           flag, apx, apy, apz, vfrac, fcx, fcy, fcz,
                                           ccc,geom[lev],m_redistribution_type);
                }
                if (advect_enthalpy) {

                  ncomp = 1;
                  redistribution::redistribute_data( bx,ncomp, icomp,
                                           ld.h_g->array(mfi), ld.h_go->array(mfi),
                                           flag, apx, apy, apz, vfrac, fcx, fcy, fcz,
                                           ccc,geom[lev],m_redistribution_type);
                }
                if (advect_tracer) {

                  ncomp = ntrac;
                  redistribution::redistribute_data( bx,ncomp, icomp,
                                           ld.trac->array(mfi), ld.trac_o->array(mfi),
                                           flag, apx, apy, apz, vfrac, fcx, fcy, fcz,
                                           ccc,geom[lev],m_redistribution_type);
                }
                if (advect_fluid_species) {

                  ncomp = fluid.nspecies;
                  redistribution::redistribute_data( bx,ncomp, icomp,
                                           ld.X_gk->array(mfi), ld.X_gko->array(mfi),
                                           flag, apx, apy, apz, vfrac, fcx, fcy, fcz,
                                           ccc,geom[lev],m_redistribution_type);
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

        if(advect_fluid_species)
          ld.X_gk->FillBoundary();
    }
  }
}
