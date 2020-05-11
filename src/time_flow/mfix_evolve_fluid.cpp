#include <mfix_F.H>
#include <mfix.H>

#include <AMReX_VisMF.H>
#include <MFIX_MFHelpers.H>
#include <MFIX_DEM_Parms.H>
#include <MFIX_PIC_Parms.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

void
mfix::EvolveFluid (int nstep, Real& dt,  Real& time, Real stop_time, Real& coupling_timing)
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

    // Extrapolate boundary values for ro_g, temperature tracer, ep_g and mu_g
    // The subsequent call to mfix_set_scalar_bcs will only overwrite
    // ep_g ghost values for PINF and POUT
    for (int lev = 0; lev <= finest_level; lev++)
    {
      m_leveldata[lev]->ro_g->FillBoundary(geom[lev].periodicity());
      m_leveldata[lev]->T_g->FillBoundary(geom[lev].periodicity());
      m_leveldata[lev]->trac->FillBoundary(geom[lev].periodicity());
      m_leveldata[lev]->ep_g->FillBoundary(geom[lev].periodicity());
      m_leveldata[lev]->mu_g->FillBoundary(geom[lev].periodicity());
    }

    // Fill ghost nodes and reimpose boundary conditions
    //mfix_set_velocity_bcs(time, vel_g, 0);

    mfix_set_density_bcs(time, get_ro_g());
    //mfix_set_temperature_bcs(time, get_T_g()); // TODO do we need it?
    //mfix_set_scalar_bcs(time, trac, mu_g);

    //
    // Start loop: if we are not seeking a steady state solution,
    // the loop will execute only once
    //
    int keep_looping = 1;
    int iter = 1;

    // Create temporary multifabs to hold the old-time conv and divtau
    //    so we don't have to re-compute them in the corrector
    Vector< MultiFab* >  conv_u_old;
    Vector< MultiFab* >  conv_s_old;
    Vector< MultiFab* >  divtau_old;
    Vector< MultiFab* >    laps_old;
    Vector< MultiFab* > laptemp_old;

     conv_u_old.resize(finest_level+1);
     conv_s_old.resize(finest_level+1);
     divtau_old.resize(finest_level+1);
       laps_old.resize(finest_level+1);
    laptemp_old.resize(finest_level+1);

    for (int lev = 0; lev <= finest_level; lev++)
    {
       conv_u_old[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
       // TODO: check this, I set 3 components since we have density, tracer and
       // temperature
       conv_s_old[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
       divtau_old[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
         laps_old[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), *ebfactory[lev]);
       laptemp_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);

        conv_u_old[lev]->setVal(0.0);
        conv_s_old[lev]->setVal(0.0);
        divtau_old[lev]->setVal(0.0);
          laps_old[lev]->setVal(0.0);
       laptemp_old[lev]->setVal(0.0);
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

          MultiFab& T_g = *m_leveldata[lev]->T_g;
          MultiFab& T_go = *m_leveldata[lev]->T_go;

          MultiFab& trac = *m_leveldata[lev]->trac;
          MultiFab& trac_o = *m_leveldata[lev]->trac_o;

          MultiFab& vel_g = *m_leveldata[lev]->vel_g;
          MultiFab& vel_go = *m_leveldata[lev]->vel_go;

          // Back up field variables to old
          MultiFab::Copy(ep_go, ep_g, 0, 0, ep_g.nComp(), ep_go.nGrow());
          MultiFab::Copy(p_go, p_g, 0, 0, p_g.nComp(), p_go.nGrow());
          MultiFab::Copy(ro_go, ro_g, 0, 0, ro_g.nComp(), ro_go.nGrow());
          MultiFab::Copy(T_go, T_g, 0, 0, T_g.nComp(), T_go.nGrow());
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
        if (DEM::solve or PIC::solve) {
          Real start_drag = ParallelDescriptor::second();
          mfix_calc_drag_fluid(time);
          coupling_timing += ParallelDescriptor::second() - start_drag;
        }

        // Predictor step
        bool proj_2_pred = true;
        mfix_apply_predictor(conv_u_old, conv_s_old, divtau_old, laps_old, laptemp_old, time, dt, proj_2_pred);

        // Calculate drag coefficient
        if (DEM::solve or PIC::solve)
        {
          Real start_drag = ParallelDescriptor::second();
          amrex::Print() << "\nRecalculating drag ..." << std::endl;
          mfix_calc_drag_fluid(new_time);
          coupling_timing += ParallelDescriptor::second() - start_drag;
        }

        bool proj_2_corr = true;
        // Corrector step
        if (!steady_state)
           mfix_apply_corrector(conv_u_old, conv_s_old, divtau_old, laps_old, laptemp_old, time, dt, proj_2_corr);

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
       delete  conv_u_old[lev];
       delete  conv_s_old[lev];
       delete  divtau_old[lev];
       delete    laps_old[lev];
       delete laptemp_old[lev];
    }

    BL_PROFILE_REGION_STOP("mfix::EvolveFluid");
}
