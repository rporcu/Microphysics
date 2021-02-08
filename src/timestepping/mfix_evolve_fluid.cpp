#include <mfix.H>

#include <AMReX_VisMF.H>
#include <mfix_mf_helpers.H>
#include <mfix_eb_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_des_heterogeneous_rates_K.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

void
mfix::EvolveFluid (int nstep,
                   Real& dt,
                   Real& prev_dt,
                   Real& time,
                   Real stop_time,
                   Real& coupling_timing)
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
      m_leveldata[lev]->trac->FillBoundary(geom[lev].periodicity());
      m_leveldata[lev]->ep_g->FillBoundary(geom[lev].periodicity());
      m_leveldata[lev]->mu_g->FillBoundary(geom[lev].periodicity());
      m_leveldata[lev]->MW_g->FillBoundary(geom[lev].periodicity());

      if (advect_enthalpy) {
        m_leveldata[lev]->T_g->FillBoundary(geom[lev].periodicity());
        m_leveldata[lev]->cp_g->FillBoundary(geom[lev].periodicity());
        m_leveldata[lev]->k_g->FillBoundary(geom[lev].periodicity());
        m_leveldata[lev]->h_g->FillBoundary(geom[lev].periodicity());
      }

      if (advect_fluid_species) {
        m_leveldata[lev]->D_gk->FillBoundary(geom[lev].periodicity());
        m_leveldata[lev]->X_gk->FillBoundary(geom[lev].periodicity());
      }

      if (advect_enthalpy and advect_fluid_species) {
        m_leveldata[lev]->cp_gk->FillBoundary(geom[lev].periodicity());
        m_leveldata[lev]->h_gk->FillBoundary(geom[lev].periodicity());
      }
    }

    // Fill ghost nodes and reimpose boundary conditions
    //mfix_set_velocity_bcs(time, vel_g, 0);

    mfix_set_density_bcs(time, get_ro_g());

    // TODO: commenting the following makes BENCH03 GPU to pass
    //mfix_set_scalar_bcs(time, get_mu_g(), get_cp_g(), get_k_g(), get_MW_g());
    //mfix_set_tracer_bcs(time, get_trac());

    if (advect_enthalpy) {
      mfix_set_temperature_bcs(time, get_T_g());
      mfix_set_enthalpy_bcs(time, get_h_g());

      if (EB::fix_temperature) {
        const Vector< MultiFab* >& T_g_on_eb = get_T_g_on_eb();
        const Vector< MultiFab* >& k_g = get_k_g();
        const Vector< MultiFab* >& k_g_on_eb = get_k_g_on_eb();

        for (int lev(0); lev <= finest_level; ++lev) {
          T_g_on_eb[lev]->setVal(0);
          MultiFab::Copy(*k_g_on_eb[lev], *k_g[lev], 0, 0, 1, k_g_on_eb[lev]->nGrow());
        }

        mfix_set_eb_temperature_bcs(get_T_g_on_eb(), get_k_g_on_eb());
      }
    }

    if (advect_fluid_species)
      mfix_set_species_bcs(time, get_X_gk(), get_D_gk(), get_cp_gk(), get_h_gk());

    //
    // Start loop: if we are not seeking a steady state solution,
    // the loop will execute only once
    //
    int keep_looping = 1;
    int iter = 1;

    // Create temporary multifabs to hold the old-time conv and vel_RHS
    //    so we don't have to re-compute them in the corrector
    Vector< MultiFab* > conv_s_old(finest_level+1);
    Vector< MultiFab* > conv_X_old(finest_level+1);
    Vector< MultiFab* > conv_u_old(finest_level+1);
    Vector< MultiFab* > lap_T_old(finest_level+1);
    Vector< MultiFab* > lap_T_star(finest_level+1);
    Vector< MultiFab* > lap_X_old(finest_level+1);
    Vector< MultiFab* > lap_X_star(finest_level+1);
    Vector< MultiFab* > lap_trac_old(finest_level+1);
    Vector< MultiFab* > divtau_old(finest_level+1);
    Vector< MultiFab* > ro_RHS_old(finest_level+1);
    Vector< MultiFab* > enthalpy_RHS_old(finest_level+1);
    Vector< MultiFab* > species_RHS_old(finest_level+1);

    for (int lev = 0; lev <= finest_level; lev++)
    {
       // 3 components since we have density, tracer and enthalpy
       conv_s_old[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
       conv_u_old[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
       lap_T_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
       lap_T_star[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
       lap_trac_old[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), *ebfactory[lev]);
       divtau_old[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
       ro_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
       enthalpy_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);

       conv_s_old[lev]->setVal(0.0);
       conv_u_old[lev]->setVal(0.0);
       lap_T_old[lev]->setVal(0.0);
       lap_T_star[lev]->setVal(0.0);
       lap_trac_old[lev]->setVal(0.0);
       divtau_old[lev]->setVal(0.0);
       ro_RHS_old[lev]->setVal(0.0);
       enthalpy_RHS_old[lev]->setVal(0.0);

       if (advect_fluid_species) {
         conv_X_old[lev] = new MultiFab(grids[lev], dmap[lev], FLUID::nspecies, 0, MFInfo(), *ebfactory[lev]);
         lap_X_old[lev] = new MultiFab(grids[lev], dmap[lev], FLUID::nspecies, 0, MFInfo(), *ebfactory[lev]);
         lap_X_star[lev] = new MultiFab(grids[lev], dmap[lev], FLUID::nspecies, 0, MFInfo(), *ebfactory[lev]);
         species_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], FLUID::nspecies, 0, MFInfo(), *ebfactory[lev]);

         conv_X_old[lev]->setVal(0.0);
         lap_X_old[lev]->setVal(0.0);
         lap_X_star[lev]->setVal(0.0);
         species_RHS_old[lev]->setVal(0.0);
       }
    }

    do
    {
        mfix_compute_dt(nstep, time, stop_time, dt, prev_dt);

        // Set new and old time to correctly use in fillpatching
        for (int lev = 0; lev <= finest_level; lev++)
        {
            t_old[lev] = time;
            t_new[lev] = time+dt;
        }

        if (m_steady_state)
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

          MultiFab& X_gk = *m_leveldata[lev]->X_gk;
          MultiFab& X_gko = *m_leveldata[lev]->X_gko;

          MultiFab& T_g = *m_leveldata[lev]->T_g;
          MultiFab& T_go = *m_leveldata[lev]->T_go;

          MultiFab& h_g = *m_leveldata[lev]->h_g;
          MultiFab& h_go = *m_leveldata[lev]->h_go;

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

          if (advect_enthalpy) {
            MultiFab::Copy(T_go, T_g, 0, 0, T_g.nComp(), T_go.nGrow());
            MultiFab::Copy(h_go, h_g, 0, 0, h_g.nComp(), h_go.nGrow());
          }

          if (advect_fluid_species)
            MultiFab::Copy(X_gko, X_gk, 0, 0, X_gk.nComp(), X_gko.nGrow());

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
          mfix_calc_txfr_fluid(time);

          if (REACTIONS::solve) {
            mfix_calc_chem_txfr(time, get_ep_g(), get_ro_g_old(), get_X_gk_old());
          }

          coupling_timing += ParallelDescriptor::second() - start_drag;
        }

        // Predictor step
        bool proj_2_pred = true;
        mfix_apply_predictor(conv_u_old, conv_s_old, conv_X_old, ro_RHS_old,
            divtau_old, lap_trac_old, lap_T_old, lap_T_star, enthalpy_RHS_old, species_RHS_old,
            lap_X_old, lap_X_star, time, dt, prev_dt, proj_2_pred);

        // Calculate drag coefficient
        if (DEM::solve or PIC::solve)
        {
          Real start_drag = ParallelDescriptor::second();
          amrex::Print() << "\nRecalculating drag ..." << std::endl;
          mfix_calc_txfr_fluid(new_time);

          // TODO NOW: not sure we need to do this again
          if (REACTIONS::solve) {
            mfix_calc_chem_txfr(time, get_ep_g(), get_ro_g(), get_X_gk());
          }

          coupling_timing += ParallelDescriptor::second() - start_drag;
        }

        bool proj_2_corr = true;
        // Corrector step
        if (advection_type() == AdvectionType::MOL && !m_steady_state) {
           mfix_apply_corrector(conv_u_old, conv_s_old, conv_X_old, ro_RHS_old,
               divtau_old, lap_trac_old, lap_T_old, lap_T_star, enthalpy_RHS_old, species_RHS_old,
               lap_X_old, lap_X_star, time, dt, prev_dt, proj_2_corr);
        }

        //
        // Check whether to exit the loop or not
        //
        if (m_steady_state) {
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
       delete ro_RHS_old[lev];
       delete divtau_old[lev];
       delete lap_trac_old[lev];
       delete enthalpy_RHS_old[lev];
       delete lap_T_old[lev];
       delete lap_T_star[lev];

       if (advect_fluid_species) {
         delete conv_X_old[lev];
         delete lap_X_old[lev];
         delete lap_X_star[lev];
         delete species_RHS_old[lev];
       }
    }

    BL_PROFILE_REGION_STOP("mfix::EvolveFluid");
}
