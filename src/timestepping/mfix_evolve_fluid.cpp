#include <mfix.H>

#include <mfix_mf_helpers.H>
#include <mfix_eb_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_utils.H>
#include <mfix_monitors.H>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

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

    // Fill ghost nodes and reimpose boundary conditions
    //mfix_set_velocity_bcs(time, vel_g, 0);

    mfix_set_density_bcs(time, get_ro_g());

    // TODO: commenting the following makes BENCH03 GPU to pass
    //mfix_set_tracer_bcs(time, get_trac());

    if (fluid.solve_enthalpy) {
      mfix_set_temperature_bcs(time, get_T_g());
      mfix_set_enthalpy_bcs(time, get_h_g());

      if (EB::fix_temperature) {
        mfix_set_eb_temperature_bcs(get_T_g_on_eb());
      }
    }

    if (fluid.solve_species)
      mfix_set_species_bcs(time, get_X_gk());

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
    Vector< MultiFab* > lap_trac_old(finest_level+1);
    Vector< MultiFab* > ro_RHS_old(finest_level+1);
    Vector< MultiFab* > enthalpy_RHS_old(finest_level+1);
    Vector< MultiFab* > species_RHS_old(finest_level+1);
    Vector< MultiFab* > vel_RHS_old(finest_level+1);
    Vector< MultiFab* > div_J_old(finest_level+1);
    Vector< MultiFab* > div_hJ_old(finest_level+1);
    Vector< Real > rhs_pressure_g_old(finest_level+1, 0.);
    Vector< Real > rhs_pressure_g(finest_level+1, 0.);
    Vector< MultiFab* > eb_flow_vel(finest_level+1, nullptr);
    Vector< MultiFab* > eb_flow_scalars(finest_level+1, nullptr);
    Vector< MultiFab* > eb_flow_species(finest_level+1, nullptr);

    for (int lev = 0; lev <= finest_level; lev++)
    {
       conv_u_old[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
       conv_u_old[lev]->setVal(0.0);

       // 2+ntrac components since we have density, tracers and enthalpy
       conv_s_old[lev] = new MultiFab(grids[lev], dmap[lev], 2+ntrac, 0, MFInfo(), *ebfactory[lev]);
       conv_s_old[lev]->setVal(0.0);

       if (fluid.solve_enthalpy) {
         lap_T_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
         lap_T_old[lev]->setVal(0.0);
         enthalpy_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
         enthalpy_RHS_old[lev]->setVal(0.0);
       }

       if (fluid.solve_tracer) {
         lap_trac_old[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), *ebfactory[lev]);
         lap_trac_old[lev]->setVal(0.0);
       }

       if (fluid.solve_density) {
         ro_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
         ro_RHS_old[lev]->setVal(0.0);
       }

       if (fluid.solve_species) {
         conv_X_old[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies, 0, MFInfo(), *ebfactory[lev]);
         conv_X_old[lev]->setVal(0.0);

         species_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies, 0, MFInfo(), *ebfactory[lev]);
         species_RHS_old[lev]->setVal(0.0);
       }

       if (fluid.solve_species) {
         div_J_old[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies, 0, MFInfo(), *ebfactory[lev]);
         div_J_old[lev]->setVal(0.0);
       }

       if (fluid.solve_enthalpy && fluid.solve_species) {
         div_hJ_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
         div_hJ_old[lev]->setVal(0.0);
       }

       if (reactions.solve) {
         vel_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
         vel_RHS_old[lev]->setVal(0.0);
       }

       if (EB::has_flow) {

         eb_flow_vel[lev] = new MultiFab(grids[lev], dmap[lev], 3, nghost_state(), MFInfo(), *ebfactory[lev]);
         eb_flow_vel[lev]->setVal(0.0);

         eb_flow_scalars[lev] = new MultiFab(grids[lev], dmap[lev], 2+ntrac, nghost_state(), MFInfo(), *ebfactory[lev]);
         eb_flow_scalars[lev]->setVal(0.0);

         if (fluid.solve_species) {
           eb_flow_species[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies, nghost_state(), MFInfo(), *ebfactory[lev]);
           eb_flow_species[lev]->setVal(0.0);
         }
       }
    }

    // We do this call after the bc values have been set
    if (EB::has_flow) {
      mfix_set_eb_velocity_bcs(time, eb_flow_vel);
      mfix_set_eb_scalar_bcs(eb_flow_scalars, eb_flow_species);
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

        //
        // Time integration step
        //
        Real new_time = time+dt;

        // NOTE: it is important to do fillpatch in HERE (before swapping
        // pointers), and on the NEW variables (not on the OLD ones)
        fillpatch_all(get_vel_g(), get_ro_g(), get_h_g(), get_trac(), get_X_gk(), new_time);

        for (int lev = 0; lev <= finest_level; lev++)
        {
          // Swap old and new variables
          std::swap(m_leveldata[lev]->p_g, m_leveldata[lev]->p_go);
          std::swap(m_leveldata[lev]->ro_g, m_leveldata[lev]->ro_go);
          std::swap(m_leveldata[lev]->trac, m_leveldata[lev]->trac_o);
          std::swap(m_leveldata[lev]->vel_g, m_leveldata[lev]->vel_go);

          if (fluid.solve_enthalpy) {
            std::swap(m_leveldata[lev]->T_g, m_leveldata[lev]->T_go);
            std::swap(m_leveldata[lev]->h_g, m_leveldata[lev]->h_go);
          }

          if (fluid.solve_species)
            std::swap(m_leveldata[lev]->X_gk, m_leveldata[lev]->X_gko);

          if (reactions.solve)
            std::swap(m_leveldata[lev]->thermodynamic_p_g, m_leveldata[lev]->thermodynamic_p_go);

          // User hooks
          for (MFIter mfi(*m_leveldata[lev]->ep_g, false); mfi.isValid(); ++mfi)
             mfix_usr2();
        }

        // Calculate drag coefficient
        if (DEM::solve || PIC::solve) {
          Real start_drag = ParallelDescriptor::second();
          mfix_calc_txfr_fluid(get_txfr(), get_chem_txfr(), get_ep_g(),
                               get_ro_g_old(), get_vel_g_old(), get_T_g_old(),
                               get_X_gk_old(), get_thermodynamic_p_g_old(), time);

          coupling_timing += ParallelDescriptor::second() - start_drag;
        }

        // Predictor step
        bool proj_2_pred = true;
        mfix_apply_predictor(conv_u_old, conv_s_old, conv_X_old, ro_RHS_old,
            lap_trac_old, lap_T_old, enthalpy_RHS_old, species_RHS_old,
            vel_RHS_old, div_J_old, div_hJ_old, rhs_pressure_g_old, rhs_pressure_g,
            eb_flow_vel, eb_flow_scalars, eb_flow_species,
            time, dt, prev_dt, proj_2_pred, coupling_timing);

        // Corrector step
        if (advection_type() == AdvectionType::MOL && !m_steady_state) {

          // Calculate drag coefficient
          if (DEM::solve || PIC::solve) {

            Real start_drag = ParallelDescriptor::second();
            amrex::Print() << "\nRecalculating drag ..." << std::endl;
            mfix_calc_txfr_fluid(get_txfr(), get_chem_txfr(), get_ep_g(),
                                 get_ro_g(), get_vel_g(), get_T_g(), get_X_gk(),
                                 get_thermodynamic_p_g(), new_time);

            coupling_timing += ParallelDescriptor::second() - start_drag;
          }

          //  Do fillpatch in here
          fillpatch_all(get_vel_g(), get_ro_g(), get_h_g(), get_trac(), get_X_gk(), new_time);

          bool proj_2_corr = true;

          mfix_apply_corrector(conv_u_old, conv_s_old, conv_X_old, ro_RHS_old,
              lap_trac_old, lap_T_old, enthalpy_RHS_old, species_RHS_old,
              vel_RHS_old, div_J_old, div_hJ_old, rhs_pressure_g_old,
              rhs_pressure_g, eb_flow_vel, eb_flow_scalars, eb_flow_species,
              time, dt, prev_dt, proj_2_corr, coupling_timing);
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
    while (keep_looping);

    if (test_tracer_conservation)
    {
       amrex::Print() << "Sum tracer volume wgt = "
                      << Utils::volWgtSum(0, *m_leveldata[0]->trac, 0, ebfactory)
                      << " " << Utils::volEpsWgtSum(0, *m_leveldata[0]->trac, *m_leveldata[0]->ep_g, 0, ebfactory)
                      << std::endl;

//      const Box& domain = geom[0].Domain();
//
//      const auto& trac = *m_leveldata[0]->trac;
//      const auto& epg = *m_leveldata[0]->ep_g;
//
//      auto monitor = EulerianMonitor::VolumeIntegral(domain, EBFactory(0));
//
//      amrex::Print() << "Sum tracer volume wgt = "
//                     << monitor.volume_weighted_sum(trac)
//                     << " " << monitor.volume_weighted_sum(trac, epg)
//                     << std::endl;
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

       if (fluid.solve_density)
         delete ro_RHS_old[lev];

       if (fluid.solve_tracer)
         delete lap_trac_old[lev];

       if (fluid.solve_enthalpy) {
       delete enthalpy_RHS_old[lev];
       delete lap_T_old[lev];
       }

       if (fluid.solve_enthalpy && fluid.solve_species) {
         delete div_hJ_old[lev];
       }

       if (fluid.solve_species) {
         delete conv_X_old[lev];
         delete species_RHS_old[lev];
         delete div_J_old[lev];
       }

       if (reactions.solve)
         delete vel_RHS_old[lev];

       if (EB::has_flow) {
         delete eb_flow_vel[lev];
         delete eb_flow_scalars[lev];
         if (fluid.solve_species) {
           delete eb_flow_species[lev];
         }
      }

    }

    BL_PROFILE_REGION_STOP("mfix::EvolveFluid");
}
