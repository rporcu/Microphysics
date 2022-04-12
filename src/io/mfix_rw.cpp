#include <mfix_rw.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>

#include <AMReX_ParmParse.H>

using namespace amrex;


namespace MfixIO {

MfixRW::MfixRW (int nlev_in,
                amrex::Vector<amrex::BoxArray>& grids_in,
                amrex::Vector<amrex::Geometry>& geom_in,
                MFIXParticleContainer* pc_in,
                FluidPhase& fluid_in,
                amrex::Vector<std::unique_ptr<LevelData>>& m_leveldata_in,
                amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>& ebfactory_in,
                amrex::Vector<amrex::DistributionMapping>& dmap_in,
                bool ooo_debug_in,
                amrex::Vector<std::unique_ptr<amrex::MultiFab>>& level_sets_in,
                const amrex::Vector<amrex::BoxArray>& box_array_in,
                int levelset_refinement_in,
                int levelset_pad_in,
                int levelset_eb_refinement_in,
                int levelset_eb_pad_in,
                SolidsPhase& solids_in,
                Reactions& reactions_in,
                amrex::Vector<amrex::MultiFab*> particle_cost_in,
                amrex::Vector<amrex::MultiFab*> particle_proc_in,
                amrex::Vector<amrex::MultiFab*> fluid_proc_in,
                amrex::Real covered_val_in,
                const amrex::Vector<amrex::IntVect>& ref_ratio_in,
                amrex::Vector<const amrex::EB2::Level*>& eb_levels_in,
                int nghost_eb_basic_in,
                int nghost_eb_volume_in,
                int nghost_eb_full_in,
                amrex::EBSupport& m_eb_support_level_in,
                std::string load_balance_type_in,
                BCList& bc_list_in,
                Vector<std::unique_ptr<EBFArrayBoxFactory>>& particle_ebfactory_in)
  : finest_level(nlev_in-1)
  , nlev(nlev_in)
  , grids(grids_in)
  , geom(geom_in)
  , pc(pc_in)
  , fluid(fluid_in)
  , m_leveldata(m_leveldata_in)
  , ebfactory(ebfactory_in)
  , dmap(dmap_in)
  , ooo_debug(ooo_debug_in)
  , level_sets(level_sets_in)
  , box_array(box_array_in)
  , levelset_refinement(levelset_refinement_in)
  , levelset_pad(levelset_pad_in)
  , levelset_eb_refinement(levelset_eb_refinement_in)
  , levelset_eb_pad(levelset_eb_pad_in)
  , solids(solids_in)
  , reactions(reactions_in)
  , particle_cost(particle_cost_in)
  , particle_proc(particle_proc_in)
  , fluid_proc(fluid_proc_in)
  , covered_val(covered_val_in)
  , ref_ratio(ref_ratio_in)
  , eb_levels(eb_levels_in)
  , nghost_eb_basic(nghost_eb_basic_in)
  , nghost_eb_volume(nghost_eb_volume_in)
  , nghost_eb_full(nghost_eb_full_in)
  , m_eb_support_level(m_eb_support_level_in)
  , load_balance_type(load_balance_type_in)
  , bc_list(bc_list_in)
  , particle_ebfactory(particle_ebfactory_in)
{
  readParameters();
}


void MfixRW::readParameters ()
{
  {
     ParmParse pp("amr");

     // Checkpoint output control
     pp.query("checkpoint_files_output", checkpoint_files_output);
     pp.query("check_file", check_file);
     pp.query("check_int", check_int);


     // Plot output control
     pp.query("plot_file", plot_file);
     pp.query("plotfile_on_restart", plotfile_on_restart);

     pp.query("plot_int", plot_int);
     //pp.query("plot_per_exact", plot_per_exact);
     pp.query("plot_per_approx", plot_per_approx);

     if ((plot_int       > 0 && plot_per_exact  > 0) ||
         (plot_int       > 0 && plot_per_approx > 0) /*||
         (plot_per_exact > 0 && plot_per_approx > 0) */ )
       amrex::Abort("Must choose only one of plot_int or plot_per_exact or plot_per_approx");


     // Ascent output control
     pp.query("ascent_on_restart", ascent_on_restart);

     pp.query("ascent_int", ascent_int);
     pp.query("ascent_per_approx", ascent_per_approx);

     if ((ascent_int > 0 && ascent_per_approx > 0) )
       amrex::Abort("Must choose only one of ascent_int or ascent_per_approx");

     pp.query("avg_int", avg_int );
     pp.query("avg_file", avg_file);

     pp.query("par_ascii_file", par_ascii_file);
     pp.query("par_ascii_int", par_ascii_int);

     pp.query("restart", restart_file);

     pp.query("repl_x", repl_x);
     pp.query("repl_y", repl_y);
     pp.query("repl_z", repl_z);
     pp.query("regrid_int",regrid_int);

     if ( regrid_int == 0 )
       amrex::Abort("regrid_int must be > 0 or < 0");

     pp.queryarr("avg_p_g", avg_p_g);
     pp.queryarr("avg_ep_g", avg_ep_g);
     pp.queryarr("avg_vel_g", avg_vel_g);
     pp.queryarr("avg_T_g", avg_T_g);

     pp.queryarr("avg_vel_p", avg_vel_p);

     pp.queryarr("avg_T_p", avg_T_p);

     // Regions geometry
     pp.queryarr("avg_region_x_e", avg_region_x_e);
     pp.queryarr("avg_region_x_w", avg_region_x_w);
     pp.queryarr("avg_region_y_n", avg_region_y_n);
     pp.queryarr("avg_region_y_s", avg_region_y_s);
     pp.queryarr("avg_region_z_t", avg_region_z_t);
     pp.queryarr("avg_region_z_b", avg_region_z_b);
  }

  {
     ParmParse pp("mfix");

     pp.query("stop_time", stop_time);
     pp.query("max_step", max_step);

     pp.query("write_eb_surface", write_eb_surface);
     pp.query("write_ls", write_ls);
     pp.query("stop_for_unused_inputs", stop_for_unused_inputs);
  }

#ifdef MFIX_CATALYST
  {
    ParmParse pp("catalyst");
    pp.query("script", catalyst_script);
  }
#endif
}


void MfixRW::writeNow (int nstep, Real time, Real dt, bool first, bool last)
{


/*--------------------------------------------------------------------------------------------------
 *
 *                                     AMReX Plot File Output Control
 *
 *------------------------------------------------------------------------------------------------*/
    int plot_test = 0;

    if ( first ) {
        if ( (restart_file.empty() || plotfile_on_restart) &&
          (plot_int > 0 /*|| plot_per_exact > 0*/ || plot_per_approx > 0) )
          plot_test = 1;
    }

    else if ( last && plot_int > 0 ) {
        plot_test = 1;
    }

    else if (plot_per_approx > 0.0)
    {
        // Check to see if we've crossed a plot_per_approx interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = static_cast<int>( (time-dt) / plot_per_approx );
        int num_per_new = static_cast<int>( (time   ) / plot_per_approx );

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next plot_per_approx interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * amrex::Math::abs(time);
        const Real next_plot_time = (num_per_old + 1) * plot_per_approx;

        if ((num_per_new == num_per_old) && amrex::Math::abs(time - next_plot_time) <= eps)
        {
            num_per_new += 1;
        }

        // Similarly, we have to account for the case where the old time is within
        // machine epsilon of the beginning of this interval, so that we don't double
        // count that time threshold -- we already plotted at that time on the last timestep.

        if ((num_per_new != num_per_old) && amrex::Math::abs((time - dt) - next_plot_time) <= eps)
            num_per_old += 1;

        if (num_per_old != num_per_new)
            plot_test = 1;

    }/*
    else if ( plot_per_exact  > 0 && (amrex::Math::abs(remainder(time, plot_per_exact)) < 1.e-12) )
    {
        plot_test = 1;
    }*/


    if ( (plot_test == 1) || ( ( plot_int > 0) && ( nstep %  plot_int == 0 ) ) )
    {
        if (fluid.solve)
          ComputeVort();

        WritePlotFile(plot_file, nstep, time);
    }


/*--------------------------------------------------------------------------------------------------
 *
 *                                       Ascent Output Control
 *
 *------------------------------------------------------------------------------------------------*/

#ifdef AMREX_USE_ASCENT
    int ascent_test = 0;

    if ( first )
    {
        if ((restart_file.empty() || ascent_on_restart) &&
            (ascent_int > 0 || ascent_per_approx > 0) )
            ascent_test = 1;
    }
    else if (ascent_per_approx > 0.0)
    {
        // Check to see if we've crossed a ascent_per_approx interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = static_cast<int>( (time-dt) / ascent_per_approx );
        int num_per_new = static_cast<int>( (time   ) / ascent_per_approx );

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next ascent_per_approx interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * amrex::Math::abs(time);
        const Real next_ascent_time = (num_per_old + 1) * ascent_per_approx;

        if ((num_per_new == num_per_old) && amrex::Math::abs(time - next_ascent_time) <= eps)
        {
            num_per_new += 1;
        }

        // Similarly, we have to account for the case where the old time is within
        // machine epsilon of the beginning of this interval, so that we don't double
        // count that time threshold -- we already plotted at that time on the last timestep.

        if ((num_per_new != num_per_old) && amrex::Math::abs((time - dt) - next_ascent_time) <= eps)
            num_per_old += 1;

        if (num_per_old != num_per_new)
            ascent_test = 1;

    }

    if ( (ascent_test == 1) || ( ( ascent_int > 0) && ( nstep %  ascent_int == 0 ) ) )
    {
        const int myProc = ParallelDescriptor::MyProc();
        WriteAscentFile(nstep, time);
    }
#endif


/*--------------------------------------------------------------------------------------------------
 *
 *                               AMReX checkpoint file output control
 *
 *------------------------------------------------------------------------------------------------*/

    if (checkpoint_files_output && check_int > 0) {

        // We automatically write checkpoint files with the initial data
        if ( first ) {
            if ( restart_file.empty() ) {
                WriteCheckPointFile(check_file, nstep, dt, time);
                last_chk = nstep;
            }
        }
        // We automatically write checkpoint files with the final data
        else if (last) {
            if ( nstep != last_chk) {
                WriteCheckPointFile(check_file, nstep, dt, time);
                last_chk = nstep;
            }
        }
        else if ( nstep %  check_int == 0 ) {
            WriteCheckPointFile( check_file, nstep, dt, time );
            last_chk = nstep;
        }
    }

/*--------------------------------------------------------------------------------------------------
 *
 *                               AMReX particle ASCII output control
 *
 *------------------------------------------------------------------------------------------------*/
    if ( par_ascii_int > 0) {
        if ( first || last ) {
            WriteParticleAscii(par_ascii_file, nstep);
            last_par_ascii = nstep;
        }
        else if ( nstep %  par_ascii_int == 0 ) {
            WriteParticleAscii( par_ascii_file, nstep );
            last_par_ascii = nstep;
        }
    }


/*--------------------------------------------------------------------------------------------------
 *
 *                               MFIX averaging region output control
 *
 *------------------------------------------------------------------------------------------------*/

    if ( avg_int > 0 ) {
        if ( first || last ) {
            WriteAverageRegions(avg_file, nstep, time);
            last_avg = nstep;
        }
        else if ( nstep %  avg_int == 0 ) {
            WriteAverageRegions( avg_file, nstep, time );
            last_avg = nstep;
        }
    }


/*--------------------------------------------------------------------------------------------------
 *
 *                                  MFIX mass balance output control
 *
 *------------------------------------------------------------------------------------------------*/
    int mass_balance_report_test = 0;
    if (mass_balance_report_per_approx > 0.0)
      {
        // Check to see if we've crossed a mass_balance_report_per_approx interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = static_cast<int>( (time-dt) / mass_balance_report_per_approx);
        int num_per_new = static_cast<int>( (time   ) / mass_balance_report_per_approx);

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next 
        // mass_balance_report_per_approx interval at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * amrex::Math::abs(time);
        const Real next_mass_balance_report_time = (num_per_old + 1) * mass_balance_report_per_approx;

        if ((num_per_new == num_per_old) && amrex::Math::abs(time - next_mass_balance_report_time) <= eps)
        {
            num_per_new += 1;
        }

        // Similarly, we have to account for the case where the old time is within
        // machine epsilon of the beginning of this interval, so that we don't double
        // count that time threshold -- we already plotted at that time on the last timestep.

        if ((num_per_new != num_per_old) && amrex::Math::abs((time - dt) - next_mass_balance_report_time) <= eps)
            num_per_old += 1;

        if (num_per_old != num_per_new)
            mass_balance_report_test = 1;

    }

    if ( (mass_balance_report_test == 1) ||
         (( mass_balance_report_int > 0) &&
          ( nstep %  mass_balance_report_int == 0 ) ) ) {
      WriteMassBalanceReport(time);
    }


}

void MfixRW::writeEBSurface() const
{
   if(write_eb_surface)
     WriteMyEBSurface();
}


void MfixRW::writeStaticPlotFile() const
{
   if ((DEM::solve || PIC::solve) && write_ls)
      WriteStaticPlotFile(static_plt_file);
}

void MfixRW::reportGridStats() const
{
   if (fluid.solve)
     ReportGridStats();
}

} // end of namespace MfixIO
