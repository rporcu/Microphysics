#include <mfix_rw.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>

#include <AMReX_ParmParse.H>

namespace MfixIO {

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

void MfixRW::writeNow (mfix& mfix, int nstep, Real time, Real dt, bool first, bool last)
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
      if (mfix.fluid.solve)
           mfix.mfix_compute_vort();
        mfix.WritePlotFile( plot_file, nstep, time );
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
        mfix.WriteAscentFile(nstep, time);
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
                mfix.WriteCheckPointFile(check_file, nstep, dt, time);
                last_chk = nstep;
            }
        }
        // We automatically write checkpoint files with the final data
        else if (last) {
            if ( nstep != last_chk) {
                mfix.WriteCheckPointFile(check_file, nstep, dt, time);
                last_chk = nstep;
            }
        }
        else if ( nstep %  check_int == 0 ) {
            mfix.WriteCheckPointFile( check_file, nstep, dt, time );
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
            mfix.WriteParticleAscii(par_ascii_file, nstep);
            last_par_ascii = nstep;
        }
        else if ( nstep %  par_ascii_int == 0 ) {
            mfix.WriteParticleAscii( par_ascii_file, nstep );
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
            mfix.WriteAverageRegions(avg_file, nstep, time);
            last_avg = nstep;
        }
        else if ( nstep %  avg_int == 0 ) {
            mfix.WriteAverageRegions( avg_file, nstep, time );
            last_avg = nstep;
        }
    }


/*--------------------------------------------------------------------------------------------------
 *
 *                                  MFIX mass balance output control
 *
 *------------------------------------------------------------------------------------------------*/
    int mass_balance_report_test = 0;
    if (mfix::mass_balance_report_per_approx > 0.0)
      {
        // Check to see if we've crossed a mfix::mass_balance_report_per_approx interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = static_cast<int>( (time-dt) / mfix::mass_balance_report_per_approx );
        int num_per_new = static_cast<int>( (time   ) / mfix::mass_balance_report_per_approx );

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next mfix::mass_balance_report_per_approx interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * amrex::Math::abs(time);
        const Real next_mass_balance_report_time = (num_per_old + 1) * mfix::mass_balance_report_per_approx;

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
         (( mfix::mass_balance_report_int > 0) &&
          ( nstep %  mfix::mass_balance_report_int == 0 ) ) ) {
      mfix.WriteMassBalanceReport(time);
    }


}

void MfixRW::writeEBSurface(mfix &mfix) const
{
   if(write_eb_surface) mfix.WriteMyEBSurface();
}


void MfixRW::writeStaticPlotFile(const mfix &mfix) const
{
   if ((DEM::solve || PIC::solve) && write_ls)
      mfix.WriteStaticPlotFile(static_plt_file);
}

void MfixRW::reportGridStats(const mfix &mfix) const
{
   if (mfix.fluid.solve) mfix.ReportGridStats();
}

} // end of namespace MfixIO
