#include <mfix_rw.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>

#include <AMReX_ParmParse.H>

namespace MfixIO {

void MfixRW::readParameters ()
{
  {
     ParmParse pp("amr");

     pp.query("checkpoint_files_output", checkpoint_files_output);
     pp.query("check_file", check_file);
     pp.query("check_int", check_int);

     pp.query("plot_file", plot_file);

     pp.query("plotfile_on_restart", plotfile_on_restart);
     pp.query("ascent_on_restart", ascent_on_restart);

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

     pp.query("input_deck", mfix_dat);
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

void MfixRW::writeNow (int nstep, Real time, Real dt, mfix& mfix)
{
    int plot_test = 0;
    if (mfix::plot_per_approx > 0.0)
    {
        // Check to see if we've crossed a mfix::plot_per_approx interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = static_cast<int>( (time-dt) / mfix::plot_per_approx );
        int num_per_new = static_cast<int>( (time   ) / mfix::plot_per_approx );

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next mfix::plot_per_approx interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * amrex::Math::abs(time);
        const Real next_plot_time = (num_per_old + 1) * mfix::plot_per_approx;

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

    }
    else if ( mfix::plot_per_exact  > 0 && (amrex::Math::abs(remainder(time, mfix::plot_per_exact)) < 1.e-12) )
    {
        plot_test = 1;
    }

    if ( (plot_test == 1) || ( ( mfix::plot_int > 0) && ( nstep %  mfix::plot_int == 0 ) ) )
    {
      if (mfix.fluid.solve)
           mfix.mfix_compute_vort();
        mfix.WritePlotFile( plot_file, nstep, time );
    }
#ifdef AMREX_USE_ASCENT
    int ascent_test = 0;
    if (mfix::ascent_per_approx > 0.0)
    {
        // Check to see if we've crossed a mfix::ascent_per_approx interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = static_cast<int>( (time-dt) / mfix::ascent_per_approx );
        int num_per_new = static_cast<int>( (time   ) / mfix::ascent_per_approx );

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next mfix::ascent_per_approx interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * amrex::Math::abs(time);
        const Real next_ascent_time = (num_per_old + 1) * mfix::ascent_per_approx;

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

    if ( (ascent_test == 1) || ( ( mfix::ascent_int > 0) && ( nstep %  mfix::ascent_int == 0 ) ) )
    {
        mfix.WriteAscentFile();
    }
#endif

    if ( ( check_int > 0) && ( nstep %  check_int == 0 ) )
    {
        if (checkpoint_files_output) mfix.WriteCheckPointFile( check_file, nstep, dt, time );
        last_chk = nstep;
    }

    if ( ( par_ascii_int > 0) && ( nstep %  par_ascii_int == 0 ) )
    {
        mfix.WriteParticleAscii( par_ascii_file, nstep );
        last_par_ascii = nstep;
    }


    if ( ( avg_int > 0) && ( nstep %  avg_int == 0 ) )
    {
      mfix.WriteAverageRegions( avg_file, nstep, time );
      last_avg = nstep;
    }

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

void MfixRW::writePlotFileInitial(int nstep, Real time, mfix &mfix)
{
   // Write checkpoint and plotfiles with the initial data
   if ( (restart_file.empty() || plotfile_on_restart) &&
         (mfix::plot_int > 0 || mfix::plot_per_exact > 0 || mfix::plot_per_approx > 0) )
   {
      if (mfix.fluid.solve)
         mfix.mfix_compute_vort();
      mfix.WritePlotFile(plot_file, nstep, time);
   }
}

void MfixRW::writePlotFileFinal(int nstep, Real time, mfix &mfix)
{
   if ( mfix::plot_int > 0)
      mfix.WritePlotFile(plot_file, nstep, time);
}

void MfixRW::writeCheckPointFile(int nstep, Real dt, Real time, 
      mfix &mfix)
{
   // We automatically write checkpoint files with the initial data
   //    if check_int > 0
   if ( restart_file.empty() && check_int > 0 )
   {
      if (checkpoint_files_output) mfix.WriteCheckPointFile(check_file, nstep, dt, time);
      last_chk = nstep;
   }
}

void MfixRW::writeCheckPointFileFinal(int nstep, Real dt, Real time, 
      mfix &mfix)
{
   if ( check_int > 0 && nstep != last_chk)
      if (checkpoint_files_output) mfix.WriteCheckPointFile(check_file, nstep, dt, time);
}

void MfixRW::writeParticleAscii(int nstep, const mfix &mfix)
{
   // We automatically write ASCII files with the particle data
   //    if par_ascii_int > 0
   if ( par_ascii_int > 0 )
   {
      mfix.WriteParticleAscii(par_ascii_file, nstep);
      last_par_ascii = nstep;
   }
}

void MfixRW::writeParticleAsciiFinal(int nstep, const mfix &mfix)
{
   if ( par_ascii_int > 0  && nstep != last_par_ascii)
      mfix.WriteParticleAscii(par_ascii_file, nstep);
}

void MfixRW::writeStaticPlotFile(const mfix &mfix) const
{
   if ((DEM::solve || PIC::solve) && write_ls)
      mfix.WriteStaticPlotFile(static_plt_file);
}

void MfixRW::writeAverageRegions(int nstep, Real time, const mfix &mfix)
{
   if ( avg_int > 0 )
   {
      mfix.WriteAverageRegions(avg_file, nstep, time);
      last_avg = nstep;
   }

}

void MfixRW::writeAsentFile(mfix &mfix) const
{
   amrex::ignore_unused(mfix);
#ifdef AMREX_USE_ASCENT
   if ( (restart_file.empty() || ascent_on_restart) &&
         (mfix::ascent_int > 0 || mfix::ascent_per_approx > 0) ) {
      mfix.WriteAscentFile();
   }
#endif
}

void MfixRW::reportGridStats(const mfix &mfix) const
{
   if (mfix.fluid.solve) mfix.ReportGridStats();
}

} // end of namespace MfixIO
