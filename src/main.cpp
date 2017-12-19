#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include <mfix_level.H>
#include <mfix_F.H>

int   max_step    = -1;
int   verbose     = -1;
int   regrid_int  = -1;
Real stop_time = -1.0;

std::string restart_file {""};

int repl_x = 1;
int repl_y = 1;
int repl_z = 1;

int  check_int = -1;
std::string check_file {"chk"};

int   plot_int = -1;
std::string plot_file {"plt"};

bool plotfile_on_restart = false;

int par_ascii_int = -1;
std::string par_ascii_file {"par"};

void ReadParameters ()
{
  // Traditionally, max_step and stop_time do not have prefix.
  {
  ParmParse pp;
  pp.query("max_step", max_step);
  pp.query("stop_time", stop_time);
  }

  // Traditionally, these have prefix "amr", but we will
  // give them prefix mfix to make it clear that they affect the
  // behavior of the solver and not amr (even thought they are read
  // via BoxLib
  ParmParse pp("amr");

  pp.add("blocking_factor",1);

  pp.query("check_file", check_file);
  pp.query("check_int", check_int);

  pp.query("plot_file", plot_file);
  pp.query("plot_int", plot_int);

  pp.query("plotfile_on_restart", plotfile_on_restart);

  pp.query("par_ascii_file", par_ascii_file);
  pp.query("par_ascii_int", par_ascii_int);

  pp.query("restart", restart_file);
  pp.query("repl_x", repl_x);
  pp.query("repl_y", repl_y);
  pp.query("repl_z", repl_z);
  pp.query("verbose", verbose);

  pp.query("regrid_int",regrid_int);
}

int main (int argc, char* argv[])
{
    // Issue an error if AMR input file is not given
    if ( argc < 2 )
       amrex::Abort("AMReX input file missing");

    // AMReX will now read the inputs file and the command line arguments, but the
    //        command line arguments are in mfix-format so it will just ignore them.
    amrex::Initialize(argc,argv);

    // Copy arguments into MFIX -- note that the first argument is now the name of the
    //      inputs file to be read by AMReX, so we only pass the arguments after that
    for(int i=2; i < argc; i++) {
       int nlen = strlen(argv[i]);

       // If-statement avoids passing the name of the mfix input file if it is
       // specified on the command line or any AMReX command.
       if ( strstr(argv[i], "input_file") == NULL && strstr(argv[i], "amr") == NULL)
         mfix_add_argument(argv[i], &nlen);
    }

    Real strt_time = ParallelDescriptor::second();

    ReadParameters();

    int max_nit;
    int solve_fluid;
    int solve_dem;
    int steady_state;
    int call_udf;
    Real dt, dt_min, dt_max, tstop;
    Real time=0.0L;
    int nstep = 0;  // Which time step are we on
    Real normg;
    int set_normg;

    mfix_get_data( &solve_fluid,
       &solve_dem,
       &steady_state,
       &dt, &dt_min, &dt_max, &tstop, &max_nit,
       &normg, &set_normg, &call_udf);

    if ( ParallelDescriptor::IOProcessor() )
       check_inputs(&dt);

    int lev = 0;

    // Note that the constructor constructs the Geometry object now.
    mfix_level my_mfix;

    my_mfix.InitParams(solve_fluid,solve_dem,max_nit,call_udf);

    my_mfix.Init(lev,dt,time);

    // Either init from scratch or from the checkpoint file
    int restart_flag = 0;
    if (restart_file.empty())
    {
       my_mfix.InitLevelData(lev,dt,time);
    }
    else
    {
       restart_flag = 1;
       IntVect Nrep(repl_x,repl_y,repl_z);
       my_mfix.Restart( restart_file, &nstep, &dt, &time, Nrep);

       // This call checks if we want to regrid using the
       //   max_grid_size just read in from the inputs file used to restart
       //   (only relevant if load_balance_type = "FixedSize" or "KnapSack")

       // Note that this call does not depend on regrid_int
       my_mfix.RegridOnRestart(lev);
    }

    // We move this to after restart and/or regrid so we make the EB data structures with the correct 
    //    BoxArray and DistributionMapping
    bool hourglass = false;
    if (hourglass)
       my_mfix.make_eb_hourglass(lev);
    else
       my_mfix.make_eb_geometry(lev);

    // This checks if we want to regrid using the KDTree approach
    //    (only if load_balance_type = "KDTree")
    my_mfix.Regrid(lev,nstep);

    my_mfix.PostInit( lev, dt, time, nstep, restart_flag );

    // Write out EB sruface
    my_mfix.WriteEBSurface(lev);

    Real end_init = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(end_init, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
       std::cout << "Time spent in init      " << end_init << std::endl;

    int finish  = 0;
    int estatus = 0;

    // Call to output before entering time march loop
    if (solve_fluid && ParallelDescriptor::IOProcessor()  && solve_dem )
       my_mfix.output(lev,estatus,finish,nstep,dt,time);

    // Initialize prev_dt here; it will be re-defined by call to evolve_fluid but
    // only if solve_fluid = T
    Real prev_dt = dt;

    // We automatically write checkpoint and plotfiles with the initial data
    //    if plot_int > 0
    if ( (restart_file.empty() || plotfile_on_restart) && plot_int > 0 )
       my_mfix.WritePlotFile( plot_file, nstep, dt, time );

    // We automatically write checkpoint files with the initial data
    //    if check_int > 0
    if ( restart_file.empty() && check_int > 0 )
       my_mfix.WriteCheckPointFile( check_file, nstep, dt, time );

    // We automatically write ASCII files with the particle data
    //    if par_ascii_int > 0
    if ( par_ascii_int > 0 )
       my_mfix.WriteParticleAscii( par_ascii_file, nstep );

    if (time <  tstop)
    {
       while (finish == 0)
       {
          mfix_usr1();

          Real strt_step = ParallelDescriptor::second();

          if (!steady_state && regrid_int > -1 && nstep%regrid_int == 0)
             my_mfix.Regrid(lev,nstep);

          my_mfix.Evolve(lev,nstep,set_normg,dt,prev_dt,time,normg);

          Real end_step = ParallelDescriptor::second() - strt_step;
          ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
          if (ParallelDescriptor::IOProcessor())
             std::cout << "Time per step        " << end_step << std::endl;

          if (!steady_state)
          {
             time += prev_dt;
             nstep++;

             if ( ( plot_int > 0) && ( nstep %  plot_int == 0 ) )
                my_mfix.WritePlotFile( plot_file, nstep, dt, time );

             if ( ( check_int > 0) && ( nstep %  check_int == 0 ) )
                my_mfix.WriteCheckPointFile( check_file, nstep, dt, time );

             if ( ( par_ascii_int > 0) && ( nstep %  par_ascii_int == 0 ) )
                my_mfix.WriteParticleAscii( par_ascii_file, nstep );
          }

          if (ParallelDescriptor::IOProcessor() && solve_dem )
             my_mfix.output(lev,estatus,finish,nstep,dt,time);

          // Mechanism to terminate MFIX normally.
          if (steady_state || (time + 0.1*dt >= tstop)) finish = 1;
       }
    }

    // Dump plotfile at the end if enabled for steady state
    if (steady_state) {
        nstep = 1;
        if ( check_int > 0)
           my_mfix.WriteCheckPointFile( check_file    , nstep, dt, time );
        if ( plot_int > 0 )
           my_mfix.WritePlotFile      ( plot_file     , nstep, dt, time );
        if ( par_ascii_int > 0 )
           my_mfix.WriteParticleAscii ( par_ascii_file, nstep );
    }

    my_mfix.usr3(0);

    Real end_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(end_time, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
       std::cout << "Time spent in main      " << end_time << std::endl;
       std::cout << "Time spent in main-init " << end_time-end_init << std::endl;
    }

    amrex::Finalize();
    return 0;
}
