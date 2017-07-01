#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include <mfix_level.H>
#include <mfix_F.H>

int   max_step = -1;
int   verbose  = -1;
Real stop_time = -1.0;

std::string restart_file {""};

int  check_int = -1;
std::string check_file {"chk"};

int   plot_int = -1;
std::string plot_file {"plt"};

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

  pp.query("par_ascii_file", par_ascii_file);
  pp.query("par_ascii_int", par_ascii_int);

  pp.query("restart_file", restart_file);
  pp.query("verbose", verbose);
}

int main (int argc, char* argv[])
{
    // Issue an error if AMR input file is not given
    if ( argc < 2 )
       amrex::Abort("AMReX input file missing");

    // AMReX will now read the inputs file and the command line arguments, but the
    //        command line arguments are in mfix-format so it will just ignore them.
    amrex::Initialize(argc,argv);

    BL_PROFILE("mfix_level::main()");

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
    if (restart_file.empty())
    {
       my_mfix.InitLevelData(lev,dt,time);
    } else {
       my_mfix.Restart( restart_file, &nstep, &dt, &time);
       my_mfix.InitLevelDataFromRestart( lev, dt, time );
    }

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
    if ( restart_file.empty() && plot_int > 0 )
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

          my_mfix.Evolve(lev,nstep,set_normg,dt,prev_dt,time,normg);

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
          if (steady_state || (time + 0.1*dt >= tstop) || (solve_dem && !solve_fluid)) finish = 1;
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

    if (ParallelDescriptor::IOProcessor())
       std::cout << "Time spent in main " << end_time << std::endl;

    amrex::Finalize();
    return 0;
}
