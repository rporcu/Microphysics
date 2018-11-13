#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include <mfix.H>
#include <mfix_F.H>

int   max_step    = -1;
int   regrid_int  = -1;
Real stop_time    = -1.0;

bool write_user   = false;
bool write_eb_surface = false;

std::string restart_file {""};

int repl_x = 1;
int repl_y = 1;
int repl_z = 1;

int  check_int = -1;
int  last_chk  = -1;
std::string check_file {"chk"};

int   plot_int = -1;
int   last_plt = -1;
std::string plot_file {"plt"};

bool plotfile_on_restart = false;

int par_ascii_int = -1;
int  last_par_ascii  = -1;
std::string par_ascii_file {"par"};

void ReadParameters ()
{
  {
     ParmParse pp("amr");

     pp.query("stop_time", stop_time);
     pp.query("max_step", max_step);

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
     pp.query("regrid_int",regrid_int);
  }

  {
     ParmParse pp("mfix");

     pp.query("write_user", write_user);
     pp.query("write_eb_surface", write_eb_surface);
  }
}

int main (int argc, char* argv[])
{
    // Issue an error if AMR input file is not given
    if ( argc < 2 )
       amrex::Abort("AMReX input file missing");

    // AMReX will now read the inputs file and the command line arguments, but the
    //        command line arguments are in mfix-format so it will just ignore them.
    amrex::Initialize(argc,argv);
    BL_PROFILE_VAR("main()", pmain)
    BL_PROFILE_REGION_START("mfix::main()");

    // Warn that this is the xp branch
    amrex::Print() << "\n\nWARNING:\nMFIX was configured in Approximate Projection Mode\n\n";

    // Setting format to NATIVE rather than default of NATIVE_32
    FArrayBox::setFormat(FABio::FAB_NATIVE);

    // Copy arguments into MFIX -- note that the first argument is now the name of the
    //      inputs file to be read by AMReX, so we only pass the arguments after that
    for(int i=2; i < argc; i++) {
       int nlen = strlen(argv[i]);

       // If-statement avoids passing the name of the mfix input file if it is
       // specified on the command line or any AMReX command.
       if ( (strstr(argv[i], "input_file") == NULL) && (strstr(argv[i], "amr") == NULL)
                                                    && (strstr(argv[i], "mfix") == NULL) )
         mfix_add_argument(argv[i], &nlen);
    }

    Real strt_time = ParallelDescriptor::second();

    ReadParameters();

    int solve_fluid;
    int solve_dem;
    int steady_state;
    int call_udf;
    Real dt, dt_min, dt_max;
    Real time=0.0L;
    int nstep = 0;  // Current time step

    // Loads parameters (data) from fortran backend. Most notably this
    // subroutine loads the parameters from the `mfix.dat` file:
    //     mfix_get_data -> get_data -> read_namelist
    //                                        |
    //      (loads `mfix.dat`) ---------------+
    mfix_get_data( &solve_fluid, &solve_dem, &steady_state,
                   &dt, &dt_min, &dt_max,
                   &stop_time, &call_udf
                  );

    if ( ParallelDescriptor::IOProcessor() )
       check_inputs(&dt);

    // Default constructor. Note inheritance: mfix : AmrCore : AmrMesh
    //                                                                  |
    //  => Geometry is constructed here:  (constructs Geometry) --------+
    mfix my_mfix;

    // Initialize internals from ParamParse database
    my_mfix.InitParams(solve_fluid, solve_dem, call_udf);

    // Initialize memory for data-array internals
    // Note: MFIXParticleContainer is created here
    my_mfix.ResizeArrays();

    // Initialize derived internals
    my_mfix.Init(dt,time);

    // Either init from scratch or from the checkpoint file
    int restart_flag = 0;
    if (restart_file.empty())
    {
        // NOTE: this also builds ebfactories and level-set
        my_mfix.InitLevelData(dt,time);
    }
    else
    {
        restart_flag = 1;
        // NOTE: mfix::levelset__restart == true loading level-set from a
        // checkpoint file. However, if this is a replicating restart,
        // mfix::levelset__restart is set to false again (so that the level-sets
        // are recomputed for the replicated system).
        my_mfix.levelset__restart = true;

        // NOTE: 1) this also builds ebfactories and level-set 2) this can
        // change the grids (during replication)
        IntVect Nrep(repl_x,repl_y,repl_z);
        my_mfix.Restart(restart_file, &nstep, &dt, &time, Nrep);
    }

    // amrex::Print() << "Finished INIT Level Data " << std::endl;
    // exit(0);

    // This checks if we want to regrid using the KDTree or KnapSack approach
    amrex::Print() << "Regridding at step " << nstep << std::endl;
    my_mfix.Regrid();

    my_mfix.PostInit( dt, time, nstep, restart_flag, stop_time, steady_state );

    // Write out EB sruface
    if(write_eb_surface)
      my_mfix.WriteEBSurface();

    Real end_init = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(end_init, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
       std::cout << "Time spent in init      " << end_init << std::endl;

    int finish  = 0;
    int estatus = 0;

    // Call to output before entering time march loop
    if (solve_fluid && ParallelDescriptor::IOProcessor()  && solve_dem )
       my_mfix.output(estatus,finish,nstep,dt,time);

    // Initialize prev_dt here; it will be re-defined by call to evolve_fluid but
    // only if solve_fluid = T
    Real prev_dt = dt;

    // We automatically write checkpoint and plotfiles with the initial data
    //    if plot_int > 0
    if ( (restart_file.empty() || plotfile_on_restart) && plot_int > 0 )
    {
       if (solve_fluid)
          my_mfix.mfix_compute_vort();
       my_mfix.WritePlotFile( plot_file, nstep, dt, time );
    }

    // We automatically write checkpoint files with the initial data
    //    if check_int > 0
    if ( restart_file.empty() && check_int > 0 )
    {
       my_mfix.WriteCheckPointFile( check_file, nstep, dt, time );
       last_chk = nstep;
    }

    // We automatically write ASCII files with the particle data
    //    if par_ascii_int > 0
    if ( par_ascii_int > 0 )
    {
       my_mfix.WriteParticleAscii( par_ascii_file, nstep );
       last_par_ascii = nstep;
    }

    bool do_not_evolve =  !steady_state && ( (max_step == 0) ||
                     ( (stop_time >= 0.) && (time >  stop_time) ) || 
                     ( (stop_time <= 0.) && (max_step <= 0) ) );

    { // Start profiling solve here

        BL_PROFILE("mfix_solve");
        BL_PROFILE_REGION("mfix_solve");

        if ( !do_not_evolve) 
        {
            while (finish == 0)
            {
                mfix_usr1(&time);

                Real strt_step = ParallelDescriptor::second();

                if (!steady_state && regrid_int > -1 && nstep%regrid_int == 0)
                {
                   amrex::Print() << "Regridding at step " << nstep << std::endl;
                   my_mfix.Regrid();
                }

                my_mfix.Evolve(nstep,steady_state,dt,prev_dt,time,stop_time);

                Real end_step = ParallelDescriptor::second() - strt_step;
                ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
                if (ParallelDescriptor::IOProcessor())
                    std::cout << "Time per step        " << end_step << std::endl;

                if (!steady_state)
                {
                    time += prev_dt;
                    nstep++;

                    if ( ( plot_int > 0) && ( nstep %  plot_int == 0 ) )
                    {
                        if (solve_fluid)
                           my_mfix.mfix_compute_vort();
                        my_mfix.WritePlotFile( plot_file, nstep, dt, time );
                        last_plt = nstep;
                    }

                    if ( ( check_int > 0) && ( nstep %  check_int == 0 ) )
                    {
                        my_mfix.WriteCheckPointFile( check_file, nstep, dt, time );
                        last_chk = nstep;
                    }

                    if ( ( par_ascii_int > 0) && ( nstep %  par_ascii_int == 0 ) )
                    {
                        my_mfix.WriteParticleAscii( par_ascii_file, nstep );
                        last_par_ascii = nstep;
                    }

                    if (write_user) my_mfix.WriteUSER(dt, time);
                }

                if (ParallelDescriptor::IOProcessor() && solve_dem )
                    my_mfix.output(estatus,finish,nstep,dt,time);

                // Mechanism to terminate MFIX normally.
                do_not_evolve =  steady_state || (
                     ( (stop_time >= 0.) && (time+0.1*dt >= stop_time) ) ||  
                     ( max_step >= 0 && nstep >= max_step ) );
                if ( do_not_evolve ) finish = 1;
            }
        }
    }

    if (steady_state)
        nstep = 1;

    // Dump plotfile at the final time
    if ( check_int > 0 && nstep != last_chk)
        my_mfix.WriteCheckPointFile( check_file    , nstep, dt, time );
    if ( plot_int > 0  && nstep != last_plt)
        my_mfix.WritePlotFile      ( plot_file     , nstep, dt, time );
    if ( par_ascii_int > 0  && nstep != last_par_ascii)
        my_mfix.WriteParticleAscii ( par_ascii_file, nstep );

    my_mfix.usr3();

    Real end_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(end_time, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Time spent in main      " << end_time << std::endl;
        std::cout << "Time spent in main-init " << end_time-end_init << std::endl;
    }

    BL_PROFILE_REGION_STOP("mfix::main()");
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
    return 0;
}
