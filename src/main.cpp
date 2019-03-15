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
bool write_ls         = false;

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
std::string static_plt_file {"const_plt"};

bool plotfile_on_restart = false;

int par_ascii_int = -1;
int last_par_ascii  = -1;
std::string par_ascii_file {"par"};

int avg_int = -1;
int last_avg  = -1;
std::string avg_file {"avg_region"};

std::string mfix_dat {"mfix.dat"};

void set_ptr_to_mfix(mfix& my_mfix);

void ReadParameters ()
{
  {
     ParmParse pp("amr");

     pp.query("stop_time", stop_time);
     pp.query("max_step", max_step);

     pp.query("check_file", check_file);
     pp.query("check_int", check_int);

     pp.query("plot_file", plot_file);
     pp.query("plot_int", plot_int);

     pp.query("plotfile_on_restart", plotfile_on_restart);

     pp.query("avg_int", avg_int );
     pp.query("avg_file", avg_file);

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

     pp.query("input_deck", mfix_dat);
     pp.query("write_user", write_user);
     pp.query("write_eb_surface", write_eb_surface);
     pp.query("write_ls", write_ls);
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
    { // This start bracket and the end bracket before Finalize are essential so
      // that the mfix object is deleted before Finalize

    BL_PROFILE_VAR("main()", pmain)
    BL_PROFILE_REGION_START("mfix::main()");

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

    const char *cmfix_dat = mfix_dat.c_str();
    int name_len=mfix_dat.length();

    // Loads parameters (data) from fortran backend. Most notably this
    // subroutine loads the parameters from the `mfix.dat` file:
    //     mfix_get_data -> get_data -> read_namelist
    //                                        |
    //      (loads `mfix.dat`) ---------------+
    mfix_get_data( &solve_fluid, &solve_dem, &steady_state,
                   &dt, &dt_min, &dt_max,
                   &stop_time, &call_udf, &name_len, cmfix_dat
                  );

    // Default constructor. Note inheritance: mfix : AmrCore : AmrMesh
    //                                                             |
    //  => Geometry is constructed here: (constructs Geometry) ----+
    mfix my_mfix;

    my_mfix.get_input_bcs();

    if ( ParallelDescriptor::IOProcessor() )
      check_inputs(&dt);

    my_mfix.SetParameters(steady_state);

    // Set global static pointer to mfix object. Used by fill-patch utility
    set_ptr_to_mfix(my_mfix);

    // Initialize internals from ParamParse database
    my_mfix.InitParams(solve_fluid, solve_dem, call_udf);

    // Initialize memory for data-array internals
    my_mfix.ResizeArrays();

    // Initialize EB geometry. This needs to be done before grid creation (in
    // mfix::Init), as the grids are created using each EB-level's volfrac.
    my_mfix.make_eb_geometry();

    // Initialize derived internals
    my_mfix.Init(dt, time);

    // Create EB factories on new grids
    my_mfix.make_eb_factories();

    if (solve_dem)
    {
        // Fill level-sets on each level
        my_mfix.fill_eb_levelsets();
    }

    // Either init from scratch or from the checkpoint file
    int restart_flag = 0;
    if (restart_file.empty())
    {
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

        // NOTE: during replication 1) this also re-builds ebfactories and
        // level-set 2) this can change the grids
        IntVect Nrep(repl_x,repl_y,repl_z);
        my_mfix.Restart(restart_file, &nstep, &dt, &time, Nrep);
    }

    // This checks if we want to regrid using the KDTree or KnapSack approach
    amrex::Print() << "Regridding at step " << nstep << std::endl;
    my_mfix.Regrid();

    if (solve_dem && write_ls)
        my_mfix.WriteStaticPlotFile(static_plt_file);

    my_mfix.PostInit(dt, time, nstep, restart_flag, stop_time);

    // Write out EB sruface
    if(write_eb_surface)
      my_mfix.WriteEBSurface();

    Real end_init = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(end_init, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
       std::cout << "Time spent in init      " << end_init << std::endl;

    int finish  = 0;
    int estatus = 0;


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

    if ( avg_int > 0 )
      {
        my_mfix.WriteAverageRegions( avg_file, nstep, time );
        last_avg = nstep;
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

                my_mfix.Evolve(nstep,dt,prev_dt,time,stop_time);

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


                    if ( ( avg_int > 0) && ( nstep %  avg_int == 0 ) )
                      {
                        my_mfix.WriteAverageRegions( avg_file, nstep, time );
                        last_avg = nstep;
                      }

                }

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

    } // This end bracket and the start bracket after Initialize are essential so
      // that the mfix object is deleted before Finalize

    amrex::Finalize();
    return 0;
}
