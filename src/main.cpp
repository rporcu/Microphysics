#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include <AMReX_buildInfo.H>

#include <mfix.H>

#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_fluid_parms.H>

int  max_step   = -1;
int  regrid_int = -1;
Real stop_time  = -1.0;

bool write_eb_surface = false;
bool write_ls         = false;

std::string restart_file {""};

int repl_x = 1;
int repl_y = 1;
int repl_z = 1;

int check_int = -1;
int last_chk  = -1;
std::string check_file {"chk"};

std::string plot_file {"plt"};
std::string static_plt_file {"plt_ls"};

bool plotfile_on_restart = false;
bool ascent_on_restart = false;

int par_ascii_int = -1;
int last_par_ascii  = -1;
std::string par_ascii_file {"par"};

int avg_int = -1;
int last_avg = -1;
std::string avg_file {"avg_region"};

std::string mfix_dat {"mfix.dat"};

// Set the extend domain flag by default, since the mfix default
// is different (true) from the amrex default (false)
// only if its not already specified in the inputs file
void add_par () {
   ParmParse pp("eb2");
   if(!pp.contains("extend_domain_face")) {
      pp.add("extend_domain_face",true);
   }
} 

void writeBuildInfo ();

void ReadParameters ()
{
  {
     ParmParse pp("amr");

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
  }

}

void writeNow (int nstep, Real time, Real dt, mfix& mfix)
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
        mfix.WriteCheckPointFile( check_file, nstep, dt, time );
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
}

int main (int argc, char* argv[])
{


    // check to see if it contains --describe
    if (argc >= 2) {
        for (auto i = 1; i < argc; i++) {
            if (std::string(argv[i]) == "--describe") {
                writeBuildInfo();
                return 0;
            }
        }
    }

    // Issue an error if AMR input file is not given
    if ( argc < 2 ) {
       std::cerr << "AMReX input file missing" << std::endl << std::endl;
       std::cerr << "Usage:  " << argv[0] << " inputs [--describe]" << std::endl;
       return -1;
    }

    // AMReX will now read the inputs file and the command line arguments, but the
    //        command line arguments are in mfix-format so it will just ignore them.
    amrex::Initialize(argc,argv,true,MPI_COMM_WORLD,add_par);
    { // This start bracket and the end bracket before Finalize are essential so
      // that the mfix object is deleted before Finalize
    BL_PROFILE_VAR("main()", pmain)
    BL_PROFILE_REGION_START("mfix::main()");

    // Write out the MFIX git hash (the AMReX git hash is already written)
    const char* githash_mfix = buildInfoGetGitHash(1);
    amrex::Print() << "MFiX git hash: " << githash_mfix<< "\n";

    // Setting format to NATIVE rather than default of NATIVE_32
    FArrayBox::setFormat(FABio::FAB_NATIVE);

    Real strt_time = ParallelDescriptor::second();

    Real time=0.0L;
    int nstep = 0;  // Current time step

    Real dt = -1.;

    // Default constructor. Note inheritance: mfix : AmrCore : AmrMesh
    //                                                             |
    //  => Geometry is constructed here: (constructs Geometry) ----+
    mfix mfix;

    ReadParameters();

    // Initialize internals from ParamParse database
    mfix.InitParams();

    // Initialize memory for data-array internals
    mfix.ResizeArrays();

    // Initialize EB geometry. This needs to be done before grid creation (in
    // mfix::Init), as the grids are created using each EB-level's volfrac.
    mfix.make_eb_geometry();

    // Initialize derived internals
    mfix.Init(time);

    // Create EB factories on new grids
    mfix.make_eb_factories();

    // Write out EB sruface
    if(write_eb_surface)
      mfix.WriteMyEBSurface();

    if (DEM::solve || PIC::solve)
    {
        // Fill level-sets on each level
        mfix.fill_eb_levelsets();
    }

    // Either init from scratch or from the checkpoint file
    int restart_flag = 0;
    if (restart_file.empty())
    {
        mfix.InitLevelData(time);
    }
    else
    {
        restart_flag = 1;
        // NOTE: mfix::levelset_restart == true loading level-set from a
        // checkpoint file. However, if this is a replicating restart,
        // mfix::levelset_restart is set to false again (so that the level-sets
        // are recomputed for the replicated system).
        mfix.levelset_restart = true;

        // NOTE: during replication 1) this also re-builds ebfactories and
        // level-set 2) this can change the grids
        IntVect Nrep(repl_x,repl_y,repl_z);
        mfix.Restart(restart_file, &nstep, &dt, &time, Nrep);
    }

    if (mfix.fluid.solve){
      mfix.init_advection();
    
      //amrex::Abort("111");

      mfix.mfix_init_solvers();
    }

    // This checks if we want to regrid
    if (!mfix.IsSteadyState() && regrid_int > -1 && nstep%regrid_int == 0)
    {
        amrex::Print() << "Regridding at step " << nstep << std::endl;
        mfix.Regrid();
    }

    if ((DEM::solve || PIC::solve) && write_ls)
        mfix.WriteStaticPlotFile(static_plt_file);

    mfix.PostInit(dt, time, restart_flag, stop_time);

    if (mfix.fluid.solve)
      mfix.ReportGridStats();

    Real end_init = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(end_init, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
       std::cout << "Time spent in init      " << end_init << std::endl;

    int finish  = 0;

    // Initialize prev_dt here; it will be re-defined by call to evolve_fluid but
    // only if fluid.solve = T
    Real prev_dt = dt;

    // Write checkpoint and plotfiles with the initial data
    if ( (restart_file.empty() || plotfile_on_restart) &&
         (mfix::plot_int > 0 || mfix::plot_per_exact > 0 || mfix::plot_per_approx > 0) )
    {
      if (mfix.fluid.solve)
          mfix.mfix_compute_vort();
       mfix.WritePlotFile(plot_file, nstep, time);
    }

#ifdef AMREX_USE_ASCENT
    if ( (restart_file.empty() || ascent_on_restart) &&
         (mfix::ascent_int > 0 || mfix::ascent_per_approx > 0) ) {
      mfix.WriteAscentFile();
    }
#endif

    // We automatically write checkpoint files with the initial data
    //    if check_int > 0
    if ( restart_file.empty() && check_int > 0 )
    {
       mfix.WriteCheckPointFile(check_file, nstep, dt, time);
       last_chk = nstep;
    }

    // We automatically write ASCII files with the particle data
    //    if par_ascii_int > 0
    if ( par_ascii_int > 0 )
    {
       mfix.WriteParticleAscii(par_ascii_file, nstep);
       last_par_ascii = nstep;
    }

    if ( avg_int > 0 )
      {
        mfix.WriteAverageRegions(avg_file, nstep, time);
        last_avg = nstep;
      }

    bool do_not_evolve = !mfix.IsSteadyState() && ( (max_step == 0) ||
                     ( (stop_time >= 0.) && (time >  stop_time) ) ||
                     ( (stop_time <= 0.) && (max_step <= 0) ) );

    if (restart_file.empty())
    {
        amrex::Print() << " " << std::endl;
        bool unused_inputs = ParmParse::QueryUnusedInputs();
        if (unused_inputs)
           amrex::Print() << "We should think about aborting here due to unused inputs" << std::endl;
    }


    for (int lev = 0; lev <= mfix.finestLevel(); lev++)
    {
      if (DEM::restart_from_PIC) {

        PIC::solve = false;
        DEM::solve = true;

        mfix.PIC_to_DEM(lev);
      }
    }


    { // Start profiling solve here

        BL_PROFILE("mfix_solve");
        BL_PROFILE_REGION("mfix_solve");

        if ( !do_not_evolve)
        {
            while (finish == 0)
            {
                mfix.mfix_usr1(time);

                Real strt_step = ParallelDescriptor::second();

                if (!mfix.IsSteadyState() && regrid_int > -1 && nstep%regrid_int == 0)
                {
                   amrex::Print() << "Regridding at step " << nstep << std::endl;
                   mfix.Regrid();
                }

                mfix.Evolve(nstep, dt, prev_dt, time, stop_time);

                Real end_step = ParallelDescriptor::second() - strt_step;
                ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
                if (ParallelDescriptor::IOProcessor())
                    std::cout << "   Time per step        " << end_step << std::endl;

                if (!mfix.IsSteadyState())
                {
                    time += prev_dt;
                    nstep++;

                    writeNow(nstep, time, prev_dt, mfix);
                }

                // Mechanism to terminate MFIX normally.
                do_not_evolve =  mfix.IsSteadyState() || (
                     ( (stop_time >= 0.) && (time+0.1*dt >= stop_time) ) ||
                     ( max_step >= 0 && nstep >= max_step ) );
                if ( do_not_evolve ) finish = 1;
            }
        }
    }

    if (mfix.IsSteadyState())
        nstep = 1;

    // Dump plotfile at the final time
    if ( check_int > 0 && nstep != last_chk)
        mfix.WriteCheckPointFile(check_file, nstep, dt, time);
    if ( mfix::plot_int > 0)
        mfix.WritePlotFile(plot_file, nstep, time);
    if ( par_ascii_int > 0  && nstep != last_par_ascii)
        mfix.WriteParticleAscii(par_ascii_file, nstep);

    mfix.mfix_usr3();

    Real end_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(end_time, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Time spent in main (after init) " << end_time-end_init << std::endl;
        std::cout << "Time spent in main      " << end_time << std::endl;
    }

    amrex::Print() << " " << std::endl;
    //bool unused_inputs = ParmParse::QueryUnusedInputs(); UNUSED VARIABLE

    BL_PROFILE_REGION_STOP("mfix::main()");
    BL_PROFILE_VAR_STOP(pmain);

    } // This end bracket and the start bracket after Initialize are essential so
      // that the mfix object is deleted before Finalize

    amrex::Finalize();
    return 0;
}
