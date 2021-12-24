#include <fstream>
#include <iomanip>

#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include <AMReX_buildInfo.H>

#include <mfix.H>

#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_rw.H>

#ifdef MFIX_CATALYST
#include "catalyst.hpp"
#include "AMReX_Conduit_Blueprint.H"
#endif

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

    MfixIO::MfixRW mfixRW;

    // Initialize internals from ParamParse database
    mfix.InitParams();

    // Initialize memory for data-array internals
    mfix.ResizeArrays();

    // Initialize EB geometry. This needs to be done before grid creation (in
    // mfix::Init), as the grids are created using each EB-level's volfrac.
    mfix.make_eb_geometry();

    // Initialize derived internals
    mfix.Init(time);

#ifdef MFIX_CATALYST
    conduit_cpp::Node params;
    params["catalyst/scripts/script0"].set_string(catalyst_script);
    params["catalyst_load/implementation"] = "paraview";
    params["catalyst_load/search_paths/paraview"] = "/home/corey/Builds/pvsb-dev/install/lib/catalyst";
    catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&params));
    if (err != catalyst_status_ok)
    {
        std::cerr << "Failed to initialize Catalyst: " << err << std::endl;
        return 1;
    }
#endif

    // Create EB factories on new grids
    mfix.make_eb_factories();

    // Write out EB sruface
    mfixRW.writeEBSurface(mfix);

    if (DEM::solve || PIC::solve)
    {
        // Fill level-sets on each level
        mfix.fill_eb_levelsets();
    }

    // Either init from scratch or from the checkpoint file
    int restart_flag = 0;
    if (mfixRW.restart_file.empty())
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
        IntVect Nrep(mfixRW.repl_x, mfixRW.repl_y, mfixRW.repl_z);
        mfix.Restart(mfixRW.restart_file, &nstep, &dt, &time, Nrep);
    }

    if (mfix.fluid.solve){
      mfix.init_advection();

      //amrex::Abort("111");

      mfix.mfix_init_solvers();
    }

    // This checks if we want to regrid
    if (!mfix.IsSteadyState() && mfixRW.regrid_int > -1 && nstep%mfixRW.regrid_int == 0)
    {
        amrex::Print() << "Regridding at step " << nstep << std::endl;
        mfix.Regrid();
    }

    mfixRW.writeStaticPlotFile(mfix);

    mfix.PostInit(dt, time, restart_flag, mfixRW.stop_time);

    mfixRW.reportGridStats(mfix);

    Real end_init = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(end_init, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
       std::cout << "Time spent in init      " << end_init << std::endl;

    int finish  = 0;

    // Initialize prev_dt here; it will be re-defined by call to evolve_fluid but
    // only if fluid.solve = T
    Real prev_dt = dt;

    mfixRW.writePlotFileInitial(nstep, time, mfix);
    mfixRW.writeCheckPointFile(nstep, dt, time, mfix);
    mfixRW.writeParticleAscii(nstep, mfix);
    mfixRW.writeAverageRegions(nstep, time, mfix);

    bool do_not_evolve = !mfix.IsSteadyState() && ( (mfixRW.max_step == 0) ||
                     ( (mfixRW.stop_time >= 0.) && (time >  mfixRW.stop_time) ) ||
                     ( (mfixRW.stop_time <= 0.) && (mfixRW.max_step <= 0) ) );

    mfix.ComputeMassAccum(0);

    if (mfixRW.restart_file.empty())
    {
        amrex::Print() << " " << std::endl;
        bool unused_inputs = ParmParse::QueryUnusedInputs();
        if (mfixRW.stop_for_unused_inputs && unused_inputs)
           amrex::Abort("Aborting here due to unused inputs");
        else if (unused_inputs)
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

                if (!mfix.IsSteadyState() && mfixRW.regrid_int > -1 && nstep%mfixRW.regrid_int == 0)
                {
                   amrex::Print() << "Regridding at step " << nstep << std::endl;
                   mfix.Regrid();
                }

                mfix.Evolve(nstep, dt, prev_dt, time, mfixRW.stop_time);

                Real end_step = ParallelDescriptor::second() - strt_step;
                ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
                if (ParallelDescriptor::IOProcessor())
                    std::cout << "   Time per step        " << end_step << std::endl;

                if (!mfix.IsSteadyState())
                {
                    time += prev_dt;
                    nstep++;

                    mfixRW.writeNow(nstep, time, prev_dt, mfix);
#ifdef MFIX_CATALYST
                    mfix.RunCatalystAdaptor(nstep, time);
#endif
                }

                // Mechanism to terminate MFIX normally.
                do_not_evolve =  mfix.IsSteadyState() || (
                     ( (mfixRW.stop_time >= 0.) && (time+0.1*dt >= mfixRW.stop_time) ) ||
                     ( mfixRW.max_step >= 0 && nstep >= mfixRW.max_step ) );
                if ( do_not_evolve ) finish = 1;
            }
        }
    }

    if (mfix.IsSteadyState())
        nstep = 1;

    // Dump plotfile at the final time
    mfixRW.writeCheckPointFileFinal(nstep, dt, time, mfix);
    mfixRW.writePlotFileFinal(nstep, time, mfix);
    mfixRW.writeParticleAsciiFinal(nstep, mfix);

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
#ifdef MFIX_CATALYST
    conduit_node* f_params = conduit_node_create();
    catalyst_finalize(f_params);
#endif
    amrex::Finalize();
    return 0;
}
