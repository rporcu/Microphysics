#include <fstream>
#include <iomanip>

#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_buildInfo.H>

#include <mfix.H>
#include <restarter.H>

//#include <build_info.H>
#include <mfix_dem.H>
#include <mfix_fluid.H>
#include <mfix_rw.H>


// Set defaults that are different that what ARMeX uses.  We only
// add them if they are not already specified in the inputs file.
void add_par () {

  {// Set the extend domain flag by default
    ParmParse pp("eb2");
    if(!pp.contains("extend_domain_face")) {
       pp.add("extend_domain_face",true);
    }
  }

  {
    // Disable managed memory for GPUs
    ParmParse pp("amrex");
    if(!pp.contains("the_arena_is_managed")) {
      pp.add("the_arena_is_managed",0);
    }
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
  amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, add_par);

  {
    // Write out the MFIX git hash (the AMReX git hash is already written)
    const char* githash_mfix = buildInfoGetGitHash(1);
    amrex::Print() << "   MFIX git describe: " << githash_mfix<< "\n";
//    amrex::Print() << "     CSG-EB git hash: " << CsgEbGitHash() << "\n";

    // Setting format to NATIVE rather than default of NATIVE_32
    FArrayBox::setFormat(FABio::FAB_NATIVE);

    Real strt_time = ParallelDescriptor::second();

    Real time=0.0L;
    int nstep = 0;
    Real dt = -1.;

    // Default constructor. Inheritance: mfix : AmrCore : AmrMesh
    //                                                             |
    //  => Geometry is constructed here: (constructs Geometry) ----+
    mfix mfix_coarse;
    mfix mfix_fine;
    MFIXRestarter mfix_restarter(mfix_coarse.nlev);

    auto& rw_coarse = *(mfix_coarse.mfixRW);
    auto& rw_fine = *(mfix_fine.mfixRW);

    // Initialize internals from ParamParse database
    mfix_coarse.InitParams();
    mfix_fine.InitParams();

    // Initialize memory for data-array internals
    mfix_coarse.ResizeArrays();

    // Initialize EB geometry. This needs to be done before grid creation (in
    // mfix::Init), as the grids are created using each EB-level's volfrac.
    mfix_coarse.make_eb_geometry();

    // Initialize derived internals
    bool init_fluid_grids = true;
    mfix_coarse.Init(time, init_fluid_grids);

    // Create EB factories on new grids
    mfix_coarse.make_eb_factories();

//    // Write out EB surface
//    mfix_coarse.mfixRW->writeEBSurface();

    if (mfix_coarse.m_dem.solve())
    {
      // Fill level-sets on each level
      mfix_coarse.fill_eb_levelsets();
    }

    // NOTE: mfix::levelset_restart == true loading level-set from a
    // checkpoint file. However, if this is a replicating restart,
    // mfix::levelset_restart is set to false again (so that the level-sets
    // are recomputed for the replicated system).
    mfix_coarse.levelset_restart = true;
    mfix_fine.levelset_restart = false;
    rw_fine.restart_file.clear();

    // NOTE: during replication 1) this also re-builds ebfactories and
    // level-set 2) this can change the grids
    amrex::IntVect Nrep(1,1,1);
    mfix_coarse.Restart(rw_coarse.restart_file, &nstep, &dt, &time, Nrep);

//    // This checks if we want to regrid
//    if (rw_coarse.regrid_int > -1)
//    {
//      amrex::Print() << "Regridding" << std::endl;
//      mfix_coarse.Regrid();
//    }

    mfix_restarter.allocate_coarse_arrays(&mfix_coarse);

    // Set fine mesh objects
    mfix_restarter.set_fine_objects(&mfix_fine, &mfix_coarse);

    // Initialize memory for data-array internals
    mfix_fine.ResizeArrays();

    // Initialize EB geometry. This needs to be done before grid creation (in
    // mfix::Init), as the grids are created using each EB-level's volfrac.
    mfix_fine.make_eb_geometry();

    // Initialize derived internals
    init_fluid_grids = false;
    mfix_fine.Init(time, init_fluid_grids);

    // Create EB factories on new grids
    mfix_fine.make_eb_factories();

    rw_fine.writeEBSurface();

    if (mfix_fine.m_dem.solve())
    {
      mfix_fine.fill_eb_levelsets();
    }

    const int dem_solve_flag = mfix_fine.m_dem.solve();
    AMREX_ALWAYS_ASSERT(mfix_fine.m_pic.solve() == 0);

    mfix_fine.m_dem.set_solve(0);
    mfix_fine.InitLevelData(time);

    // Reset dem_solve flag to the previous value
    mfix_fine.m_dem.set_solve(dem_solve_flag);

    mfix_restarter.allocate_fine_arrays(&mfix_fine);

    mfix_restarter.calc_txfr(&mfix_coarse, mfix_restarter.avgdPIC_coarse, time);

    // Free memory
    delete mfix_coarse.pc;
    mfix_coarse.pc = nullptr;

    // Here starts the part with new stuff
    // Set mfix_fine dem_solve to false so we initialize only the fluid data
    // add here the copy of fluid's coarse to fine variables in here
    mfix_restarter.txfr_fluid_data(&mfix_coarse, &mfix_fine);

    mfix_restarter.generate_particles(&mfix_coarse, &mfix_fine);

    // Free coarse fluid data
    for (int lev(0); lev < mfix_coarse.nlev; ++lev) {
      auto obj =  mfix_coarse.m_leveldata[lev].release();
      delete obj;
    }

    int restart_flag(1);

    mfix_fine.Restart(rw_fine.restart_file, &nstep, &dt, &time, Nrep);

    if (mfix_fine.fluid.solve()) {
      mfix_fine.init_advection();
      mfix_fine.mfix_init_solvers();
    }

    ///////// Post Init
    mfix_fine.PostInit(dt, time, restart_flag, time);

    mfix_restarter.init_particles_data(&mfix_fine);

    mfix_fine.Evolve(nstep, dt, dt, time, time+dt);

    rw_fine.reportGridStats();

    if (rw_fine.stop_for_unused_inputs && ParmParse::QueryUnusedInputs())
      amrex::Abort("Aborting here due to unused inputs");

    rw_fine.writeNow(nstep, time, dt, /*first=*/false, /*last=*/true);

    Real end_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(end_time, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
      std::cout << "Time spent in restarter " << end_time << std::endl;
    }

    amrex::Print() << " " << std::endl;

  } // This end bracket and the start bracket after Initialize are essential so
    // that the mfix object is deleted before Finalize

  amrex::Finalize();

  return 0;
}
