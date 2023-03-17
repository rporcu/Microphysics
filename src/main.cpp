#include <fstream>
#include <iomanip>

#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include <AMReX_buildInfo.H>

#include <mfix.H>

#include "build_info.H"
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_fluid.H>
#include <mfix_rw.H>

#ifdef MFIX_CATALYST
#include "catalyst.hpp"
#include "AMReX_Conduit_Blueprint.H"
#endif

// Set defaults that are different that what ARMeX uses.  We only
// add them if they are not already specified in the inputs file.
void add_par () {

   {
      // Set the extend domain flag by default
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

class FixInputs {
  public:
    // Constructor
    FixInputs (int is_IOProc)
      : m_is_IOProc(is_IOProc)
    {}

    template <class T>
    void fix (ParmParse& ppA,
              ParmParse& ppB,
              const char* input_string,
              bool deprecated = false)
    {
      if (ppA.contains(input_string)) {

        if (deprecated && m_is_IOProc) {
          std::string message = "amr." + std::string(input_string) +
            " is deprecated and support will be removed in the future. Use mfix." +
            std::string(input_string);
          amrex::Warning(message.c_str());
        }

        T input_value;
        ppA.get(input_string, input_value);

        if(ppB.contains(input_string))
          ppB.remove(input_string);

        ppB.add(input_string, input_value);
      }
    }

    template <class T>
    void fix_arr (ParmParse& ppA,
                  ParmParse& ppB,
                  const char* input_string,
                  bool deprecated = false)
    {
      if (ppA.contains(input_string)) {

        if (deprecated && m_is_IOProc) {
          std::string message = "amr." + std::string(input_string) +
            " is deprecated and support will be removed in the future. Use mfix." +
            std::string(input_string);
          amrex::Warning(message.c_str());
        }

        T input_value;
        ppA.getarr(input_string, input_value);

        if(ppB.contains(input_string))
          ppB.remove(input_string);

        ppB.addarr(input_string, input_value);
      }
    }

  private:
    const int m_is_IOProc;
};

void fix_par ()
{
  ParmParse pp_amr("amr");
  ParmParse pp_mfix("mfix");

  FixInputs fix_inputs(ParallelDescriptor::IOProcessor());

  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_regtest");
  fix_inputs.fix<bool>(pp_amr, pp_mfix, "checkpoint_files_output");
  fix_inputs.fix<std::string>(pp_amr, pp_mfix, "check_file");
  fix_inputs.fix<int>(pp_amr, pp_mfix, "check_int");
  fix_inputs.fix<std::string>(pp_amr, pp_mfix, "plot_file");
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plot_int");
  fix_inputs.fix<std::string>(pp_amr, pp_mfix, "restart");
}

void fix_par_for_backward_compatibility ()
{
  ParmParse pp_amr("amr");
  ParmParse pp_mfix("mfix");

  FixInputs fix_inputs(ParallelDescriptor::IOProcessor());

  fix_inputs.fix<Real>(pp_amr, pp_mfix, "repl_x", true);
  fix_inputs.fix<Real>(pp_amr, pp_mfix, "repl_y", true);
  fix_inputs.fix<Real>(pp_amr, pp_mfix, "repl_z", true);
  fix_inputs.fix<bool>(pp_amr, pp_mfix, "plotfile_on_restart", true);
  fix_inputs.fix<bool>(pp_amr, pp_mfix, "dual_grid", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "regrid_int", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "par_ascii_int", true);
  fix_inputs.fix<std::string>(pp_amr, pp_mfix, "par_ascii_file", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "avg_int", true);
  fix_inputs.fix<std::string>(pp_amr, pp_mfix, "avg_file", true);
  fix_inputs.fix_arr<Vector<Real>>(pp_amr, pp_mfix, "avg_vel_p", true);
  fix_inputs.fix_arr<Vector<Real>>(pp_amr, pp_mfix, "avg_p_g", true);
  fix_inputs.fix_arr<Vector<Real>>(pp_amr, pp_mfix, "avg_ep_g", true);
  fix_inputs.fix_arr<Vector<Real>>(pp_amr, pp_mfix, "avg_vel_g", true);
  fix_inputs.fix_arr<Vector<Real>>(pp_amr, pp_mfix, "avg_region_x_w", true);
  fix_inputs.fix_arr<Vector<Real>>(pp_amr, pp_mfix, "avg_region_x_e", true);
  fix_inputs.fix_arr<Vector<Real>>(pp_amr, pp_mfix, "avg_region_y_s", true);
  fix_inputs.fix_arr<Vector<Real>>(pp_amr, pp_mfix, "avg_region_y_n", true);
  fix_inputs.fix_arr<Vector<Real>>(pp_amr, pp_mfix, "avg_region_z_b", true);
  fix_inputs.fix_arr<Vector<Real>>(pp_amr, pp_mfix, "avg_region_z_t", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "ascent_int", true);
  fix_inputs.fix<Real>(pp_amr, pp_mfix, "ascent_per_approx", true);
  fix_inputs.fix_arr<Vector<std::string>>(pp_amr, pp_mfix, "regions", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_vel_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_ep_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_p_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_ro_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_MW_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_h_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_T_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_trac", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_cp_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_k_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_mu_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_diveu", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_vort", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_volfrac", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_gradp_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_X_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_D_g", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_cp_gk", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_h_gk", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_txfr", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_chem_txfr", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_proc", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_proc_p", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_cost_p", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_radius", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_volume", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_mass", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_ro_p", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_omoi", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_vel_p", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_omega_p", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_statwt", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_drag_p", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_cp_s", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_T_p", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_convection", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_X_s", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_vel_s_txfr", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_h_s_txfr", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_mass_sn_txfr", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_phase", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_state", true);
  fix_inputs.fix<int>(pp_amr, pp_mfix, "plt_ptype", true);
}

const char* HypreVersion ();
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
    amrex::Print() << "   MFIX git describe: " << githash_mfix<< "\n";
    amrex::Print() << "AMReX-Hydro git hash: " << HydroGitHash() << "\n";
    amrex::Print() << "     CSG-EB git hash: " << CsgEbGitHash() << "\n";
#ifdef AMREX_USE_HYPRE
    amrex::Print() << "       HYPRE Version: " << HypreVersion() << "\n";
#endif

    // Setting format to NATIVE rather than default of NATIVE_32
    FArrayBox::setFormat(FABio::FAB_NATIVE);

    // Call parameters fixing
    fix_par();
    fix_par_for_backward_compatibility();

    Real strt_time = ParallelDescriptor::second();

    Real time=0.0L;
    int nstep = 0;  // Current time step

    Real dt = -1.;

    // Default constructor. Note inheritance: mfix : AmrCore : AmrMesh
    //                                                             |
    //  => Geometry is constructed here: (constructs Geometry) ----+
    mfix mfix;

    MfixIO::MfixRW& mfixRW = *(mfix.mfixRW);

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
    params["catalyst/scripts/script0"].set_string(mfixRW.catalyst_script);
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
    mfixRW.writeEBSurface();

    if (mfixRW.only_print_grid_report)
    {
        mfix.InitLevelData(time);
        mfixRW.reportGridStats();

        Real end_init = ParallelDescriptor::second() - strt_time;
        ParallelDescriptor::ReduceRealMax(end_init, ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
           std::cout << "Time spent in init      " << end_init << std::endl;
    }
    else
    {

       if (mfix.m_dem.solve() || mfix.m_pic.solve())
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
           mfix.Restart(mfixRW.restart_file, nstep, dt, time, Nrep);
       }

       mfixRW.setReportTime(time);

       if (mfix.fluid.solve()){
         mfix.init_advection();

         mfix.mfix_init_solvers();
       }

       mfixRW.writeStaticPlotFiles();

       mfix.PostInit(dt, time, restart_flag, mfixRW.stop_time);

       mfixRW.reportGridStats();

       Real end_init = ParallelDescriptor::second() - strt_time;
       ParallelDescriptor::ReduceRealMax(end_init, ParallelDescriptor::IOProcessorNumber());

       if (ParallelDescriptor::IOProcessor())
          std::cout << "Time spent in init      " << end_init << std::endl;

       int finish  = 0;

       // Initialize prev_dt here; it will be re-defined by call to evolve_fluid but
       // only if fluid.solve() = T
       Real prev_dt = dt;

       if (mfixRW.restart_file.empty())
       {
           amrex::Print() << " " << std::endl;
           bool unused_inputs = ParmParse::QueryUnusedInputs();
           if (mfixRW.stop_for_unused_inputs && unused_inputs)
              amrex::Abort("Aborting here due to unused inputs");
           else if (unused_inputs)
              amrex::Print() << "We should think about aborting here due to unused inputs" << std::endl;
       }

       mfixRW.writeNow(nstep, time, dt, /*first=*/true, /*last=*/false);

       bool do_not_evolve = !mfix.IsSteadyState() && ( (mfixRW.max_step == 0) ||
                        ( (mfixRW.stop_time >= 0.) && (time >  mfixRW.stop_time) ) ||
                        ( (mfixRW.stop_time <= 0.) && (mfixRW.max_step <= 0) ) );

       mfixRW.ComputeMassAccum(0);


       for (int lev = 0; lev <= mfix.finestLevel(); lev++)
       {
         if (mfix.m_dem.restart_from_PIC()) {

           mfix.m_pic.set_solve(false);
           mfix.m_dem.set_solve(true);

           mfix.PIC_to_DEM(lev);
         }
       }

       { // Start profiling solve here

           BL_PROFILE("mfix_solve");
           //BL_PROFILE_REGION("mfix_solve");

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

                       mfixRW.writeNow(nstep, time, prev_dt);
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

       mfixRW.writeNow(nstep, time, dt, /*first=*/false, /*last=*/true);

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
    }

    } // This end bracket and the start bracket after Initialize are essential so
      // that the mfix object is deleted before Finalize
#ifdef MFIX_CATALYST
    conduit_node* f_params = conduit_node_create();
    catalyst_finalize(f_params);
#endif
    amrex::Finalize();
    return 0;
}
