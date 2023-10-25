#include <mfix.H>
#include <mfix_rw.H>
#include <mfix_dem.H>
#include <mfix_pic.H>

#include <AMReX_ParmParse.H>
#include <AMReX_EBFArrayBox.H>

#include <string>
#include <sstream>

#ifdef MFIX_CATALYST
#include "catalyst.hpp"
#include "AMReX_Conduit_Blueprint.H"
#endif


using namespace amrex;


MFIXReadWrite::MFIXReadWrite (int nlev_in,
                amrex::Vector<amrex::BoxArray>& grids_in,
                amrex::Vector<amrex::Geometry>& geom_in,
                MFIXParticleContainer* pc_in,
                MFIXFluidPhase& fluid_in,
                amrex::Vector<std::unique_ptr<LevelData>>& m_leveldata_in,
                amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>& ebfactory_in,
                amrex::Vector<amrex::DistributionMapping>& dmap_in,
                bool ooo_debug_in,
                amrex::Vector<std::unique_ptr<amrex::MultiFab>>& level_sets_in,
                int levelset_refinement_in,
                int levelset_pad_in,
                int levelset_eb_refinement_in,
                int levelset_eb_pad_in,
                MFIXSolidsPhase& solids_in,
                MFIXDEM& dem,
                MFIXPIC& pic,
                MFIXReactions& reactions_in,
                amrex::Vector<amrex::MultiFab*>& particle_cost_in,
                amrex::Vector<amrex::MultiFab*>& particle_proc_in,
                amrex::Vector<amrex::MultiFab*>& fluid_proc_in,
                const amrex::Vector<amrex::IntVect>& ref_ratio_in,
                BCList& bc_list_in,
                Vector<std::unique_ptr<EBFArrayBoxFactory>>& particle_ebfactory_in,
                MFIXRegions& regions_in)
  : finest_level(nlev_in-1)
  , nlev(nlev_in)
  , grids(grids_in)
  , geom(geom_in)
  , pc(pc_in)
  , fluid(fluid_in)
  , m_leveldata(m_leveldata_in)
  , ebfactory(ebfactory_in)
  , dmap(dmap_in)
  , ooo_debug(ooo_debug_in)
  , level_sets(level_sets_in)
  , levelset_refinement(levelset_refinement_in)
  , levelset_pad(levelset_pad_in)
  , levelset_eb_refinement(levelset_eb_refinement_in)
  , levelset_eb_pad(levelset_eb_pad_in)
  , solids(solids_in)
  , m_dem(dem)
  , m_pic(pic)
  , reactions(reactions_in)
  , particle_cost(particle_cost_in)
  , particle_proc(particle_proc_in)
  , fluid_proc(fluid_proc_in)
  , ref_ratio(ref_ratio_in)
  , bc_list(bc_list_in)
  , particle_ebfactory(particle_ebfactory_in)
  , regions(regions_in)
  , m_ascent_actions_yaml("")
{
  readParameters();

#ifdef MFIX_CATALYST
  if (catalyst_enabled) {

    conduit_cpp::Node params;

    params["catalyst/scripts/script0"].set_string(catalyst_script);
    params["catalyst_load/implementation"].set_string(catalyst_impl);
    params["catalyst_load/search_paths/paraview"].set_string(catalyst_library_path);

    catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&params));

    if (err != catalyst_status_ok) {
      std::string message = " Error: Failed to initialize Catalyst!\n";
      std::cerr << message << err << std::endl;
      amrex::Print() << message;
      amrex::Abort(message);
    }
  }
#endif
}


MFIXReadWrite::~MFIXReadWrite ()
{
#ifdef MFIX_CATALYST
  if (catalyst_enabled) {
      conduit_node* f_params = conduit_node_create();
      catalyst_finalize(f_params);
  }
#endif
}



void MFIXReadWrite::readParameters ()
{
  {
     const std::string pp_root = "mfix";

     ParmParse pp(pp_root.c_str());

     // Checkpoint output control
     pp.query("checkpoint_files_output", checkpoint_files_output);
     pp.query("check_file", check_file);
     pp.query("check_int", check_int);

     //pp.query("check_per_exact", check_per_exact);
     pp.query("check_per_approx", check_per_approx);

     std::string walltime_interval_in;
     int has_walltime_interval = pp.query("check_walltime_interval", walltime_interval_in);

     if (has_walltime_interval) {
       int HH(0), MM(0), SS(0);
       if (sscanf(walltime_interval_in.c_str(), "%d:%d:%d", &HH, &MM, &SS) >= 2) {
         check_walltime_interval = static_cast<Real>(HH*3600 + MM*60 + SS);
       } else {
         std::string message =
           " Error: Unable to correctly parse checkpoint walltime interval "
           + walltime_interval_in + "\n" + " The correct format is HH:MM:SS\n";
         amrex::Print() << message;
         amrex::Abort(message);
       }
     }

     if (/*(check_int        > 0 && check_per_exact  > 0) ||*/
         (check_int        > 0 && check_per_approx        > 0) ||
         (check_int        > 0 && check_walltime_interval > 0) ||
         (check_per_approx > 0 && check_walltime_interval > 0) /*||
         (check_per_exact  > 0 && check_per_approx > 0) */ )
       amrex::Abort("Must choose only one of check_int, check_per_approx, or "
           "check_walltime_interval");

     pp.query("geom_chk_file", geom_chk_file);
     pp.query("geom_levelset_chk_file", geom_levelset_chk_file);
     pp.query("geom_chk_write", geom_chk_write);
     pp.query("geom_chk_read", geom_chk_read);
     pp.query("geom_chk_ccse_regtest", geom_chk_ccse_regtest);

     // Plot output control
     pp.query("plot_file", plot_file);
     pp.query("plotfile_on_restart", plotfile_on_restart);
     pp.query("plot_int", plot_int);

     //pp.query("plot_per_exact", plot_per_exact);
     pp.query("plot_per_approx", plot_per_approx);

     if (/*(plot_int       > 0 && plot_per_exact  > 0) ||*/
         (plot_int       > 0 && plot_per_approx > 0) /*||
         (plot_per_exact > 0 && plot_per_approx > 0) */ )
       amrex::Abort("Must choose only one of plot_int or plot_per_approx");

     // Plot solids only in specific regions
     {
       const std::string solids_pp_root = pp_root + ".solids";
       ParmParse ppSolids(solids_pp_root.c_str());

       std::vector<std::string> solids_regions;
       ppSolids.queryarr("regions", solids_regions);

       if (solids_regions.size() > 0 &&
           amrex::toLower(solids_regions[0]).compare("none") != 0) {

         m_solids_plot_regions.resize(solids_regions.size());

         // Loop over input plot regions
         for (size_t n(0); n < solids_regions.size(); n++) {

           auto& plot_region = m_solids_plot_regions[n];

           const std::string& region_name = solids_regions[n];
           plot_region.m_region_name = region_name;

           ppSolids.queryarr(region_name.c_str(), plot_region.m_plot_names);

           const std::string region_pp_root = solids_pp_root + "." + region_name;
           ParmParse ppSolidsRegion(region_pp_root.c_str());

           plot_region.m_pp_string = region_pp_root;

           int ppint = ppSolidsRegion.query("plot_int", plot_region.m_plot_int);
           int ppapprox = ppSolidsRegion.query("plot_per_approx", plot_region.m_plot_per_approx);

           AMREX_ALWAYS_ASSERT(ppint || ppapprox);

           ppSolidsRegion.queryarr("plt_fluid_vars", plot_region.m_plot_fluid_vars);
         }
       }
     }

     // Monitors
     m_monitors.initialize(pp_root, regions, m_leveldata, ebfactory, fluid,
         pc, particle_ebfactory, solids);

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
     pp.query("par_ascii_per_approx", par_ascii_per_approx);

     pp.query("restart", restart_file);

     ParmParse("pic2dem").query("convert", restart_file);

     pp.query("repl_x", repl_x);
     pp.query("repl_y", repl_y);
     pp.query("repl_z", repl_z);
     pp.query("regrid_int",regrid_int);

     if ( regrid_int == 0 )
       amrex::Abort("regrid_int must be > 0 or < 0");

     pp.queryarr("avg_p_g", avg_p_g);
     pp.queryarr("avg_ep_g", avg_ep_g);
     pp.queryarr("avg_vel_g", avg_vel_g);
     pp.queryarr("avg_T_g", avg_T_g);

     pp.queryarr("avg_ro_p", avg_ro_p);
     pp.queryarr("avg_vel_p", avg_vel_p);
     pp.queryarr("avg_T_p", avg_T_p);

     // MFIXRegions geometry
     pp.queryarr("avg_region_x_e", avg_region_x_e);
     pp.queryarr("avg_region_x_w", avg_region_x_w);
     pp.queryarr("avg_region_y_n", avg_region_y_n);
     pp.queryarr("avg_region_y_s", avg_region_y_s);
     pp.queryarr("avg_region_z_t", avg_region_z_t);
     pp.queryarr("avg_region_z_b", avg_region_z_b);
  }

  {
     ParmParse pp("mfix");

     pp.query("write_eb_surface", write_eb_surface);
     pp.query("write_ls", write_ls);
     pp.query("stop_for_unused_inputs", stop_for_unused_inputs);
     pp.query("only_print_grid_report", only_print_grid_report);
     pp.query("plt_geom", plt_geom);
  }

  {
     ParmParse pp("ascent");

     pp.query("actions", m_ascent_actions_yaml);
  }

#ifdef MFIX_CATALYST
  {
    ParmParse pp("catalyst");
    pp.query("catalyst_on_restart", catalyst_on_restart);
    pp.query("script", catalyst_script);
    pp.query("implementation", catalyst_impl);
    pp.query("library_path", catalyst_library_path);
    pp.query("enabled", catalyst_enabled);
  }
#endif
}


void
MFIXReadWrite::Initialize ()
{
  real_comp_names.clear();
  int_comp_names.clear();

  // Vectors of names for solids plot
  {
    real_comp_names.push_back("radius");
    real_comp_names.push_back("volume");
    real_comp_names.push_back("mass");
    real_comp_names.push_back("density");

    if (m_dem.solve()) {
      real_comp_names.push_back("omoi");
    } else {
      real_comp_names.push_back("ep_s");
    }

    real_comp_names.push_back("velx");
    real_comp_names.push_back("vely");
    real_comp_names.push_back("velz");

    if (m_dem.solve()){
      real_comp_names.push_back("omegax");
      real_comp_names.push_back("omegay");
      real_comp_names.push_back("omegaz");
    } else {
      real_comp_names.push_back("grad_tau_x");
      real_comp_names.push_back("grad_tau_y");
      real_comp_names.push_back("grad_tau_z");
    }

    real_comp_names.push_back("statwt");
    real_comp_names.push_back("dragcoeff");
    real_comp_names.push_back("dragx");
    real_comp_names.push_back("dragy");
    real_comp_names.push_back("dragz");

    real_comp_names.push_back("c_ps");
    real_comp_names.push_back("temperature");
    real_comp_names.push_back("convection");

    if (solids.solve_species())
      for (auto species: solids.species_names())
        real_comp_names.push_back("X_"+species);

    if (solids.solve_species() && reactions.solve()) {
      for (int n(0); n < solids.nspecies(); ++n) {
        auto species = solids.species_names(n);
        real_comp_names.push_back("chem_mass_txfr_"+species);
      }
    }

    if (reactions.solve()) {
      real_comp_names.push_back("chem_velx_txfr");
      real_comp_names.push_back("chem_vely_txfr");
      real_comp_names.push_back("chem_velz_txfr");
    }

    if (reactions.solve())
      real_comp_names.push_back("chem_h_txfr");

    int_comp_names.push_back("phase");
    int_comp_names.push_back("state");
#if MFIX_POLYDISPERSE
    int_comp_names.push_back("ptype");
#endif
  }

  // Finalize initialization of solids plot regions
  for (int n(0); n < m_solids_plot_regions.size(); ++n) {

    auto& plot_region = m_solids_plot_regions[n];

    const int N_names = plot_region.m_plot_names.size();
    if (N_names > 0) {

      plot_region.m_h_plot_types.clear();
      plot_region.m_h_plot_types.resize(N_names);

      for (int i(0); i < N_names; ++i) {
        const std::string& name = plot_region.m_plot_names[i];

        int found(0);

        for (int j(0); j < solids.names().size(); ++j) {
          const std::string& solids_name = solids.names(j);

          if (name.compare(solids_name) == 0) {
            plot_region.m_h_plot_types[i] = j+1;
            found = 1;
            break;
          }
        }

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(found, "Solid type not found in solids");
      }

      plot_region.m_d_plot_types.resize(N_names);
      Gpu::copy(Gpu::hostToDevice, plot_region.m_h_plot_types.begin(),
                plot_region.m_h_plot_types.end(), plot_region.m_d_plot_types.begin());
    }

    const RealBox* region_extents = regions.getRegion(plot_region.m_region_name);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(region_extents != nullptr, "Invalid solids plot region!");

    plot_region.m_region_extents = *region_extents;

    ParmParse pp(plot_region.m_pp_string.c_str());
    GetSolidsIOPltFlags(plot_region.m_write_real_comp, plot_region.m_write_int_comp);
  }

  // Initialize the mass balance variables
  m_mass_accum.resize(fluid.nspecies()*2, 0.);
  m_mass_inflow.resize(fluid.nspecies(), 0.);
  m_mass_outflow.resize(fluid.nspecies(), 0.);
  m_mass_prod.resize(fluid.nspecies(), 0.);


}


void
MFIXReadWrite::writeNow (MFIXTimer& timer,
                         Real dt,
                         bool first,
                         bool last)
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
        plot_test = test_per_approx(timer.time(), dt, plot_per_approx);

    }/*
    else if ( plot_per_exact  > 0 && (amrex::Math::abs(remainder(timer.time(), plot_per_exact)) < 1.e-12) )
    {
        plot_test = 1;
    }*/


    if ( (plot_test == 1) || ( ( plot_int > 0) && ( timer.nstep() %  plot_int == 0 ) ) )
    {
        if (fluid.solve() && mfix::m_run_type != RunType::PIC2DEM)
          ComputeVort();

        WritePlotFile(plot_file, timer.nstep(), timer.time());
    }


/*--------------------------------------------------------------------------------------------------
 *
 *                                     AMReX Solids Plot File Output Control
 *
 *------------------------------------------------------------------------------------------------*/
    if ((m_dem.solve() || m_pic.solve()) && (solids_plot_regions() == true)) {

      BL_PROFILE("mfix::WriteSolidsPlotFile()");

      for (int n(0); n <m_solids_plot_regions.size(); ++n) {

        auto& plot_region = m_solids_plot_regions[n];

        plot_test = 0;

        if (first) {
          if (//(restart_file.empty() || plotfile_on_restart) &&
              (plot_region.m_plot_int > 0 /*|| plot_region.m_plot_per_exact > 0*/ ||
               plot_region.m_plot_per_approx > 0))

            plot_test = 1;
        }

        else if (last && plot_region.m_plot_int > 0) {
          plot_test = 1;
        }

        else if (plot_region.m_plot_per_approx > 0.0)
        {
          plot_test = test_per_approx(timer.time(), dt, plot_region.m_plot_per_approx);

        }/*
        else if (plot_region.m_plot_per_exact  > 0 &&
                 (amrex::Math::abs(remainder(timer.time(), plot_region.m_plot_per_exact)) < 1.e-12) )
        {
            plot_test = 1;
        }*/


        if ((plot_test == 1) || ((plot_region.m_plot_int > 0) &&
            (timer.nstep() % plot_region.m_plot_int == 0)))
        {
          WriteSolidsPlotFile(plot_region, plot_solids_file, timer.nstep(), timer.time());
        }
      }
    }


/*--------------------------------------------------------------------------------------------------
 *
 *                                     MFIX Monitors Plot File Output Control
 *
 *------------------------------------------------------------------------------------------------*/
    for (int i(0); i < m_monitors.size(); ++i) {

      auto& monitor = m_monitors.get(i);

      monitor.reset_pc(pc);

      int monitor_test = 0;

      if ( first || (timer.nstep() == 0) ) {
        if ( (monitor.plot_int() > 0 || monitor.plot_per_approx() > 0) )
          monitor_test = 1;
      }

      else if ( last ) {
        monitor_test = 0;
      }

      else if (monitor.plot_per_approx() > 0.0) {
        monitor_test = test_per_approx(timer.time(), dt, monitor.plot_per_approx());
      }

      else if ((monitor.plot_int() > 0) && (timer.nstep() % monitor.plot_int() == 0)) {
        monitor_test = 1;
      }

      if ( monitor_test == 1 ) {

        monitor.write_csv(timer.time(), dt);
      }
    }


/*--------------------------------------------------------------------------------------------------
 *
 *                                       Ascent Output Control
 *
 *------------------------------------------------------------------------------------------------*/

#ifdef AMREX_USE_ASCENT
    int ascent_test = 0;

    if ( first ) {
        if ((restart_file.empty() || ascent_on_restart) &&
            (ascent_int > 0 || ascent_per_approx > 0) )
            ascent_test = 1;

    } else if (ascent_per_approx > 0.0) {
      ascent_test = test_per_approx(timer.time(), dt, ascent_per_approx);
    }

    if ( (ascent_test == 1) || ( ( ascent_int > 0) && ( timer.nstep() %  ascent_int == 0 ) ) )
    {
        const int myProc = ParallelDescriptor::MyProc();
        WriteAscentFile(timer.nstep(), timer.time());
    }
#endif


/*--------------------------------------------------------------------------------------------------
 *
 *                               AMReX checkpoint file output control
 *
 *------------------------------------------------------------------------------------------------*/

    if (checkpoint_files_output) {

      int check_test = 0;

      if ( check_int > 0 /* || check_per_exact > 0*/ ||
           check_per_approx > 0. || check_walltime_interval > 0.) {

        // We automatically write checkpoint files with the initial data
        if ( first ) {
          check_test = (restart_file.empty()) ? 1 : 0;
        }
        // We automatically write checkpoint files with the final data
        else if (last) {
          check_test = (timer.nstep() != last_chk) ? 1 : 0;
        }
        else if (check_per_approx > 0) {
          check_test = test_per_approx(timer.time(), dt, check_per_approx);
        }/*
        else if (check_per_exact > 0 &&
                 (Math::abs(remainder(timer.time(), check_per_exact)) < 1.e-12)) {
          check_test = 1;
        }*/
        else if (check_int > 0) {
          check_test = (timer.nstep() % check_int == 0) ? 1 : 0;
        }
        else if (check_walltime_interval > 0) {
          check_test = test_walltime_interval(timer);
        }

        if (check_test == 1) {

          Real time_start = timer.system_time();
          WriteCheckPointFile(check_file, timer.nstep(), dt, timer.time());
          Real chkpt_write_time = timer.elapsed_runtime(time_start);

          if (timer.walltime_limit() > 0) {
            timer.update_max_write_chkpt_time(chkpt_write_time);
          }

          last_chk = timer.nstep();
        }

        if ( timer.walltime_limit() > 0. && check_test == 0 ) {
          if ( !timer.runtime_left_is_sufficient() ) {
            WriteCheckPointFile(check_file, timer.nstep(), dt, timer.time());
          }
        }
      }

    }

/*--------------------------------------------------------------------------------------------------
 *
 *                               AMReX particle ASCII output control
 *
 *------------------------------------------------------------------------------------------------*/
    int par_ascii_test = 0;

    if ( par_ascii_int > 0) {
      if ( first || last ) {
        par_ascii_test = 1;
      } else if ( timer.nstep() %  par_ascii_int == 0 ) {
        par_ascii_test = 1;
      }

    } else if (par_ascii_per_approx > 0.0) {
      par_ascii_test = test_per_approx(timer.time(), dt,
          par_ascii_per_approx);

    }

    if( par_ascii_test == 1) {
      WriteParticleAscii(par_ascii_file, timer.nstep());
      last_par_ascii = timer.nstep();
    }

/*--------------------------------------------------------------------------------------------------
 *
 *                               MFIX averaging region output control
 *
 *------------------------------------------------------------------------------------------------*/

    int avg_region_test = 0;

    if ( avg_int > 0 ) {

      // Always write output with initial data
      if ( first ) {
        avg_region_test = 1;
      }
      // Do it for the last step
      else if ( last ) {
        avg_region_test = (timer.nstep() != last_avg) ? 1 : 0;
      }
      else {
        avg_region_test = (timer.nstep() % avg_int == 0) ? 1 : 0;
      }

      if ( avg_region_test == 1 ) {
        WriteAverageRegions( avg_file, timer.nstep(), timer.time() );
        last_avg = timer.nstep();
      }

    }


/*--------------------------------------------------------------------------------------------------
 *
 *                                  MFIX mass balance output control
 *
 *------------------------------------------------------------------------------------------------*/
    int mass_balance_report_test = 0;

    if ( mass_balance_report_int > 0 ) {

      if ( first ) { // Never write the first.
        mass_balance_report_test = 0;

      } else if (last) { // Always write the last.
        mass_balance_report_test = (timer.nstep() != last_mb_report) ? 1 :0;

      } else {
        mass_balance_report_test = (timer.nstep() % mass_balance_report_int == 0) ? 1 : 0;
      }
    }

    else if (mass_balance_report_per_approx > 0.0) {
      mass_balance_report_test = test_per_approx(timer.time(),
          dt, mass_balance_report_per_approx);
    }

    if ( mass_balance_report_test == 1) {
      WriteMassBalanceReport(timer.time());
      last_mb_report = timer.nstep();
    }

/*--------------------------------------------------------------------------------------------------
 *
 *                                  MFIX call to catalyst
 *
 *------------------------------------------------------------------------------------------------*/

#ifdef MFIX_CATALYST
    if (catalyst_enabled) {
      RunCatalystAdaptor(timer.nstep(), timer.time());
    }
#endif

}

void MFIXReadWrite::writeEBSurface() const
{
   if(write_eb_surface)
     WriteMyEBSurface();
}

void MFIXReadWrite::writeStaticPlotFiles() const
{
   if ((m_dem.solve() || m_pic.solve()) && write_ls) {
      WriteStaticPlotFileParticleLevelSet(static_plt_file_ls);
   }

   if (plt_geom) {
      WriteStaticPlotFileEBGeometry(static_plt_file_geom);
   }
}

void MFIXReadWrite::reportGridStats() const
{
   if (fluid.solve())
     ReportGridStats();
}

//
// Print the maximum values of the velocity components
//
void
MFIXReadWrite::mfix_print_max_vel (int lev,
                            const Vector<MultiFab*>& vel_g_in,
                            const Vector<MultiFab*>& p_g_in)
{
    amrex::Print() << "   max(abs(u/v/w/p))  = "
                   << vel_g_in[lev]->norm0(0,0,false,true) << "  "
                   << vel_g_in[lev]->norm0(1,0,false,true) << "  "
                   << vel_g_in[lev]->norm0(2,0,false,true) << "  "
                   << p_g_in[lev]->norm0(0,0,false,true) << std::endl;
}

//
// Print the maximum values of the pressure gradient components
//
void
MFIXReadWrite::mfix_print_max_gp (int lev,
                           const Vector<MultiFab*>& gp_g_in)
{
    amrex::Print() << "   max(abs(gpx/gpy/gpz))  = "
                   << gp_g_in[lev]->norm0(0,0,false,true) << "  "
                   << gp_g_in[lev]->norm0(1,0,false,true) << "  "
                   << gp_g_in[lev]->norm0(2,0,false,true) <<  std::endl;
}


//
// Determine if it is time to write based on approximate interval
//
int
MFIXReadWrite::test_per_approx (const Real time,
                         const Real dt,
                         const Real per_approx)
{
  // Check to see if we've crossed a _per_approx interval by comparing
  // the number of intervals that have elapsed for both the current
  // time and the time at the beginning of this timestep.

  int num_per_old = static_cast<int>( (time-dt) / per_approx );
  int num_per_new = static_cast<int>( (time   ) / per_approx );

  // Before using these, however, we must test for the case where we're
  // within machine epsilon of the next interval. In that case, increment
  // the counter, because we have indeed reached the next par_ascii_per_approx interval
  // at this point.

  const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * amrex::Math::abs(time);
  const Real next_time = (num_per_old + 1) * per_approx;

  if ((num_per_new == num_per_old) && amrex::Math::abs(time - next_time) <= eps) {
      num_per_new += 1;
  }

  // Similarly, we have to account for the case where the old time is within
  // machine epsilon of the beginning of this interval, so that we don't double
  // count that time threshold -- we already plotted at that time on the last timestep.

  if ((num_per_new != num_per_old) && amrex::Math::abs((time - dt) - next_time) <= eps)
      num_per_old += 1;

  return (num_per_old != num_per_new) ? 1 : 0;
}


//
// Determine if it is time to write before job is killed
//
int
MFIXReadWrite::test_walltime_interval (const MFIXTimer& timer)
{
  const Real walltime = m_interval_nb * check_walltime_interval;

  if (timer.elapsed_runtime() > walltime) {
    m_interval_nb++;
    return 1;
  }

  return 0;
}


//
//
//
void
MFIXReadWrite::ReportGridStats () const
{
  std::vector<long> counts(6,0);

  int lev = 0;

  const MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());

  // Count the number of regular cells
  counts[0] = static_cast<int>(amrex::ReduceSum(*volfrac, *(m_leveldata[lev]->ep_g), 0,
    [=] AMREX_GPU_HOST_DEVICE (Box const & bx,
                               Array4<const Real> const & vfrc,
                               Array4<const Real> const & ep) -> int
    {
      int dm = 0;

      amrex::Loop(bx, [vfrc,ep,&dm] (int i, int j, int k) noexcept
      {if(vfrc(i,j,k)==1.0) dm += 1;});

      return dm;
    }));

  // Count the number of covered cells
  counts[1] = static_cast<int>(amrex::ReduceSum( *volfrac, *(m_leveldata[lev]->ep_g), 0,
    [=] AMREX_GPU_HOST_DEVICE (Box const & bx,
                               Array4<const Real> const & vfrc,
                               Array4<const Real> const & ep) -> int
    {
      int dm = 0;

      amrex::Loop(bx, [vfrc,ep,&dm] (int i, int j, int k) noexcept
      {if(vfrc(i,j,k)==0.0) dm += 1;});

      return dm;
    }));

  // Count the number of cut cells
  counts[2] = static_cast<int>(amrex::ReduceSum( *volfrac, *(m_leveldata[lev]->ep_g), 0,
    [=] AMREX_GPU_HOST_DEVICE (Box const & bx,
                               Array4<const Real> const & vfrc,
                               Array4<const Real> const & ep) -> int
    {
      int dm = 0;

      amrex::Loop(bx, [vfrc,ep,&dm] (int i, int j, int k) noexcept
      {if(0.0 < vfrc(i,j,k) && vfrc(i,j,k) < 1.0) dm += 1;});

      return dm;
    }));

  int regular(0), covered(0), cut(0);

#ifdef _OPENMP
#pragma omp parallel reduction(+:regular, covered, cut) if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(*m_leveldata[lev]->ep_g,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const auto& epg_fab   =
      static_cast<EBFArrayBox const&>((*m_leveldata[lev]->ep_g)[mfi]);

    const Box& bx     = mfi.tilebox();
    const auto& flags = epg_fab.getEBCellFlagFab();

    // Count number of regular grids
    if (flags.getType(amrex::grow(bx,1)) == FabType::regular ) {
      regular += 1;
    } else if (flags.getType(amrex::grow(bx,1)) == FabType::covered ) {
      covered += 1;
    } else {
      cut += 1;
    }
  }

  counts[3] = regular;
  counts[4] = covered;
  counts[5] = cut;

  ParallelDescriptor::ReduceLongSum(counts.data(), 6);

  if(ParallelDescriptor::IOProcessor()){
    printf("\n\n****************************************\n");
    printf("  Coverage report:  Grids        Cells\n");
    printf("          regular:  %5ld   %10ld\n", counts[3], counts[0]);
    printf("          covered:  %5ld   %10ld\n", counts[4], counts[1]);
    printf("              cut:  %5ld   %10ld\n", counts[5], counts[2]);
    printf("****************************************\n\n");
  }
}

//
// Print the minimum volume fraction and cell location.
//
IntVect
MFIXReadWrite::mfix_print_min_epg ()
{

  ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
  ReduceData<int, int, int, int> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (int lev = 0; lev <= finest_level; lev++) {

    const Real tolerance = std::numeric_limits<Real>::epsilon();
    auto& ld = *m_leveldata[lev];
    const Real min_epg = ld.ep_g->min(0);

    for (MFIter mfi(*ld.vel_g,false); mfi.isValid(); ++mfi) {
      Box const& bx = mfi.tilebox();
      Array4<Real const> const& epg = ld.ep_g->const_array(mfi);

      reduce_op.eval(bx, reduce_data, [epg, min_epg, tolerance]
      AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        int found(0);
        int iloc(0);
        int jloc(0);
        int kloc(0);
        if( amrex::Math::abs(epg(i,j,k) - min_epg) < tolerance ){
          iloc = i;
          jloc = j;
          kloc = k;
          found = 1;
        }
        return {found, iloc, jloc, kloc};
      });

      ReduceTuple htuple = reduce_data.value();

      const int found(amrex::get<0>(htuple));
      if(found > 0){

        IntVect epg_cell = {amrex::get<1>(htuple),
                            amrex::get<2>(htuple),
                            amrex::get<3>(htuple)};

        amrex::Print(Print::AllProcs)
          << std::endl << std::endl << "min epg "  << min_epg
          << "  at " << epg_cell[0] << "  " << epg_cell[1] << "  " << epg_cell[2]
          << "   total found " << found << std::endl << std::endl;

        return epg_cell;

      }

      //AMREX_ALWAYS_ASSERT(min_epg > 0.275);

    } // mfi
  } // lev

  IntVect fake = {0,0,0};
  return fake;
}
