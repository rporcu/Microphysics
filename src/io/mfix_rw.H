#ifndef _MFIX_RW_H_
#define _MFIX_RW_H_

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_BoxArray.H>

#include <mfix_timer.H>
#include <mfix_pc.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_regions.H>
#include <mfix_monitors.H>


class MFIXReadWrite
{
  public:
   MFIXReadWrite (int nlev_in,
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
                  amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>& particle_ebfactory_in,
                  MFIXRegions& regions_in);

   ~MFIXReadWrite ();

   // Class for specifying a subregion for ParticleContainer::WritePlotFile
   // and add a filter to plot only a selected group of particles
   struct SolidsPlotRegion {
     public:
       SolidsPlotRegion ()
         : m_region_name("")
         , m_region_extents(amrex::RealBox())
         , m_pp_string("")
         , m_plot_names(0)
         , m_h_plot_types(0)
         , m_d_plot_types(0)
         , m_plot_int(-1)
         , m_plot_per_approx(-1.)
         , m_plot_fluid_vars(0)
       {}

       std::string m_region_name;
       amrex::RealBox m_region_extents;

       std::string m_pp_string;

       amrex::Vector<std::string> m_plot_names;
       amrex::Gpu::HostVector<int> m_h_plot_types;
       amrex::Gpu::DeviceVector<int> m_d_plot_types;

       int m_plot_int;
       amrex::Real m_plot_per_approx;
       amrex::Vector<std::string> m_plot_fluid_vars;

       int m_last_solids_plt = -1;

       amrex::Vector<int> m_write_real_comp = amrex::Vector<int>(SoArealData::count,1);
       amrex::Vector<int> m_write_int_comp = amrex::Vector<int>(SoAintData::count,1);
   };

   void set_pc (MFIXParticleContainer* pc_in) { this->pc = pc_in; }

   void output (int estatus,
                int finish,
                int nstep,
                amrex::Real dt,
                amrex::Real time)
   { pc->output(estatus, finish, nstep, dt, time); }

   void InitIOChkData ();

   void ResetIOChkData ();

   void InitIOPltData ();

   void writeNow (MFIXTimer& timer,
                  bool first=false, /* first call*/
                  bool last=false)  /* last call*/
   { writeNow(timer, timer.dt(), first, last); }

   void writeNow (MFIXTimer& timer,
                  amrex::Real dt,
                  bool first=false, /* first call*/
                  bool last=false); /* last call*/

   void writeStaticPlotFiles () const;

   void writeEBSurface () const;

   void WriteMyEBSurface () const;

   void reportGridStats () const;

   void WriteCheckPointFile (std::string& check_file_name,
                             int nstep,
                             amrex::Real dt,
                             amrex::Real time);

   void WritePlotFile (std::string& plot_file_name,
                       int nstep = 0,
                       amrex::Real time = 0.0);

   void WriteSolidsPlotFile (SolidsPlotRegion& plot_region,
                             std::string& plot_file_name,
                             int nstep = 0,
                             amrex::Real time = 0.0);

   void WriteAscentFile (int nstep,
                         amrex::Real time) const;

#ifdef MFIX_CATALYST
    void RunCatalystAdaptor (int nstep,
                             amrex::Real time) const;
#endif

   void WriteParticleAscii (std::string& par_ascii_file_name,
                            int nstep = 0) const;

   void WriteAverageRegions (std::string& avg_file,
                             int nstep,
                             amrex::Real time = 0.) const;

   void ComputeVort ();

   void ComputeAverageFluidVars (const int lev,
                                 const amrex::Real time,
                                 const std::string& basename) const;

   void InitMassBalance (const int nspecies);

   void ComputeMassAccum (const int offset = 1);

   void ComputeMassProduction (const amrex::Real dt,
                               amrex::Vector<amrex::MultiFab const*> const& txfr);

   void ComputeMassFlux (amrex::Vector<amrex::MultiFab const*> const& flux_x,
                         amrex::Vector<amrex::MultiFab const*> const& flux_y,
                         amrex::Vector<amrex::MultiFab const*> const& flux_z,
                         const int scomp,
                         const int ncomp,
                         const bool fluxes_are_area_weighted,
                         const int eb_has_flow,
                         amrex::Vector<amrex::MultiFab const*> const& eb_vel,
                         amrex::Vector<amrex::MultiFab const*> const& eb_species,
                         const amrex::Real dt);

   void ReportGridStats () const;

   void WriteMassBalanceReport (const amrex::Real time);

   void mfix_print_max_vel (int lev,
                            const amrex::Vector<amrex::MultiFab*>& vel_g_in,
                            const amrex::Vector<amrex::MultiFab*>& p_g_in);

   void mfix_print_max_gp (int lev,
                           const amrex::Vector<amrex::MultiFab*>& gp_g_in);

   amrex::IntVect mfix_print_min_epg ();

   void GotoNextLine (std::istream& is);

   int repl_x = 1;
   int repl_y = 1;
   int repl_z = 1;
   int regrid_int = -1;
   bool stop_for_unused_inputs = false;
   bool only_print_grid_report = false;
   std::string restart_file {""};

   // Geometry checkpoint file handling inputs
   std::string geom_chk_file {"geom_chk"};
   std::string geom_levelset_chk_file {"geom_levelset_chk"};
   bool geom_chk_write = false;
   bool geom_chk_read = false;
   bool geom_chk_ccse_regtest = false; // Used when running CCSE test for EB checkpoint files

   int report_mass_balance = 0;
   int mass_balance_report_int = -1;
   int last_mb_report  = -1;
   amrex::Real mass_balance_report_per_approx = -1;
   amrex::Real mass_balance_report_time;

   void setReportTime (amrex::Real const a_time)
   { mass_balance_report_time = a_time; }

#ifdef MFIX_CATALYST
   bool catalyst_on_restart = false;
   std::string catalyst_script {""};
   std::string catalyst_impl {""};
   std::string catalyst_library_path {""};
   bool catalyst_enabled = false;
#endif

   // Variables to simplify checkpoint IO
   // amrex::Vector< amrex::Vector< amrex::MultiFab* > > vectorVars;
   amrex::Vector<std::string> vecVarsName;

   amrex::Vector<amrex::Vector<amrex::MultiFab*>> chkScalarVars;
   amrex::Vector<std::string> chkscaVarsName;

   amrex::Vector<amrex::Vector<amrex::MultiFab*>> chkTVars;
   amrex::Vector<std::string> chkTVarsName;

   amrex::Vector<amrex::Vector<amrex::MultiFab*>> chkSpeciesVars;
   amrex::Vector<std::string> chkSpeciesVarsName;

   void Initialize ();

    bool solids_plot_regions () const
    { return m_solids_plot_regions.size() > 0; }

private:
   // Private members
   int finest_level;
   int nlev;

   amrex::Vector<amrex::BoxArray>& grids;

   amrex::Vector<amrex::Geometry>& geom;

   MFIXParticleContainer* pc;

   MFIXFluidPhase& fluid;

   // Options to control the computing of particle eulerian velocities
   amrex::Vector<amrex::Real> avg_p_g;
   amrex::Vector<amrex::Real> avg_ep_g;
   amrex::Vector<amrex::Real> avg_vel_g;
   amrex::Vector<amrex::Real> avg_T_g;

   amrex::Vector<amrex::Real> avg_ro_p;
   amrex::Vector<amrex::Real> avg_vel_p;
   amrex::Vector<amrex::Real> avg_T_p;

   amrex::Vector<amrex::Real> avg_region_x_e;
   amrex::Vector<amrex::Real> avg_region_x_w;
   amrex::Vector<amrex::Real> avg_region_y_s;
   amrex::Vector<amrex::Real> avg_region_y_n;
   amrex::Vector<amrex::Real> avg_region_z_b;
   amrex::Vector<amrex::Real> avg_region_z_t;

   amrex::Vector<std::unique_ptr<LevelData>>& m_leveldata;
   amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>& ebfactory;
   amrex::Vector<amrex::DistributionMapping>& dmap;

   bool ooo_debug;

   amrex::Vector<std::unique_ptr<amrex::MultiFab>>& level_sets;

   int levelset_refinement;
   int levelset_pad;
   int levelset_eb_refinement;
   int levelset_eb_pad;

   MFIXSolidsPhase& solids;

   MFIXDEM& m_dem;

   MFIXPIC& m_pic;

   MFIXReactions& reactions;

   amrex::Vector<amrex::MultiFab*>& particle_cost;
   amrex::Vector<amrex::MultiFab*>& particle_proc;
   amrex::Vector<amrex::MultiFab*>& fluid_proc;

   const amrex::Vector<amrex::IntVect>& ref_ratio;

   BCList& bc_list;

   amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>& particle_ebfactory;

   MFIXRegions& regions;

   amrex::Vector<amrex::Real> m_mass_accum;
   amrex::Vector<amrex::Real> m_mass_inflow;
   amrex::Vector<amrex::Real> m_mass_outflow;
   amrex::Vector<amrex::Real> m_mass_prod;

   void readParameters ();

   int test_per_approx(const amrex::Real time,
                       const amrex::Real dt,
                       const amrex::Real per_approx);

   int test_walltime_interval (const MFIXTimer& timer);

   bool write_eb_surface = false;

   // Flag for writing static particle levelset function in plot file
   bool write_ls = false;
   std::string static_plt_file_ls {"plt_ls"};

   // Flag for writing static EB geometry data in plot file
   bool plt_geom = false;
   std::string static_plt_file_geom {"plt_eb_geom"};

   // Checkpoint file controls
   bool checkpoint_files_output = true;
   std::string check_file {"chk"};

   int check_int = -1;
   amrex::Real check_per_approx = -1.;
   amrex::Real check_walltime_interval = -1.;
   int m_interval_nb = 1;
   //amrex::Real check_per_exact  = -1.;
   int last_chk  = -1;

   // Plot file controls
   std::string plot_file {"plt"};
   std::string plot_solids_file {"solids"};
   bool plotfile_on_restart = false;

   int plot_int = -1;
   amrex::Real plot_per_approx = -1.;
   //amrex::Real plot_per_exact  = -1.;

   // Ascent controls
   bool ascent_on_restart = false;

   int ascent_int = -1;
   amrex::Real ascent_per_approx = -1.;

   // Particle ASCII output files (debug only)
   int par_ascii_int = -1;
   int last_par_ascii  = -1;
   std::string par_ascii_file {"par"};
   amrex::Real par_ascii_per_approx = -1;

   // Averaging regions
   int avg_int = -1;
   int last_avg = -1;
   std::string avg_file {"avg_region"};

   // Last step at which we wrote a plotfile
   int last_plt = -1;

   // Flags for saving fluid data in plot files
   int plt_vel_g     = 1;
   int plt_ep_g      = 1;
   int plt_p_g       = 0;
   int plt_ro_g      = 0;
   int plt_MW_g      = 0;
   int plt_trac      = 0;
   int plt_cp_g      = 0;
   int plt_T_g       = 0;
   int plt_h_g       = 0;
   int plt_k_g       = 0;
   int plt_mu_g      = 0;
   int plt_X_gk      = 0;
   int plt_D_gk      = 0;
   int plt_cp_gk     = 0;
   int plt_h_gk      = 0;
   int plt_diveu     = 0;
   int plt_volfrac   = 0;
   int plt_gradp_g   = 0;
   int plt_vort      = 0;
   int plt_txfr      = 0;
   int plt_chem_txfr = 0;
   int plt_proc      = 0;
   int plt_proc_p    = 0;
   int plt_cost_p    = 0;
   // Geometry data
   int plt_impf      = 0;
   int plt_centroid  = 0;
   int plt_bndryarea = 0;
   int plt_bndrycent = 0;
   int plt_bndrynorm = 0;
   int plt_areafrac  = 0;
   int plt_facecent  = 0;
   int plt_edgecent  = 0;

   // Total number of variables to write in plot file
   int pltVarCount = 0;

   // Flags for saving particle data. By default, we have all flags on,
   // we turn off what we don't want in the init IO routine. This is somewhat
   // different from what we do with the fluid.
   amrex::Vector<int> write_real_comp = amrex::Vector<int>(SoArealData::count,1);
   amrex::Vector<int> write_int_comp = amrex::Vector<int>(SoAintData::count,1);

   amrex::Vector<std::string> real_comp_names;
   amrex::Vector<std::string>  int_comp_names;

   amrex::Vector<SolidsPlotRegion> m_solids_plot_regions;

   void GetSolidsIOPltFlags (amrex::Vector<int>& write_real_comp_out,
                             amrex::Vector<int>& write_int_comp_out);

   void WriteCheckHeader (const std::string& name,
                          int nstep,
                          amrex::Real dt,
                          amrex::Real time) const;

   void WriteJobInfo (const std::string& dir) const;

   //! Save variables that don't change to plot file. The idea is that they can
   //! be saved _once_ per simulations and not repeatedly every Nth time step.
   void WriteStaticPlotFileParticleLevelSet (const std::string& plotfilename) const;
   void WriteStaticPlotFileEBGeometry (const std::string& plotfilename) const;

   Monitors m_monitors;

   std::string m_ascent_actions_yaml;
};

#endif
