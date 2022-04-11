#ifndef _MFIX_RW_H_
#define _MFIX_RW_H_

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_BoxArray.H>

#include <mfix_pc.H>
#include <mfix_fluid_parms.H>
#include <mfix_solids_parms.H>
#include <mfix_reactions_parms.H>


namespace MfixIO {

class MfixRW {

public:
   MfixRW (int nlev_in,
           amrex::Vector<amrex::BoxArray>& grids_in,
           amrex::Vector<amrex::Geometry>& geom_in,
           MFIXParticleContainer* pc_in,
           FluidPhase& fluid_in,
           amrex::Vector<std::unique_ptr<LevelData>>& m_leveldata_in,
           amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>& ebfactory_in,
           amrex::Vector<amrex::DistributionMapping>& dmap_in,
           bool ooo_debug_in,
           amrex::Vector<std::unique_ptr<amrex::MultiFab>>& level_sets_in,
           const amrex::Vector<amrex::BoxArray>& box_array_in,
           int levelset_refinement_in,
           int levelset_pad_in,
           int levelset_eb_refinement_in,
           int levelset_eb_pad_in,
           SolidsPhase& solids_in,
           Reactions& reactions_in,
           amrex::Vector<amrex::MultiFab*> particle_cost_in,
           amrex::Vector<amrex::MultiFab*> particle_proc_in,
           amrex::Vector<amrex::MultiFab*> fluid_proc_in,
           amrex::Real covered_val_in,
           const amrex::Vector<amrex::IntVect>& ref_ratio_in,
           amrex::Vector<const amrex::EB2::Level*>& eb_levels_in,
           int nghost_eb_basic_in,
           int nghost_eb_volume_in,
           int nghost_eb_full_in,
           amrex::EBSupport& m_eb_support_level_in,
           std::string load_balance_type_in,
           BCList& bc_list_in,
           amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>& particle_ebfactory_in);

   // TODO TODO TODO
   void set_pc (MFIXParticleContainer* pc) { this->pc = pc; };

   void InitIOChkData ();

   void ResetIOChkData ();

   void InitIOPltData ();

   void writeNow (int nstep,
                  amrex::Real time,
                  amrex::Real dt,
                  bool first=false, /* first call*/
                  bool last=false   /* last call*/);

   void writeStaticPlotFile () const;

   //! Save variables that don't change to plot file. The idea is that they can
   //! be saved _once_ per simulations and not repeatedly every Nth time step.
   void WriteStaticPlotFile (const std::string& plot_file_name) const;

   void writeEBSurface () const;

   void WriteMyEBSurface () const;

   void reportGridStats() const;

   void WriteCheckPointFile (std::string& check_file_name,
                             int nstep = 0,
                             amrex::Real dt = 0.0,
                             amrex::Real time = 0.0);

   void WritePlotFile (std::string& plot_file_name,
                       int nstep = 0,
                       amrex::Real time = 0.0);

   void WriteAscentFile (int nstep,
                         amrex::Real time) const;

   void WriteParticleAscii (std::string& par_ascii_file_name,
                            int nstep = 0) const;

   void WriteAverageRegions (std::string& avg_file,
                             int nstep,
                             amrex::Real time = 0.) const;

   void Restart (std::string& restart_chkfile,
                 int* nstep,
                 amrex::Real* dt,
                 amrex::Real* time,
                 amrex::IntVect& Nrep);

   void ComputeAverageFluidVars (const int lev,
                                 const amrex::Real time,
                                 const std::string& basename) const;

   void ComputeMassAccum (const int offset = 1);

   void ComputeMassProduction (const amrex::Real dt,
                               amrex::Vector<amrex::MultiFab const*> const& chem_txfr);

   void ComputeMassFlux (amrex::Vector<amrex::MultiFab const*> const& flux_x,
                         amrex::Vector<amrex::MultiFab const*> const& flux_y,
                         amrex::Vector<amrex::MultiFab const*> const& flux_z,
                         const int scomp,
                         const int ncomp,
                         const bool fluxes_are_area_weighted,
                         const amrex::Real dt);

   void ReportGridStats () const;

   void WriteMassBalanceReport (const amrex::Real time);

   int repl_x = 1;
   int repl_y = 1;
   int repl_z = 1;
   int max_step   = -1;
   int regrid_int = -1;
   amrex::Real stop_time  = -1.0;
   bool stop_for_unused_inputs = false;
   std::string restart_file {""};

   static int report_mass_balance;
   static int mass_balance_report_int;
   static amrex::Real mass_balance_report_per_approx;
   static amrex::Real mass_balance_report_time;

#ifdef MFIX_CATALYST
   std::string catalyst_script {""};
#endif

private:
   // Private members
   int finest_level;
   int nlev;

   amrex::Vector<amrex::BoxArray>& grids;

   amrex::Vector<amrex::Geometry>& geom;

   MFIXParticleContainer* pc;

   FluidPhase& fluid;

   // Options to control the computing of particle eulerian velocities
   amrex::Vector<amrex::Real> avg_p_g;
   amrex::Vector<amrex::Real> avg_ep_g;
   amrex::Vector<amrex::Real> avg_vel_g;
   amrex::Vector<amrex::Real> avg_T_g;

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

   const amrex::Vector<amrex::BoxArray>& box_array;

   int levelset_refinement;
   int levelset_pad;
   int levelset_eb_refinement;
   int levelset_eb_pad;

   SolidsPhase& solids;
   Reactions& reactions;

   amrex::Vector<amrex::MultiFab*>& particle_cost;
   amrex::Vector<amrex::MultiFab*>& particle_proc;
   amrex::Vector<amrex::MultiFab*>& fluid_proc;

   amrex::Real covered_val;

   const amrex::Vector<amrex::IntVect>& ref_ratio;

   //! EB levels representing fluid boundary conditions
   amrex::Vector<const amrex::EB2::Level*> eb_levels;

   int nghost_eb_basic;
   int nghost_eb_volume;
   int nghost_eb_full;
   amrex::EBSupport& m_eb_support_level;

   std::string load_balance_type;

   BCList& bc_list;

   amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>& particle_ebfactory;

   amrex::GpuArray<amrex::Real,2*SPECIES::NMAX> mass_accum = {0.};
   amrex::GpuArray<amrex::Real,  SPECIES::NMAX> mass_inflow = {0.};
   amrex::GpuArray<amrex::Real,  SPECIES::NMAX> mass_outflow = {0.};
   amrex::GpuArray<amrex::Real,  SPECIES::NMAX> mass_prod = {0.};


   void readParameters();

   bool write_eb_surface = false;

   // Flag for writing static levelset function in plot file
   bool write_ls         = false;
   std::string static_plt_file {"plt_ls"};


   // Checkpoint file controls
   bool checkpoint_files_output = true;
   std::string check_file {"chk"};
   int check_int = -1;
   int last_chk  = -1;


   // Plot file controls
   std::string plot_file {"plt"};
   bool plotfile_on_restart = false;

   int  plot_int = -1;
   amrex::Real plot_per_approx = -1.;
   amrex::Real plot_per_exact  = -1.;


   // Ascent controls
   bool ascent_on_restart = false;

   int ascent_int = -1;
   amrex::Real ascent_per_approx = -1.;


   // Particle ASCII output files (debug only)
   int par_ascii_int = -1;
   int last_par_ascii  = -1;
   std::string par_ascii_file {"par"};


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

   // Total number of variables to write in plot file
   int pltVarCount = 0;

   // Flags for saving particle data. By default, we have all flags on,
   // we turn off what we don't want in the init IO routine. This is somewhat
   // different from what we do with the fluid.
   amrex::Vector<int> write_real_comp = amrex::Vector<int>(AoSrealData::count+SoArealData::count,1);
   amrex::Vector<int> write_int_comp = amrex::Vector<int>(AoSintData::count+SoAintData::count,1);

public:
   // Variables to simplify checkpoint IO
   // amrex::Vector< amrex::Vector< amrex::MultiFab* > > vectorVars;
   amrex::Vector< std::string > vecVarsName;

   amrex::Vector< amrex::Vector< amrex::MultiFab* > > chkScalarVars;
   amrex::Vector< std::string > chkscaVarsName;

   amrex::Vector< amrex::Vector< amrex::MultiFab* > > chkTVars;
   amrex::Vector< std::string > chkTVarsName;

   amrex::Vector< amrex::Vector< amrex::MultiFab* > > chkSpeciesVars;
   amrex::Vector< std::string > chkSpeciesVarsName;

private:
   void WriteCheckHeader (const std::string& name,
                          int nstep,
                          amrex::Real dt,
                          amrex::Real time) const;

   void WriteJobInfo (const std::string& dir) const;

public:
   static void GotoNextLine (std::istream& is);
};

} // end of namespace MfixIO

#endif
