#ifndef _MFIX_RW_H_
#define _MFIX_RW_H_

#include <AMReX.H>

#include <mfix.H>

namespace MfixIO {

class MfixRW {

public:
   MfixRW() { readParameters(); }


   void writeNow (mfix& mfix, int nstep,
                  amrex::Real time,
                  amrex::Real dt,
                  bool first=false, /* first call*/
                  bool last=false   /* last call*/);

   void writeStaticPlotFile(const mfix &mfix) const;
   void writeEBSurface(mfix &mfix) const;
   void reportGridStats(const mfix &mfix) const;


   int repl_x = 1;
   int repl_y = 1;
   int repl_z = 1;
   int  max_step   = -1;
   int  regrid_int = -1;
   Real stop_time  = -1.0;
   bool stop_for_unused_inputs = false;
   std::string restart_file {""};

#ifdef MFIX_CATALYST
   std::string catalyst_script {""};
#endif

private:
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


};

} // end of namespace MfixIO

#endif
