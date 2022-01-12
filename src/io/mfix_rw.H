#ifndef _MFIX_RW_H_
#define _MFIX_RW_H_

#include <AMReX.H>

#include <mfix.H>

namespace MfixIO {

class MfixRW {

public:
   MfixRW() { readParameters(); }


   void writeNow (int nstep, Real time, Real dt, mfix& mfix);

   void writeStaticPlotFile(const mfix &mfix) const;
   void writeAverageRegions(int nstep, Real time, const mfix &mfix);
   void writeAsentFile(mfix &mfix) const;
   void writeCheckPointFile(int nstep, Real dt, Real time, mfix &mfix);
   void writeParticleAscii(int nstep, const mfix &mfix);

   void writePlotFileInitial(int nstep, Real time, mfix &mfix);

   void writeCheckPointFileFinal(int nstep, Real dt, Real time, mfix &mfix);
   void writePlotFileFinal(int nstep, Real time, mfix &mfix);
   void writeParticleAsciiFinal(int nstep, const mfix &mfix);

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
   bool write_ls         = false;

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

   bool checkpoint_files_output = true;
};

} // end of namespace MfixIO

#endif
