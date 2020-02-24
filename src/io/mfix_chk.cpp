#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_AmrCore.H>
#include <AMReX_buildInfo.H>
#include <AMReX_Geometry.H>

#include <mfix.H>
#include <mfix_F.H>
#include <MFIX_FLUID_Parms.H>
#include <MFIX_DEM_Parms.H>

namespace
{
    const std::string level_prefix {"Level_"};
}

// This function initializes the attributes vecVarsName,
//                                          pltscalarVars, pltscaVarsName,
//                                          chkscalarVars, chkscaVarsName.
// If new variables need to be added to the output/checkpoint, simply add them
// here and the IO routines will automatically take care of them.
void
mfix::InitIOChkData ()
{
    if (ooo_debug) amrex::Print() << "InitIOData" << std::endl;
    // Define the list of vector variables on faces that need to be written
    // to plotfile/checkfile.
    vecVarsName = {"u_g", "v_g", "w_g", "gpx", "gpy", "gpz"};

    chkscaVarsName = {"ep_g", "p_g", "ro_g", "rop_g", "mu_g", "level_sets"};

    Vector< MultiFab* > ep_g(nlev, nullptr);
    // TODO
    chkscalarVars  = {&ep_g,  &p_g,  &ro_g,  &ep_g,  &mu_g,  &level_sets};
}


void
mfix::WriteCheckHeader (const std::string& name,
                        int nstep,
                        Real dt,
                        Real time) const
{
   bool is_checkpoint = 1;

    if (ParallelDescriptor::IOProcessor())
    {
      std::string HeaderFileName(name + "/Header");
      VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
      std::ofstream HeaderFile;

      HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

      HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                      std::ofstream::trunc |
                      std::ofstream::binary);

      if ( ! HeaderFile.good() )
          amrex::FileOpenFailed(HeaderFileName);

      HeaderFile.precision(17);

      if (is_checkpoint)
         HeaderFile << "Checkpoint version: 1\n";
      else
         HeaderFile << "HyperCLaw-V1.1\n";

      const int nlevels = finestLevel()+1;
      HeaderFile << nlevels << "\n";

      // Time stepping controls
      HeaderFile << nstep << "\n";
      HeaderFile << dt << "\n";
      HeaderFile << time << "\n";

      Geometry geometry;

      // Geometry
      for (int i = 0; i < BL_SPACEDIM; ++i)
            HeaderFile << geometry.ProbLo(i) << ' ';
      HeaderFile << '\n';

      for (int i = 0; i < BL_SPACEDIM; ++i)
         HeaderFile << geometry.ProbHi(i) << ' ';
      HeaderFile << '\n';

      // BoxArray
      for (int lev = 0; lev < nlevels; ++lev)
      {
          boxArray(lev).writeOn(HeaderFile);
          HeaderFile << '\n';
      }
    }
}

void
mfix::WriteCheckPointFile (std::string& check_file,
                           int nstep,
                           Real dt,
                           Real time) const
{
    BL_PROFILE("mfix::WriteCheckPointFile()");

    const std::string& checkpointname = amrex::Concatenate( check_file, nstep );

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "\n\t Writing checkpoint " << checkpointname << std::endl;
    }

    const int nlevels = finestLevel()+1;
    amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, nlevels, true);

    WriteCheckHeader(checkpointname, nstep, dt, time);

    WriteJobInfo(checkpointname);
    if (FLUID::solve)
    {
       for (int lev = 0; lev < nlevels; ++lev) {

          // This writes all three velocity components
          VisMF::Write( (*vel_g[lev]),
            amrex::MultiFabFileFullPrefix(lev, checkpointname,
                  level_prefix, vecVarsName[0]));

          // This writes all three pressure gradient components
          VisMF::Write( (*gp[lev]),
            amrex::MultiFabFileFullPrefix(lev, checkpointname,
                  level_prefix, vecVarsName[3]));

          // Write scalar variables
          for (int i = 0; i < chkscalarVars.size(); i++ ) {
              if ( DEM::solve || (chkscaVarsName[i] != "level_sets"))
                 VisMF::Write( *((*chkscalarVars[i])[lev]),
                   amrex::MultiFabFileFullPrefix(lev, checkpointname,
                         level_prefix, chkscaVarsName[i]));
          }
       }
    }

    if ( DEM::solve )
    {
       pc -> Checkpoint(checkpointname, "particles");
    }


    if (DEM::solve)
    {
        // The level set might have a higher refinement than the mfix level.
        //      => Current mechanism for saving checkpoint files requires the
        //         same BoxArray for all MultiFabss on the same level
        // NOTE: the unrefined level-set (and the multi-level level set) are
        // both saved with the standard checkpoint file.
        std::stringstream raw_ls_name;
        raw_ls_name << checkpointname << "/ls_raw";

        // There is always a level 1 in the level_sets array:
        //    level_sets.size() == std::max(2, nlev)
        VisMF::Write( * level_sets[1], raw_ls_name.str() );

        // Also save the parameters necessary to re-build the LSFactory
        int levelset_params[] = { levelset_refinement,
                                  levelset_pad,
                                  levelset_eb_refinement,
                                  levelset_eb_pad         };

        std::ofstream param_file;
        std::stringstream param_file_name;
        param_file_name << checkpointname << "/LSFactory_params";
        param_file.open(param_file_name.str());
        amrex::writeIntData(levelset_params, 4, param_file);
   }
}
