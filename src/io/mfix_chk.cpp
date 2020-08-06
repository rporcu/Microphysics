#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_AmrCore.H>
#include <AMReX_buildInfo.H>
#include <AMReX_Geometry.H>

#include <mfix.H>
#include <mfix_fluid_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>

namespace
{
    const std::string level_prefix {"Level_"};
}

// This function initializes the attributes vecVarsName,
//                                          pltscalarVars, pltscaVarsName,
//                                          chkScalarVars, chkscaVarsName.
// If new variables need to be added to the output/checkpoint, simply add them
// here and the IO routines will automatically take care of them.
void
mfix::InitIOChkData ()
{
    if (ooo_debug) amrex::Print() << "InitIOChkData" << std::endl;
    // Define the list of vector variables on faces that need to be written
    // to plotfile/checkfile.
    vecVarsName = {"u_g", "v_g", "w_g", "gpx", "gpy", "gpz"};

    chkscaVarsName = {"ep_g", "p_g", "ro_g", "level_sets"};
    //chkscaVarsName = {"ep_g", "p_g", "ro_g", "level_sets", "MW_g", "mu_g"};
    
    chkTVarsName = {"T_g"};
    //chkTVarsName = {"T_g", "h_g", "cp_g", "k_g"};

    chkSpeciesVarsName = {"X_gk"};
    //chkSpeciesVarsName = {"X_gk", "D_gk"};

    chkSpeciesTVarsName = {"h_gk"};
    //chkSpeciesTVarsName = {"h_gk", "cp_gk"};

    ResetIOChkData();
}


void
mfix::ResetIOChkData ()
{
  chkScalarVars.clear();
  chkScalarVars.resize(chkscaVarsName.size(), Vector< MultiFab**>(nlev));

  chkTVars.clear();
  chkTVars.resize(chkTVarsName.size(), Vector< MultiFab**>(nlev));

  chkSpeciesVars.clear();
  chkSpeciesVars.resize(chkSpeciesVarsName.size(), Vector< MultiFab**>(nlev));

  chkSpeciesTVars.clear();
  chkSpeciesTVars.resize(chkSpeciesTVarsName.size(), Vector< MultiFab**>(nlev));

  for (int lev(0); lev < nlev; ++lev) {
    chkScalarVars[0][lev] = &(m_leveldata[lev]->ep_g);
    chkScalarVars[1][lev] = &(m_leveldata[lev]->p_g);
    chkScalarVars[2][lev] = &(m_leveldata[lev]->ro_g);
    chkScalarVars[3][lev] = &level_sets[lev];
    //chkScalarVars[4][lev] = &(m_leveldata[lev]->MW_g);
    //chkScalarVars[5][lev] = &(m_leveldata[lev]->mu_g);
    
    if (advect_enthalpy) {
      chkTVars[0][lev] = &(m_leveldata[lev]->T_g);
      //chkTVars[1][lev] = &(m_leveldata[lev]->h_g);
      //chkTVars[2][lev] = &(m_leveldata[lev]->cp_g);
      //chkTVars[3][lev] = &(m_leveldata[lev]->k_g);
    }

    if (advect_fluid_species) {
      chkSpeciesVars[0][lev] = &(m_leveldata[lev]->X_gk);
      //chkSpeciesVars[2][lev] = &(m_leveldata[lev]->D_gk);
    }

    if (advect_fluid_species and advect_enthalpy) {
      chkSpeciesTVars[0][lev] = &(m_leveldata[lev]->h_gk);
      //chkSpeciesTVars[1][lev] = &(m_leveldata[lev]->cp_gk);
    }
  }
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
                           Real time)
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
       ResetIOChkData();

       for (int lev = 0; lev < nlevels; ++lev) {

          // This writes all three velocity components
          VisMF::Write( (*m_leveldata[lev]->vel_g),
            amrex::MultiFabFileFullPrefix(lev, checkpointname,
                  level_prefix, vecVarsName[0]));

          // This writes all three pressure gradient components
          VisMF::Write( (*m_leveldata[lev]->gp),
            amrex::MultiFabFileFullPrefix(lev, checkpointname,
                  level_prefix, vecVarsName[3]));

          // Write scalar variables
          for (int i = 0; i < chkScalarVars.size(); i++ ) {
            if ( DEM::solve or PIC::solve or (chkscaVarsName[i] != "level_sets"))
                 VisMF::Write( **(chkScalarVars[i][lev]),
                   amrex::MultiFabFileFullPrefix(lev, checkpointname,
                         level_prefix, chkscaVarsName[i]));
          }

          if (advect_enthalpy) {
             // Write temperature variables
             for (int i = 0; i < chkTVars.size(); i++) {
                VisMF::Write( **(chkTVars[i][lev]),
                  amrex::MultiFabFileFullPrefix(lev, checkpointname,
                        level_prefix, chkTVarsName[i]));
             }
          }

          if (advect_fluid_species) {
             // Write species variables
             for (int i = 0; i < chkSpeciesVars.size(); i++) {
                VisMF::Write( **(chkSpeciesVars[i][lev]),
                  amrex::MultiFabFileFullPrefix(lev, checkpointname,
                        level_prefix, chkSpeciesVarsName[i]));
             }
          }

          if (advect_fluid_species and advect_enthalpy) {
             // Write species energy variables
             for (int i = 0; i < chkSpeciesTVars.size(); i++) {
                VisMF::Write( **(chkSpeciesTVars[i][lev]),
                  amrex::MultiFabFileFullPrefix(lev, checkpointname,
                        level_prefix, chkSpeciesTVarsName[i]));
             }
          }
       }
    }

    if ( DEM::solve or PIC::solve )
    {
       pc->Checkpoint(checkpointname, "particles");
    }


    if (DEM::solve or PIC::solve)
    {
        // The level set might have a higher refinement than the mfix level.
        //      => Current mechanism for saving checkpoint files requires the
        //         same BoxArray for all MultiFabss on the same level
        // NOTE: the unrefined level-set (and the multi-level level set) are
        // both saved with the standard checkpoint file.
        std::stringstream raw_ls_name;
        raw_ls_name << checkpointname << "/ls_raw";

        // There is always a level 1 in the level_sets array:
        //    level_sets.size() == amrex::max(2, nlev)
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
