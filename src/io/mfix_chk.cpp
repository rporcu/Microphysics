#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_AmrCore.H>
#include <AMReX_buildInfo.H>
#include <AMReX_Geometry.H>

#include <mfix_rw.H>
#include <mfix_fluid.H>
#include <mfix_dem.H>
#include <mfix_pic.H>

using namespace amrex;

namespace
{
    const std::string level_prefix {"Level_"};
}

namespace MfixIO {

// This function initializes the attributes vecVarsName,
//                                          pltscalarVars, pltscaVarsName,
//                                          chkScalarVars, chkscaVarsName.
// If new variables need to be added to the output/checkpoint, simply add them
// here and the IO routines will automatically take care of them.
void
MfixRW::InitIOChkData ()
{
    if (ooo_debug) amrex::Print() << "InitIOChkData" << std::endl;
    // Define the list of vector variables on faces that need to be written
    // to plotfile/checkfile.
    vecVarsName = {"u_g", "v_g", "w_g", "gpx", "gpy", "gpz"};

    chkscaVarsName = {"ep_g", "p_g", "ro_g", "level_sets"};
    
    chkTVarsName = {"T_g"};

    chkSpeciesVarsName = {"X_gk"};

    ResetIOChkData();
}


void
MfixRW::ResetIOChkData ()
{
  chkScalarVars.clear();
  chkScalarVars.resize(chkscaVarsName.size(), Vector< MultiFab*>(nlev));

  chkTVars.clear();
  chkTVars.resize(chkTVarsName.size(), Vector< MultiFab*>(nlev));

  chkSpeciesVars.clear();
  chkSpeciesVars.resize(chkSpeciesVarsName.size(), Vector< MultiFab*>(nlev));

  for (int lev(0); lev < nlev; ++lev) {
    chkScalarVars[0][lev] = m_leveldata[lev]->ep_g;
    chkScalarVars[1][lev] = m_leveldata[lev]->p_g;
    chkScalarVars[2][lev] = m_leveldata[lev]->ro_g;
    chkScalarVars[3][lev] = level_sets[lev].get();
    
    if (fluid.solve_enthalpy()) {
      chkTVars[0][lev] = m_leveldata[lev]->T_g;
      //chkTVars[1][lev] = m_leveldata[lev]->h_g;
    }

    if (fluid.solve_species()) {
      chkSpeciesVars[0][lev] = m_leveldata[lev]->X_gk;
    }
  }
}


void
MfixRW::WriteCheckHeader (const std::string& name,
                          int nstep,
                          Real dt,
                          Real time) const
{
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

      HeaderFile << "Checkpoint version: 1.1\n";

      const int nlevels = finest_level+1;
      HeaderFile << nlevels << "\n";

      // Time stepping controls
      HeaderFile << nstep << "\n";
      HeaderFile << dt << "\n";
      HeaderFile << time << "\n";

      // Geometry
      HeaderFile << RealVect(geom[0].ProbLo()) << "\n";
      HeaderFile << RealVect(geom[0].ProbHi()) << "\n";
      HeaderFile << geom[0].Domain().size() << "\n";

      Real small_volfrac(0.);
      ParmParse pp("eb2");
      pp.query("small_volfrac", small_volfrac);
      HeaderFile << small_volfrac << "\n";

      // BoxArray
      for (int lev = 0; lev < nlevels; ++lev)
      {
          grids[lev].writeOn(HeaderFile);
          HeaderFile << "\n";
      }
    }
}


void
MfixRW::WriteCheckPointFile (std::string& check_file_in,
                             int nstep,
                             Real dt,
                             Real time)
{
    BL_PROFILE("mfix::WriteCheckPointFile()");

    const std::string& checkpointname = amrex::Concatenate( check_file_in, nstep );

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "\n\t Writing checkpoint " << checkpointname << std::endl;
    }

    const int nlevels = finest_level+1;
    amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, nlevels, true);

    WriteCheckHeader(checkpointname, nstep, dt, time);

    WriteJobInfo(checkpointname);
    if (fluid.solve())
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
            if ( m_dem.solve() || m_pic.solve() || (chkscaVarsName[i] != "level_sets"))
                 VisMF::Write( *(chkScalarVars[i][lev]),
                   amrex::MultiFabFileFullPrefix(lev, checkpointname,
                         level_prefix, chkscaVarsName[i]));
          }

          if (fluid.solve_enthalpy()) {
             // Write temperature variables
             for (int i = 0; i < chkTVars.size(); i++) {
                VisMF::Write( *(chkTVars[i][lev]),
                  amrex::MultiFabFileFullPrefix(lev, checkpointname,
                        level_prefix, chkTVarsName[i]));
             }
          }

          if (fluid.solve_species()) {
             // Write species variables
             for (int i = 0; i < chkSpeciesVars.size(); i++) {
                VisMF::Write( *(chkSpeciesVars[i][lev]),
                  amrex::MultiFabFileFullPrefix(lev, checkpointname,
                        level_prefix, chkSpeciesVarsName[i]));
             }
          }
       }
    }

    if ( m_dem.solve() || m_pic.solve() )
    {
       pc->Checkpoint(checkpointname, "particles");
    }


    if (m_dem.solve() || m_pic.solve())
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

}
