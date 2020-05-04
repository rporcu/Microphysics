#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_AmrCore.H>

#include <AMReX_buildInfo.H>

#include <mfix.H>
#include <mfix_F.H>
#include <MFIX_FLUID_Parms.H>
#include <MFIX_DEM_Parms.H>
#include <MFIX_PIC_Parms.H>

using namespace std;

namespace
{
    const std::string level_prefix {"Level_"};
}

void
mfix::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

void mfix::WriteJobInfo (const std::string& dir) const
{
    if (ParallelDescriptor::IOProcessor())
    {
       // job_info file with details about the run
       std::ofstream jobInfoFile;
       std::string FullPathJobInfoFile = dir;
       std::string PrettyLine = "========================================="
                                "======================================\n";

       FullPathJobInfoFile += "/job_info";
       jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

       // job information
       jobInfoFile << PrettyLine;
       jobInfoFile << " MFIx Job Information\n";
       jobInfoFile << PrettyLine;

       jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
       jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

       jobInfoFile << "\n\n";

        // build information
       jobInfoFile << PrettyLine;
       jobInfoFile << " Build Information\n";
       jobInfoFile << PrettyLine;

       jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
       jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
       jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
       jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

       jobInfoFile << "\n";

       jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
       jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";

       jobInfoFile << "\n";

       jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
       jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

       jobInfoFile << "\n";

       jobInfoFile << "Fortran comp:  " << buildInfoGetFName() << "\n";
       jobInfoFile << "Fortran flags: " << buildInfoGetFFlags() << "\n";

       jobInfoFile << "\n";

       jobInfoFile << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
       jobInfoFile << "Libraries:     " << buildInfoGetLibraries() << "\n";

       const char* githash1 = buildInfoGetGitHash(1);
       const char* githash2 = buildInfoGetGitHash(2);
       if (strlen(githash1) > 0)
       {
           jobInfoFile << "MFIX git hash: " << githash1 << "\n";
       }
       if (strlen(githash2) > 0)
       {
           jobInfoFile << "AMReX git hash: " << githash2 << "\n";
       }

       jobInfoFile << "\n\n";

       // grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for (int i = 0; i <= finest_level; i++)
        {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (int n = 0; n < BL_SPACEDIM; n++)
            {
                jobInfoFile << geom[i].Domain().length(n) << " ";
            }
            jobInfoFile << "\n\n";
        }

        jobInfoFile << "\n\n";

        // runtime parameters
        jobInfoFile << PrettyLine;
        jobInfoFile << " Inputs File Parameters\n";
        jobInfoFile << PrettyLine;

        // ParmParse::dumpTable(jobInfoFile, true);

        jobInfoFile.close();
  }
}


void
mfix::WriteParticleAscii ( std::string& par_ascii_file, int nstep ) const
{
    BL_PROFILE("mfix::WriteParticleASCII()");

    if(DEM::solve or PIC::solve) {

        const std::string& par_filename = amrex::Concatenate(par_ascii_file,nstep);
        pc->WriteAsciiFile(par_filename);
    }
}


void
mfix::WriteAverageRegions ( std::string& avg_file, int nstep, Real time ) const
{
  BL_PROFILE("mfix::WriteAverageRegions()");

  for (int lev = 0; lev < nlev; lev++)
    {
      if (FLUID::solve) {
        ComputeAverageFluidVars( lev,
                                 time,
                                 avg_file);
      }

      //  Compute Eulerian velocities in selected regions
      if(DEM::solve or PIC::solve) {
        pc->ComputeAverageVelocities ( lev,
                                       time,
                                       avg_file,
                                       avg_vel_p,
                                       avg_region_x_w, avg_region_x_e,
                                       avg_region_y_s, avg_region_y_n,
                                       avg_region_z_b, avg_region_z_t );
      }
    }

}


void
mfix::ComputeAverageFluidVars ( const int lev, const Real time,
                                const std::string&  basename) const
{

  int nregions = avg_region_x_w.size();

  const int size_p_g   = avg_p_g.size();
  const int size_ep_g  = avg_ep_g.size();
  const int size_vel_g = avg_vel_g.size();

  const Real * dx   = geom[lev].CellSize();

  //
  // Check the regions are defined correctly
  //
  if ( ( avg_region_x_e.size() != nregions ) or
       ( avg_region_y_s.size() != nregions ) or
       ( avg_region_y_n.size() != nregions ) or
       ( avg_region_z_b.size() != nregions ) or
       ( avg_region_z_t.size() != nregions )  )
  {
    amrex::Print() << "ComputeAverageVelocities: some regions are not properly"
      " defined: skipping.";
    return;
  }

  const int var_count = 6;

  // Array to hold the data for global collection.
  std::vector<Real> regions_data(var_count*nregions, 0.0);

  Box domain(geom[lev].Domain());

  const amrex::MultiFab* volfrac = &(ebfactory[lev] -> getVolFrac());

  MultiFab& ep_g = *(m_leveldata[lev]->ep_g);

  // New multiFab to hold the cell center pressure.
  std::unique_ptr<MultiFab> pg_cc(new MultiFab(
             ep_g.boxArray(), ep_g.DistributionMap(),
             ep_g.nComp(),    ep_g.nGrow(), MFInfo(), *ebfactory[lev]));

  // Create a temporary nodal pressure multifab to sum in p_g and p0_g
  MultiFab pg_nd(m_leveldata[lev]->p_g->boxArray(), dmap[lev], 1, 0);
  pg_nd.setVal(0.);
  MultiFab::Copy(pg_nd, (*m_leveldata[lev]->p_g), 0, 0, 1, 0);
  MultiFab::Add (pg_nd, (*m_leveldata[lev]->p0_g), 0, 0, 1, 0);

  // Create a cell-center version of the combined pressure
  amrex::average_node_to_cellcenter(*pg_cc, 0, pg_nd, 0, 1);


  for ( int nr = 0; nr < nregions; ++nr )
  {

    // Create Real box for this region
    RealBox avg_region ( {avg_region_x_w[nr], avg_region_y_s[nr], avg_region_z_b[nr]},
                         {avg_region_x_e[nr], avg_region_y_n[nr], avg_region_z_t[nr]} );

    // Jump to next iteration if this averaging region is not valid
    if ( !avg_region.ok () )
    {
      amrex::Print() << "ComputeAverageVelocities: region "<< nr << 
        " is invalid: skipping\n";
      continue;
    }

    Real sum_vol    = 0.;
    Real sum_velx   = 0.;
    Real sum_vely   = 0.;
    Real sum_velz   = 0.;
    Real sum_p_g    = 0.;
    Real sum_ep_g   = 0.;

#if(0)
    amrex::Print() << "\n\n  Collecting field averages for region " << nr << "\n";
    if( size_p_g   > nr ) if( avg_p_g[nr]   == 1 ) amrex::Print() << "  > Gas pressure.......(p_g)\n";
    if( size_ep_g  > nr ) if( avg_ep_g[nr]  == 1 ) amrex::Print() << "  > Volume fraction....(ep_g)\n";
    if( size_vel_g > nr ) if( avg_vel_g[nr] == 1 ) amrex::Print() << "  > Gas velocity.......(m/s)\n";
#endif

    // Not tiling this loop.
    for (MFIter mfi(ep_g, false); mfi.isValid(); ++mfi)
    {

      const Box& bx  = mfi.validbox();

      // this is to check efficiently if this grid contains any eb stuff
      const EBFArrayBox&  epg_fab = static_cast<EBFArrayBox const&>(ep_g[mfi]);
      const EBCellFlagFab&  flags = epg_fab.getEBCellFlagFab();

      RealBox box_region ( bx, Geom(lev).CellSize (), Geom(lev).ProbLo() );

      if (flags.getType(bx) != FabType::covered)
      {
        if ( box_region.intersects ( avg_region ) )
        {

          // TODO: convert this in cpp
          mfix_collect_fluid(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_BOX(domain),
                             BL_TO_FORTRAN_ANYD(         ep_g[mfi]),
                             BL_TO_FORTRAN_ANYD((     *pg_cc)[mfi]),
                             BL_TO_FORTRAN_ANYD((*m_leveldata[lev]->vel_g)[mfi]),
                             BL_TO_FORTRAN_ANYD((   *volfrac)[mfi]),
                             &avg_region_x_w[nr], &avg_region_x_e[nr],
                             &avg_region_y_s[nr], &avg_region_y_n[nr],
                             &avg_region_z_b[nr], &avg_region_z_t[nr], dx,
                             &sum_ep_g, &sum_p_g,  &sum_vol,
                             &sum_velx, &sum_vely, &sum_velz);
        }
      }
    }

    regions_data[var_count*nr + 0] = sum_vol;
    regions_data[var_count*nr + 1] = sum_velx;
    regions_data[var_count*nr + 2] = sum_vely;
    regions_data[var_count*nr + 3] = sum_velz;
    regions_data[var_count*nr + 4] = sum_p_g;
    regions_data[var_count*nr + 5] = sum_ep_g;

  }

  // Compute parallel reductions
  ParallelDescriptor::ReduceRealSum ( regions_data.data(),  var_count*nregions );

  // Only the IO processor takes care of the output
  if (ParallelDescriptor::IOProcessor())
    {
      for ( int nr = 0; nr < nregions; ++nr )
        {

          Real sum_vol  = regions_data[var_count*nr + 0];
          Real sum_velx = regions_data[var_count*nr + 1];
          Real sum_vely = regions_data[var_count*nr + 2];
          Real sum_velz = regions_data[var_count*nr + 3];
          Real sum_p_g  = regions_data[var_count*nr + 4];
          Real sum_ep_g = regions_data[var_count*nr + 5];

          std::ofstream  ofs;
          std::string    fname;

          if( size_p_g > nr )
            {
              if( avg_p_g[nr] == 1) {

                fname = basename + "_p_g_" + std::to_string(nr) + ".dat";

                std::ifstream ifile(fname.c_str());
                bool exists = (bool)ifile;

                ofs.open ( fname.c_str(), ios::out | ios::app );
                if ( !ofs.good() ) amrex::FileOpenFailed ( fname );

                if ( !exists ) ofs << "#  Time   p_g  vol" << std::endl;
                ofs << time << " " << sum_p_g/sum_vol << "  " << sum_vol << std::endl;

                ofs.close();
            }

          }

          if( size_ep_g  > nr )
            {
              if( avg_ep_g[nr] == 1)
                {
                  fname = basename + "_ep_g_" + std::to_string(nr) + ".dat";

                  std::ifstream ifile(fname.c_str());
                  bool exists = (bool)ifile;

                  ofs.open ( fname.c_str(), ios::out | ios::app );
                  if ( !ofs.good() ) amrex::FileOpenFailed ( fname );

                  if ( !exists ) ofs << "#  Time   ep_g" << std::endl;
                  ofs << time << " " << sum_ep_g / sum_vol << std::endl;

                  ofs.close();

                }
            }

          if( size_vel_g > nr )
            {
              if( avg_vel_g[nr] )
                {
                  fname = basename + "_vel_g_" + std::to_string(nr) + ".dat";

                  std::ifstream ifile(fname.c_str());
                  bool exists = (bool)ifile;

                  ofs.open ( fname.c_str(), ios::out | ios::app );
                  if ( !ofs.good() ) amrex::FileOpenFailed ( fname );

                  if ( !exists ) ofs << "#  Time   velx  vely  velz" << std::endl;
                  ofs << time << " " << sum_velx / sum_vol << " "
                      <<                sum_vely / sum_vol << " "
                      <<                sum_velz / sum_vol << std::endl;

                  ofs.close();

                }
            }

        }

    }

}
