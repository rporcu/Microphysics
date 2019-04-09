#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_AmrCore.H>

#include <AMReX_buildInfo.H>

#include <AMReX_EBMultiFabUtil.H>

#include <mfix.H>
#include <mfix_F.H>

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
mfix::InitIOData ()
{
    if (ooo_debug) amrex::Print() << "InitIOData" << std::endl;
    // Define the list of vector variables on faces that need to be written
    // to plotfile/checkfile.
    vecVarsName = {"u_g", "v_g", "w_g", "gpx", "gpy", "gpz"};

    chkscaVarsName = {"ep_g", "p_g", "ro_g", "rop_g", "mu_g", "level_sets", "implicit_functions"};
    chkscalarVars  = {&ep_g,  &p_g,  &ro_g,  &ep_g,  &mu_g,  &level_sets,  &implicit_functions};
}

void
mfix::WritePlotHeader(const std::string& name, int nstep, Real dt, Real time) const
{
   bool is_checkpoint = 0;
   WriteHeader(name, nstep, dt, time, is_checkpoint);
}

void
mfix::WriteCheckHeader(const std::string& name, int nstep, Real dt, Real time) const
{
   bool is_checkpoint = 1;
   WriteHeader(name, nstep, dt, time, is_checkpoint);
}

void
mfix::WriteHeader(const std::string& name, int nstep, Real dt, Real time, bool is_checkpoint) const
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

      // Geometry
      for (int i = 0; i < BL_SPACEDIM; ++i)
            HeaderFile << Geometry::ProbLo(i) << ' ';
      HeaderFile << '\n';

      for (int i = 0; i < BL_SPACEDIM; ++i)
         HeaderFile << Geometry::ProbHi(i) << ' ';
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
mfix::WriteCheckPointFile(std::string& check_file, int nstep, Real dt, Real time ) const
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

    if (solve_fluid)
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
              VisMF::Write( *((*chkscalarVars[i])[lev]),
                amrex::MultiFabFileFullPrefix(lev, checkpointname,
                      level_prefix, chkscaVarsName[i]));
          }
       }
    }

    if ( solve_dem )
    {
        Vector<std::string> real_comp_names;
        Vector<std::string>  int_comp_names;
        real_comp_names.push_back("radius");
        real_comp_names.push_back("volume");
        real_comp_names.push_back("mass");
        real_comp_names.push_back("density");
        real_comp_names.push_back("omoi");
        real_comp_names.push_back("velx");
        real_comp_names.push_back("vely");
        real_comp_names.push_back("velz");
        real_comp_names.push_back("omegax");
        real_comp_names.push_back("omegay");
        real_comp_names.push_back("omegaz");
        real_comp_names.push_back("dragx");
        real_comp_names.push_back("dragy");
        real_comp_names.push_back("dragz");
         int_comp_names.push_back("phase");
         int_comp_names.push_back("state");

       bool is_checkpoint = true;
       pc -> Checkpoint(checkpointname, "particles", is_checkpoint, real_comp_names, int_comp_names);
    }


    if (solve_dem)
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

        // Also save the paramters necessary to re-buid the LSFactory
        int levelset_params[] = { levelset__refinement,
                                  levelset__pad,
                                  levelset__eb_refinement,
                                  levelset__eb_pad         };

        std::ofstream param_file;
        std::stringstream param_file_name;
        param_file_name << checkpointname << "/LSFactory_params";
        param_file.open(param_file_name.str());
        amrex::writeIntData(levelset_params, 4, param_file);
   }
}

void
mfix::Restart (std::string& restart_file, int *nstep, Real *dt, Real *time,
               IntVect& Nrep)
{
    if (ooo_debug) amrex::Print() << "Restart" << std::endl;
    BL_PROFILE("mfix::Restart()");

    amrex::Print() << "  Restarting from checkpoint " << restart_file << std::endl;

    if (Nrep != IntVect::TheUnitVector()) {
        amrex::Print() << "  Replication " << Nrep << std::endl;

        // Since a replication has taken place, the level-set function needs to
        // be re-computed:
        levelset__restart = false;

        amrex::Print() << "ATTN: Due to replication, level-set will be re-calculated."
                       << std::endl;
    }

    Real prob_lo[BL_SPACEDIM];
    Real prob_hi[BL_SPACEDIM];


    /***************************************************************************
     * Load header: set up problem domain (including BoxArray)                 *
     *              load particle data                                         *
     *              allocate mfix memory (mfix::AllocateArrays)                *
     ***************************************************************************/

    {
      std::string File(restart_file + "/Header");

      VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

      Vector<char> fileCharPtr;
      ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
      std::string fileCharPtrString(fileCharPtr.dataPtr());
      std::istringstream is(fileCharPtrString, std::istringstream::in);

      std::string line, word;

      std::getline(is, line);

      int  nlevs;
      int  int_tmp;
      Real real_tmp;

      is >> nlevs;
      GotoNextLine(is);
      // finest_level = nlevs-1;

      // Time stepping controls
      is >> int_tmp;
      *nstep = int_tmp;
      GotoNextLine(is);

      is >> real_tmp;
      *dt = real_tmp;
      GotoNextLine(is);

      is >> real_tmp;
      *time = real_tmp;
      GotoNextLine(is);

        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while (lis >> word) {
               prob_lo[i++] = std::stod(word);
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while (lis >> word) {
               prob_hi[i++] = std::stod(word);
            }
        }


        if (Nrep != IntVect::TheUnitVector())
        {
           for (int d = 0; d < BL_SPACEDIM; d++)
           {
              prob_lo[d] = Nrep[d]*prob_lo[d];
              prob_hi[d] = Nrep[d]*prob_hi[d];
           }
        }
        Geometry::ProbDomain(RealBox(prob_lo,prob_hi));

        for (int lev = 0; lev < nlevs; ++lev) {

            BoxArray orig_ba,ba;
            orig_ba.readFrom(is);
            GotoNextLine(is);

            Box orig_domain(orig_ba.minimalBox());

            if (Nrep != IntVect::TheUnitVector())
            {
               amrex::Print() << " OLD BA had " << orig_ba.size()  << " GRIDS " << std::endl;
               amrex::Print() << " OLD Domain" << orig_domain      << std::endl;
            }

            // Particle data is loaded into the MFIXParticleContainer's base
            // class using amrex::NeighborParticleContainer::Restart

            if ( solve_dem && lev == 0)
              pc->Restart(restart_file, "particles");

            amrex::Print() << "  Finished reading particle data" << std::endl;

            BoxList bl;
            for (int nb = 0; nb < orig_ba.size(); nb++) {
             for (int k = 0; k < Nrep[2]; k++) {
                 for (int j = 0; j < Nrep[1]; j++) {
                   for (int i = 0; i < Nrep[0]; i++) {
                      Box b(orig_ba[nb]);
                      IntVect lo(b.smallEnd());
                      IntVect hi(b.bigEnd());
                      IntVect shift_vec(i*orig_domain.length(0),
                                        j*orig_domain.length(1),
                                        k*orig_domain.length(2));
                      b.shift(shift_vec);
                      bl.push_back(b);
                   }
                 }
               }
            }
            ba.define(bl);

            if (Nrep != IntVect::TheUnitVector())
            {
               Box new_domain(ba.minimalBox());
               geom[lev].Domain(new_domain);

               DistributionMapping dm { ba, ParallelDescriptor::NProcs() };
               ReMakeNewLevelFromScratch(lev,ba,dm);
            }

            // This is needed before initializing level MultiFabs: ebfactories
            // should not change after the eb-dependent MultiFabs are allocated.
            make_eb_geometry();
            make_eb_factories();

            // Allocate the fluid data, NOTE: this depends on the ebfactories.
            if (solve_fluid) AllocateArrays(lev);
        }
    }

    amrex::Print() << "  Finished reading header" << std::endl;

    /***************************************************************************
     * Load fluid data                                                         *
     ***************************************************************************/

    if (solve_fluid)
    {
       // Load the field data
       for (int lev = 0, nlevs=finestLevel()+1; lev < nlevs; ++lev)
       {
          // Read velocity and pressure gradients
          MultiFab mf_vel;
          VisMF::Read(mf_vel, amrex::MultiFabFileFullPrefix(lev, restart_file, level_prefix, "u_g"));

          MultiFab mf_gp;
          VisMF::Read(mf_gp, amrex::MultiFabFileFullPrefix(lev, restart_file, level_prefix, "gpx"));

          if (Nrep == IntVect::TheUnitVector())
          {
              // Simply copy mf_vel into vel_g, mf_gp into gp
              vel_g[lev] -> copy(mf_vel, 0, 0, 3, 0, 0);
                 gp[lev] -> copy(mf_gp , 0, 0, 3, 0, 0);

          } else {

               if (mf_vel.boxArray().size() > 1)
                   amrex::Abort("Replication only works if one initial grid");

               mf_vel.FillBoundary(geom[lev].periodicity());
                mf_gp.FillBoundary(geom[lev].periodicity());

               FArrayBox single_fab_vel(mf_vel.boxArray()[0],3);
               mf_vel.copyTo(single_fab_vel);

               FArrayBox single_fab_gp ( mf_gp.boxArray()[0],3);
               mf_gp.copyTo(single_fab_gp);

              // Copy and replicate mf into velocity
              for (MFIter mfi(*vel_g[lev]); mfi.isValid(); ++mfi)
              {
                int ib = mfi.index();
                (*vel_g[lev])[ib].copy(single_fab_vel,single_fab_vel.box(),0,mfi.validbox(),0,3);
                (   *gp[lev])[ib].copy(single_fab_gp , single_fab_gp.box(),0,mfi.validbox(),0,3);
              }
          }

       // Read scalar variables
       for (int i = 0; i < chkscalarVars.size(); i++ )
       {
           Print() << "Working on: " << chkscaVarsName[i] << std::endl;
           if (Nrep != IntVect::TheUnitVector())
           {
               if ((chkscaVarsName[i] == "level_sets")
                   || (chkscaVarsName[i] == "implicit_functions"))
               {
                   Print() << "Skipping!" << std::endl;
                   continue;
               }
           }

          // If we have created the walls using the domain boundary conditions and not
          //    by creating them from implicit functions, then the implicit_functions mf
          //    will be empty.  We don't want to fail when reading so we allow the code
          //    to read it in an empty multifab just for this one.
          int allow_empty_mf = 0;
          if (chkscaVarsName[i] == "implicit_functions") allow_empty_mf = 1;

          MultiFab mf;
          VisMF::Read(mf,
                      amrex::MultiFabFileFullPrefix(lev,
                                                    restart_file, level_prefix,
                                                    chkscaVarsName[i]),
                                                    nullptr,
                                                    ParallelDescriptor::IOProcessorNumber(),
                                                    allow_empty_mf
              );

          if (Nrep == IntVect::TheUnitVector()) {

              amrex::Print() << "  - loading scalar data: " << chkscaVarsName[i] << std::endl;

             // Copy from the mf we used to read in to the mf we will use going forward
             (*chkscalarVars[i])[lev]->copy(mf, 0, 0, 1, 0, 0);


          } else {

             if (mf.boxArray().size() > 1)
                 amrex::Abort("Replication only works if one initial grid");

             mf.FillBoundary(geom[lev].periodicity());

             FArrayBox single_fab(mf.boxArray()[0],1);
             mf.copyTo(single_fab);

              // Copy and replicate mf into chkscalarVars
              for (MFIter mfi( *(*chkscalarVars[i])[lev] ); mfi.isValid(); ++mfi) {
                  int ib = mfi.index();
                  (*(*chkscalarVars[i])[lev])[ib].copy(single_fab, single_fab.box(), 0, mfi.validbox(), 0, 1);
              }
          }
        }
       }
       amrex::Print() << "  Finished reading fluid data" << std::endl;
    }

    // Make sure that the particle BoxArray is the same as the mesh data -- we can
    //      create a dual grid decomposition in the regrid operation
    if ( solve_dem)
    {
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
          pc->SetParticleBoxArray       (lev, grids[lev]);
          pc->SetParticleDistributionMap(lev,  dmap[lev]);
        }
        pc->Redistribute();
    }

    int lev = 0;
    if (Nrep == IntVect::TheUnitVector())
    {
        // We need to do this on restart regardless of whether we replicate
        pc->Redistribute();
    } else {
       // This call to Replicate adds the new particles, then calls Redistribute()
       pc->Replicate(Nrep, geom[lev], dmap[lev], grids[lev]);
    }



   /****************************************************************************
    * Load level set data from checkpoint file                                 *
    *                                                                          *
    * Since the level-set data could be defined on a different BoxArray        *
    * (compared to the rest of the checkpoint data) => the level-set data is   *
    * stored in seperate ls_raw MultiFab.                                      *
    ****************************************************************************/
    if (levelset__restart) {
        // Load level-set Multifab
        std::stringstream ls_data_path;
        ls_data_path << restart_file << "/ls_raw";

        MultiFab ls_mf;
        VisMF::Read(ls_mf, ls_data_path.str());

        // Load LSFactory parameters: => in case the user has changed the inputs
        int levelset_params[4] = { levelset__refinement,
                                   levelset__pad,
                                   levelset__eb_refinement,
                                   levelset__eb_pad         };

        std::ifstream param_file;
        std::stringstream param_file_name;
        param_file_name << restart_file << "/LSFactory_params";
        param_file.open(param_file_name.str());

        amrex::readIntData(levelset_params, 4, param_file, FPC::NativeIntDescriptor());
        int ls_ref = levelset_params[0], ls_pad = levelset_params[1],
            eb_ref = levelset_params[2], eb_pad = levelset_params[3];

        amrex::Print() << "     + Loaded level-set parameters:" << std::endl
                       << "       ref = " << ls_ref << "    pad = " << ls_pad
                       << "    eb_ref = " << eb_ref << " eb_pad = " << eb_pad
                       << std::endl;

        // Inform the user if the checkpoint parameters do not match those in the
        // inputs file. The checkpoint inputs overwrite the inputs file.
        if(ls_ref != levelset__refinement)
            amrex::Print() << "     * Overwrote levelset__refinement = " << levelset__refinement
                           << " -> " << ls_ref << std::endl;
        if(ls_pad != levelset__pad)
            amrex::Print() << "     * Overwrote levelset__pad = " << levelset__pad
                           << " -> " << ls_pad << std::endl;
        if(eb_ref != levelset__eb_refinement)
            amrex::Print() << "     * Overwrote levelset__eb_refinement = " << levelset__eb_refinement
                           << " -> " << eb_ref << std::endl;
        if(eb_pad != levelset__eb_pad)
            amrex::Print() << "     * Overwrote levelset__eb_pad = " << levelset__eb_pad
                           << " -> " << eb_pad << std::endl;

        // TODO: load level-set data from checkpoint file
        // level_set->set_data(ls_mf);
    } else {
        fill_eb_levelsets();
    }

    if (solve_fluid)
    {
        ep_g[lev]->FillBoundary(geom[lev].periodicity());
       ep_go[lev]->FillBoundary(geom[lev].periodicity());

        ro_g[lev]->FillBoundary(geom[lev].periodicity());
       ro_go[lev]->FillBoundary(geom[lev].periodicity());

         mu_g[lev]->FillBoundary(geom[lev].periodicity());

       // Fill the bc's just in case
        vel_g[lev]->FillBoundary(geom[lev].periodicity());
       vel_go[lev]->FillBoundary(geom[lev].periodicity());
    }

    // used in load balancing
    if (load_balance_type == "KnapSack") {
        if (solve_dem)
        {
           particle_cost[lev].reset(new MultiFab(pc->ParticleBoxArray(lev),
                                                 pc->ParticleDistributionMap(lev), 1, 0));
           particle_cost[lev]->setVal(0.0);
        }

        if (solve_fluid)
        {
           fluid_cost[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 0));
           fluid_cost[lev]->setVal(0.0);
        }
    }
    amrex::Print() << "  Done with mfix::Restart " << std::endl;
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
       std::string PrettyLine = "===============================================================================\n";

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

    const std::string& par_filename = amrex::Concatenate(par_ascii_file,nstep);

    pc -> WriteAsciiFile(par_filename);
}


void
mfix::WriteAverageRegions ( std::string& avg_file, int nstep, Real time ) const
{
  BL_PROFILE("mfix::WriteAverageRegions()");

  for (int lev = 0; lev < nlev; lev++)
    {

      ComputeAverageFluidVars( lev,
                               time,
                               avg_file,
                               avg_p_g,
                               avg_ep_g,
                               avg_vel_g,
                               avg_region_x_w, avg_region_x_e,
                               avg_region_y_s, avg_region_y_n,
                               avg_region_z_b, avg_region_z_t );


      //  Compute Eulerian velocities in selected regions
        pc -> ComputeAverageVelocities ( lev,
                                         time,
                                         avg_file,
                                         avg_vel_p,
                                         avg_region_x_w, avg_region_x_e,
                                         avg_region_y_s, avg_region_y_n,
                                         avg_region_z_b, avg_region_z_t );

    }

}


void
mfix::ComputeAverageFluidVars ( const int lev,
                   const amrex::Real time, const string&  basename,
                   const Vector<int>& avg_p_g,
                   const Vector<int>& avg_ep_g,
                   const Vector<int>& avg_vel_g,
                   const Vector<Real>& avg_region_x_w, const Vector<Real>& avg_region_x_e,
                   const Vector<Real>& avg_region_y_s, const Vector<Real>& avg_region_y_n,
                   const Vector<Real>& avg_region_z_b, const Vector<Real>& avg_region_z_t ) const
{

  int nregions = avg_region_x_w.size();

  const int size_p_g   = avg_p_g.size();
  const int size_ep_g  = avg_ep_g.size();
  const int size_vel_g = avg_vel_g.size();

  const Real * dx   = geom[lev].CellSize();

  //
  // Check the regions are defined correctly
  //
  if (  ( avg_region_x_e.size() != nregions ) ||
        ( avg_region_y_s.size() != nregions ) ||
        ( avg_region_y_n.size() != nregions ) ||
        ( avg_region_z_b.size() != nregions ) ||
        ( avg_region_z_t.size() != nregions )  )
    {
      amrex::Print () << "ComputeAverageVelocities: some regions are not properly defined: skipping.";
      return;
    }

  const int var_count = 6;

  // Array to hold the data for global collection.
  vector<Real> regions_data(var_count*nregions, 0.0);

  Box domain(geom[lev].Domain());

  const amrex::MultiFab* volfrac = &(ebfactory[lev] -> getVolFrac());

  // New multiFab to hold the cell center pressure.
  std::unique_ptr<MultiFab> pg_cc(new MultiFab(
             ep_g[lev]->boxArray(), ep_g[lev]->DistributionMap(),
             ep_g[lev]->nComp(),    ep_g[lev]->nGrow(), MFInfo(), *ebfactory[lev]));

  // Create a temporary nodal pressure multifab to sum in p_g and p0_g
  MultiFab pg_nd(p_g[lev]->boxArray(), dmap[lev], 1, 0);
  pg_nd.setVal(0.);
  MultiFab::Copy(pg_nd, (* p_g[lev]), 0, 0, 1, 0);
  MultiFab::Add (pg_nd, (*p0_g[lev]), 0, 0, 1, 0);

  // Create a cell-center version of the combined pressure
  amrex::average_node_to_cellcenter(*pg_cc, 0, pg_nd, 0, 1);


  for ( int nr = 0; nr < nregions; ++nr )
    {

      // Create Real box for this region
      RealBox avg_region ( {AMREX_D_DECL(avg_region_x_w[nr],avg_region_y_s[nr],avg_region_z_b[nr])},
                           {AMREX_D_DECL(avg_region_x_e[nr],avg_region_y_n[nr],avg_region_z_t[nr])} );

      // Jump to next iteration if this averaging region is not valid
      if ( !avg_region.ok () )
        {
          amrex::Print() << "ComputeAverageVelocities: region "<< nr <<" is invalid: skipping\n";
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
      for (MFIter mfi(*ep_g[lev],false); mfi.isValid(); ++mfi)
        {

          const Box& bx  = mfi.validbox();

          // this is to check efficiently if this grid contains any eb stuff
          const EBFArrayBox&  epg_fab = static_cast<EBFArrayBox const&>((*ep_g[lev])[mfi]);
          const EBCellFlagFab&  flags = epg_fab.getEBCellFlagFab();

          RealBox box_region ( bx, Geom(lev).CellSize (), Geom(lev).ProbLo() );

          if (flags.getType(bx) != FabType::covered)
            {
              if ( box_region.intersects ( avg_region ) )
                {

                  mfix_collect_fluid(BL_TO_FORTRAN_BOX(bx),
                                     BL_TO_FORTRAN_BOX(domain),
                                     BL_TO_FORTRAN_ANYD(( *ep_g[lev])[mfi]),
                                     BL_TO_FORTRAN_ANYD((     *pg_cc)[mfi]),
                                     BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
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
