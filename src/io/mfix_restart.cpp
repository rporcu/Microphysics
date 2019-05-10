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
           if (chkscaVarsName[i] == "level_sets")
           {
               Print() << "Skipping!" << std::endl;
               continue;
           }

          MultiFab mf;
          VisMF::Read(mf,
                      amrex::MultiFabFileFullPrefix(lev,
                                                    restart_file, level_prefix,
                                                    chkscaVarsName[i]),
                                                    nullptr,
                                                    ParallelDescriptor::IOProcessorNumber());

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
    if (solve_dem)
    {
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
          pc->SetParticleBoxArray       (lev, grids[lev]);
          pc->SetParticleDistributionMap(lev,  dmap[lev]);
        }
        pc->Redistribute();

       int lev = 0;
       if (Nrep == IntVect::TheUnitVector())
       {
           // We need to do this on restart regardless of whether we replicate
           pc->Redistribute();
       } else {
          // This call to Replicate adds the new particles, then calls Redistribute()
          pc->Replicate(Nrep, geom[lev], dmap[lev], grids[lev]);
       }
    }

   /****************************************************************************
    * Load level set data from checkpoint file                                 *
    *                                                                          *
    * Since the level-set data could be defined on a different BoxArray        *
    * (compared to the rest of the checkpoint data) => the level-set data is   *
    * stored in seperate ls_raw MultiFab.                                      *
    ****************************************************************************/
    if (solve_dem)
    {
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
           if   (ls_pad != levelset__pad)
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
    }

    if (solve_fluid)
    {
        for (int lev = 0; lev <= finestLevel(); lev++)
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
    }

    // used in load balancing
    if (load_balance_type == "KnapSack") {
        if (solve_dem)
        {
            for (int lev = 0; lev <= finestLevel(); lev++)
            {
               particle_cost[lev].reset(new MultiFab(pc->ParticleBoxArray(lev),
                                                     pc->ParticleDistributionMap(lev), 1, 0));
               particle_cost[lev]->setVal(0.0);
            }
        }

        if (solve_fluid)
        {
            for (int lev = 0; lev <= finestLevel(); lev++)
            {
               fluid_cost[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 0));
               fluid_cost[lev]->setVal(0.0);
            }
        }
    }
    amrex::Print() << "  Done with mfix::Restart " << std::endl;
}

