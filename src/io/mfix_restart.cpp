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


namespace restart_aux
{
    const std::string level_prefix {"Level_"};
}

using namespace restart_aux;


void
mfix::Restart (std::string& restart_file,
               int *nstep, Real *dt,
               Real *time,
               IntVect& Nrep)
{
    if (ooo_debug) amrex::Print() << "Restart" << std::endl;
    BL_PROFILE("mfix::Restart()");

    amrex::Print() << "  Restarting from checkpoint " << restart_file << std::endl;


    if (Nrep != IntVect::TheUnitVector()) {
        amrex::Print() << "  Replication " << Nrep << std::endl;

        // Since a replication has taken place, the level-set function needs to
        // be re-computed:
        levelset_restart = false;

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

        for (int lev = 0; lev < nlevs; ++lev) {

            RealBox rb(prob_lo,prob_hi);
            Geom(lev).ProbDomain(rb);
            Geom(lev).ResetDefaultProbDomain(rb);

            BoxArray orig_ba,ba;
            orig_ba.readFrom(is);
            GotoNextLine(is);

            Box orig_domain(orig_ba.minimalBox());

            if (Nrep != IntVect::TheUnitVector())
            {
               amrex::Print() << " OLD BA had " << orig_ba.size()  << " GRIDS " << std::endl;
               amrex::Print() << " OLD Domain" << orig_domain      << std::endl;
            }

            BoxList bl;
            for (int nb = 0; nb < orig_ba.size(); nb++) {
             for (int k = 0; k < Nrep[2]; k++) {
                 for (int j = 0; j < Nrep[1]; j++) {
                   for (int i = 0; i < Nrep[0]; i++) {
                      Box b(orig_ba[nb]);
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

            // Particle data is loaded into the MFIXParticleContainer's base
            // class using amrex::NeighborParticleContainer::Restart

            if (DEM::solve && lev == 0) {

              pc->Restart(restart_file, "particles");

            } else if (PIC::solve && lev == 0) {

              const int PIC_restart_refinement = PIC::restart_refinement;

              if (PIC_restart_refinement > 1) {

                PIC_to_PIC(lev, restart_file, PIC_restart_refinement);

              } else {

                pc->Restart(restart_file, "particles");

              }
            }

            amrex::Print() << "  Finished reading particle data" << std::endl;

            if (fluid.solve) AllocateArrays(lev);
        }
    }

    amrex::Print() << "  Finished reading header" << std::endl;

    /***************************************************************************
     * Load fluid data                                                         *
     ***************************************************************************/
    if (fluid.solve)
    {
      // Load the field data
      for (int lev = 0, nlevs=finestLevel()+1; lev < nlevs; ++lev)
      {
        auto replicate_data = [] (MultiFab& dst, MultiFab& src) -> void
        {
          if (src.boxArray().size() > 1)
              amrex::Abort("Replication only works if one initial grid");

          const int ncomp = src.nComp();

          FArrayBox single_fab(src.boxArray()[0], ncomp);
          src.copyTo(single_fab);

          // Copy and replicate mf into velocity
          for (MFIter mfi(dst, false); mfi.isValid(); ++mfi)
          {
            int ib = mfi.index();

            dst[ib].copy<RunOn::Host>(single_fab, single_fab.box(), 0, mfi.validbox(), 0, ncomp);
          }
        };

        {
          // Read velocity
          auto prefix = amrex::MultiFabFileFullPrefix(lev, restart_file,
                                                      level_prefix, "u_g");

          MultiFab mf_vel;
          VisMF::Read(mf_vel, prefix);

          if (Nrep == IntVect::TheUnitVector())
          {
            // Simply copy mf_vel into vel_g, mf_gp into gp
            const int ng_to_copy = 0;

            m_leveldata[lev]->vel_g->ParallelCopy(mf_vel, 0, 0, 3, ng_to_copy, ng_to_copy);

          } else {

            mf_vel.FillBoundary(geom[lev].periodicity());

            replicate_data(*(m_leveldata[lev]->vel_g), mf_vel);
          }
        }

        {
          // Read pressure gradients
          auto prefix = amrex::MultiFabFileFullPrefix(lev, restart_file,
                                                      level_prefix, "gpx");

          MultiFab mf_gp;
          VisMF::Read(mf_gp, prefix);

          if (Nrep == IntVect::TheUnitVector())
          {
            // Simply copy mf_vel into vel_g, mf_gp into gp
            const int ng_to_copy = 0;

            m_leveldata[lev]->gp->ParallelCopy(mf_gp, 0, 0, 3, ng_to_copy, ng_to_copy);

          } else {

            mf_gp.FillBoundary(geom[lev].periodicity());

            replicate_data(*(m_leveldata[lev]->gp), mf_gp);
          }
        }

        // Read scalar variables
        ResetIOChkData();

        for (int i = 0; i < chkScalarVars.size(); i++ )
        {
          if (chkscaVarsName[i] == "level_sets") {

            amrex::Print() << "  Skipping " << chkscaVarsName[i] << std::endl;
            continue;

          } else {
            amrex::Print() << "  Loading " << chkscaVarsName[i] << std::endl;

            auto prefix = amrex::MultiFabFileFullPrefix(lev, restart_file,
                                                        level_prefix, chkscaVarsName[i]);

            MultiFab mf;
            VisMF::Read(mf, prefix, nullptr, ParallelDescriptor::IOProcessorNumber());

            if (Nrep == IntVect::TheUnitVector()) {

              // Copy from the mf we used to read in to the mf we will use going forward
              const int ng_to_copy = 0;

              (*(chkScalarVars[i][lev])).ParallelCopy(mf, 0, 0, 1, ng_to_copy, ng_to_copy);

            } else {

              mf.FillBoundary(geom[lev].periodicity());

              replicate_data(*(chkScalarVars[i][lev]), mf);
            }
          }
        }

        if (advect_enthalpy)
        {
          auto& fluid_parms = *fluid.parameters;

          for (int i = 0; i < chkTVars.size(); i++ )
          {
            if ( restart_from_cold_flow && chkscaVarsName[i] == "T_g")
            {
              amrex::Print() << "  Setting T_g to T_g0 = " << fluid.T_g0 << std::endl;
              m_leveldata[lev]->T_g->setVal(fluid.T_g0);
              continue;

            } else if ( restart_from_cold_flow && chkscaVarsName[i] == "h_g") {

              const Real h_g0 = fluid_parms.calc_h_g<RunOn::Cpu>(fluid.T_g0);

              amrex::Print() << "  Setting h_g to h_g(T_g0) = " << h_g0 << std::endl;

              const Real cp_g0 = fluid_parms.calc_cp_g<RunOn::Cpu>(fluid.T_g0);
              m_leveldata[lev]->h_g->setVal(fluid.T_g0*cp_g0);
              continue;
            }

            amrex::Print() << "  Loading " << chkTVarsName[i] << std::endl;

            auto prefix = amrex::MultiFabFileFullPrefix(lev, restart_file,
                                                        level_prefix, chkTVarsName[i]);

            MultiFab mf;
            VisMF::Read(mf, prefix, nullptr, ParallelDescriptor::IOProcessorNumber());

            if (Nrep == IntVect::TheUnitVector()) {

              // Copy from the mf we used to read in to the mf we will use
              // going forward
              const int ng_to_copy = 0;

              (*(chkTVars[i][lev])).ParallelCopy(mf, 0, 0, 1, ng_to_copy, ng_to_copy);

            } else {

              mf.FillBoundary(geom[lev].periodicity());

              replicate_data(*(chkTVars[i][lev]), mf);
            }
          }
        }

        if (solve_species)
        {
          for (int i = 0; i < chkSpeciesVars.size(); i++ )
          {
            amrex::Print() << "  Loading " << chkSpeciesVarsName[i] << std::endl;

            auto prefix = amrex::MultiFabFileFullPrefix(lev, restart_file,
                                                        level_prefix, chkSpeciesVarsName[i]);

            MultiFab mf;
            VisMF::Read(mf, prefix, nullptr, ParallelDescriptor::IOProcessorNumber());

            if (Nrep == IntVect::TheUnitVector()) {

              // Copy from the mf we used to read in to the mf we will use going forward
              const int ng_to_copy = 0;

              (*(chkSpeciesVars[i][lev])).ParallelCopy(mf, 0, 0, fluid.nspecies,
                  ng_to_copy, ng_to_copy);

            } else {

              mf.FillBoundary(geom[lev].periodicity());

              replicate_data(*(chkSpeciesVars[i][lev]), mf);
            }
          }
        }
      }

      amrex::Print() << "  Finished reading fluid data" << std::endl;
    }

    // Make sure that the particle BoxArray is the same as the mesh data -- we can
    //      create a dual grid decomposition in the regrid operation
    if (DEM::solve || PIC::solve)
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
    * stored in separate ls_raw MultiFab.                                      *
    ****************************************************************************/
    if (DEM::solve || PIC::solve)
    {
        if (levelset_restart) {
           // Load level-set Multifab
           std::stringstream ls_data_path;
           ls_data_path << restart_file << "/ls_raw";

           MultiFab ls_mf;
           VisMF::Read(ls_mf, ls_data_path.str());

           // Load LSFactory parameters: => in case the user has changed the inputs
           int levelset_params[4] = { levelset_refinement,
                                      levelset_pad,
                                      levelset_eb_refinement,
                                      levelset_eb_pad         };

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
           if(ls_ref != levelset_refinement)
               amrex::Print() << "     * Overwrote levelset_refinement = " << levelset_refinement
                              << " -> " << ls_ref << std::endl;
           if   (ls_pad != levelset_pad)
               amrex::Print() << "     * Overwrote levelset_pad = " << levelset_pad
                              << " -> " << ls_pad << std::endl;
           if(eb_ref != levelset_eb_refinement)
               amrex::Print() << "     * Overwrote levelset_eb_refinement = " << levelset_eb_refinement
                              << " -> " << eb_ref << std::endl;
           if(eb_pad != levelset_eb_pad)
               amrex::Print() << "     * Overwrote levelset_eb_pad = " << levelset_eb_pad
                              << " -> " << eb_pad << std::endl;

           // TODO: load level-set data from checkpoint file
           // level_set->set_data(ls_mf);
        } else {
           fill_eb_levelsets();
        }
    }
    if (fluid.solve)
    {
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
          m_leveldata[lev]->ep_g->FillBoundary(geom[lev].periodicity());

          m_leveldata[lev]->ro_g->FillBoundary(geom[lev].periodicity());
          m_leveldata[lev]->ro_go->FillBoundary(geom[lev].periodicity());

          if (advect_enthalpy) {
            m_leveldata[lev]->T_g->FillBoundary(geom[lev].periodicity());
            m_leveldata[lev]->h_g->FillBoundary(geom[lev].periodicity());
          }

          // Fill the bc's just in case
          m_leveldata[lev]->vel_g->FillBoundary(geom[lev].periodicity());
          m_leveldata[lev]->vel_go->FillBoundary(geom[lev].periodicity());

          m_leveldata[lev]->gp->FillBoundary(geom[lev].periodicity());

          // Fill the bc's just in case
          if (solve_species) {
            m_leveldata[lev]->X_gk->FillBoundary(geom[lev].periodicity());
          }
        }
    }

    if (load_balance_type == "KnapSack" || load_balance_type == "SFC" ||
        load_balance_type == "Greedy")
    {
      if (DEM::solve || PIC::solve) {
        for (int lev(0); lev < particle_cost.size(); ++lev) {
          if (particle_cost[lev] != nullptr)  delete particle_cost[lev];
          if (particle_proc[lev] != nullptr)  delete particle_proc[lev];
        }

        particle_cost.clear();
        particle_cost.resize(nlev, nullptr);
        particle_proc.clear();
        particle_proc.resize(nlev, nullptr);

        for (int lev = 0; lev <= finestLevel(); lev++)
        {
          particle_cost[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                            pc->ParticleDistributionMap(lev), 1, 0);
          particle_cost[lev]->setVal(0.0);

          const Real proc = static_cast<Real>(ParallelDescriptor::MyProc());
          particle_proc[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                            pc->ParticleDistributionMap(lev), 1, 0);
          particle_proc[lev]->setVal(proc);
        }
      }

      if (fluid.solve) {
        for (int lev(0); lev < fluid_cost.size(); ++lev) {
          if (fluid_cost[lev] != nullptr)  delete fluid_cost[lev];
          if (fluid_proc[lev] != nullptr)  delete fluid_proc[lev];
        }

        fluid_cost.clear();
        fluid_cost.resize(nlev, nullptr);
        fluid_proc.clear();
        fluid_proc.resize(nlev, nullptr);

        for (int lev = 0; lev <= finestLevel(); lev++)
        {
          fluid_cost[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0);
          fluid_cost[lev]->setVal(0.0);

          const Real proc = static_cast<Real>(ParallelDescriptor::MyProc());
          fluid_proc[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0);
          fluid_proc[lev]->setVal(proc);
        }
      }
    }

    amrex::Print() << "  Done with mfix::Restart " << std::endl;
}
