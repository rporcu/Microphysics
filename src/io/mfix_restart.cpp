#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_AmrCore.H>
#include <AMReX_buildInfo.H>
#include <AMReX_Geometry.H>

#include <mfix.H>
#include <mfix_fluid.H>
#include <mfix_dem.H>
#include <mfix_pic.H>


void
mfix::Restart (std::string& restart_file,
               int& nstep,
               Real& dt,
               Real& time,
               IntVect& Nrep)
{
  const std::string level_prefix("Level_");

  if (!restart_file.empty()) {

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

    RealVect prob_lo;
    RealVect prob_hi;
    IntVect n_cell;

    /***************************************************************************
     * Load header: set up problem domain (including BoxArray)                 *
     *              load particle data                                         *
     *              allocate mfix memory (mfix::AllocateArrays)                *
     ***************************************************************************/

    {
      std::string File(restart_file + "/Header");
      std::string particlesFile(restart_file + "/particles/Header");

///////////////////////////////////////////////////////////////////////////////
      int num_real_comps(0);
      int num_real_comps_read(0);

      // Read num_real_comps_read
      {
        std::ifstream ff;
        ff.open(particlesFile.c_str());

        std::string temp_str;
        ff >> temp_str;

        int temp_int;
        ff >> temp_int;

        ff >> num_real_comps_read;
        Print() << "num_real_comps = " << num_real_comps_read << "\n";

        ff.close();
      }

      num_real_comps = SoArealData::count + pc->m_runtimeRealData.count;
///////////////////////////////////////////////////////////////////////////////

      VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

      Vector<char> fileCharPtr;
      ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
      std::string fileCharPtrString(fileCharPtr.dataPtr());
      std::istringstream is(fileCharPtrString, std::istringstream::in);

      // Auxiliary strings
      std::string line, word;

      // Version information
      std::string version;

      is >> version; AMREX_ALWAYS_ASSERT(version == "Checkpoint");
      is >> version; AMREX_ALWAYS_ASSERT(version == "version:");
      is >> version;

      if (version == "1" && ParallelDescriptor::IOProcessor()) {
        Warning("Can't check whether inputs match Checkpoint file or not. Use Checkpoint file version > 1");
      }

      // Multi-level control
      std::getline(is, line);
      int loc_nlev;
      is >> loc_nlev;
      finest_level = loc_nlev-1;

      // Time stepping controls
      m_rw->GotoNextLine(is);
      is >> nstep;

      m_rw->GotoNextLine(is);
      is >> dt;

      m_rw->GotoNextLine(is);
      is >> time;

      const Real version_nb = std::stod(version);

      // Geometry controls
      if (version == "1.1") {
        m_rw->GotoNextLine(is);
        is >> prob_lo;

        m_rw->GotoNextLine(is);
        is >> prob_hi;

        m_rw->GotoNextLine(is);
        is >> n_cell;

        // Check that inputs geometry matches checkpoint geometry
        RealVect inputs_prob_lo(Geom(0).ProbLo());
        RealVect inputs_prob_hi(Geom(0).ProbHi());
        IntVect inputs_n_cell(Geom(0).Domain().size());

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(prob_lo == inputs_prob_lo, "Fix inputs: ProbLo mismatch");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(prob_hi == inputs_prob_hi, "Fix inputs: ProbHi mismatch");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_cell == inputs_n_cell, "Fix inputs: n_cell mismatch");

      } else if (version == "1") {

        m_rw->GotoNextLine(is);
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

      } else {
        amrex::Abort("Unknown CheckPoint file version");
      }

      // Replicate
      if (Nrep != IntVect::TheUnitVector()) {
         for (int dir(0); dir < AMREX_SPACEDIM; dir++) {
            prob_lo[dir] *= Nrep[dir];
            prob_hi[dir] *= Nrep[dir];

            if (version == "1.1")
              n_cell[dir] *= Nrep[dir];
         }
      }

      if (version == "1.1") {
        m_rw->GotoNextLine(is);
        Real small_volfrac(0.);
        is >> small_volfrac;

        Real inputs_small_volfrac(0.);
        ParmParse pp("eb2");
        pp.query("small_volfrac", inputs_small_volfrac);

        Real error = Math::abs(small_volfrac - inputs_small_volfrac);
        Real tolerance = std::numeric_limits<Real>::epsilon();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(error < tolerance, "Fix inputs: small_volfrac mismatch");
      }

      // BoxArray controls
      for (int lev = 0; lev < loc_nlev; ++lev) {

          BoxArray orig_ba, ba;
          orig_ba.readFrom(is);
          m_rw->GotoNextLine(is);

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
             RealBox rb(prob_lo.begin(), prob_hi.begin());
             Geom(lev).ProbDomain(rb);
             Geom(lev).ResetDefaultProbDomain(rb);

             Box new_domain(ba.minimalBox());
             geom[lev].Domain(new_domain);

             DistributionMapping dm { ba, ParallelDescriptor::NProcs() };
             ReMakeNewLevelFromScratch(lev,ba,dm);
             make_eb_geometry();
          }

          // Particle data is loaded into the MFIXParticleContainer's base
          // class using amrex::NeighborParticleContainer::Restart
          // at this step, pc is constructed from mfix.init and it may have
          // different ba and dm from the fluid grids.

          if (m_dem.solve() && lev == 0) {

            if (version_nb > (1.09 + std::numeric_limits<Real>::epsilon()) &&
                (num_real_comps > num_real_comps_read)) {
              restart_pc_from_new_chkpt(restart_file, "particles");
            } else {
              pc->Restart(restart_file, "particles");
            }

          } else if (m_pic.solve() && lev == 0) {

            const int PIC_restart_refinement = m_pic.restart_refinement();

            if (PIC_restart_refinement > 1) {

              PIC_to_PIC(lev, restart_file, PIC_restart_refinement);

            } else {

              if (version_nb > (1.09 + std::numeric_limits<Real>::epsilon()) &&
                  (num_real_comps > num_real_comps_read)) {
                restart_pc_from_new_chkpt(restart_file, "particles");
              } else {
                pc->Restart(restart_file, "particles");
              }

            }
          }

          amrex::Print() << "  Finished reading particle data" << std::endl;

          if (fluid.solve()) AllocateArrays(lev);
      }
    }

    amrex::Print() << "  Finished reading header" << std::endl;

    /***************************************************************************
     * Load fluid data                                                         *
     ***************************************************************************/
    if (fluid.solve())
    {
      // Load the field data
      for (int lev = 0, loc_nlev=finestLevel()+1; lev < loc_nlev; ++lev)
      {
        auto replicate_data = [] (MultiFab& dst, MultiFab& src) -> void
        {
          if (src.boxArray().size() > 1)
            amrex::Abort("Replication only works if one initial grid");

          const int ncomp = src.nComp();

          FArrayBox single_fab(src.boxArray()[0], ncomp, The_Pinned_Arena());
          src.copyTo(single_fab);

#ifdef AMREX_USE_GPU
          const auto nreals = single_fab.size();
          const auto nbytes = nreals*sizeof(Real);
          FArrayBox single_fab_d(single_fab.box(), single_fab.nComp());
          Gpu::htod_memcpy(single_fab_d.dataPtr(), single_fab.dataPtr(), nbytes);
#endif

          // Copy and replicate mf into velocity
          for (MFIter mfi(dst, false); mfi.isValid(); ++mfi)
          {
            int ib = mfi.index();

#ifdef AMREX_USE_GPU
            dst[ib].copy<RunOn::Gpu>(single_fab_d, single_fab_d.box(), 0, mfi.validbox(), 0, ncomp);
#else
            dst[ib].copy<RunOn::Host>(single_fab, single_fab.box(), 0, mfi.validbox(), 0, ncomp);
#endif
          }
        };

        {
          // Read velocity
          auto prefix = amrex::MultiFabFileFullPrefix(lev, restart_file,
                                                      level_prefix, "u_g");

          MultiFab mf_vel(The_Pinned_Arena());
          VisMF::Read(mf_vel, prefix);

          if (Nrep == IntVect::TheUnitVector())
          {
            // Simply copy mf_vel into vel_g, mf_gp into gp
            const int ng_to_copy = 0;

            m_leveldata[lev]->vel_g->ParallelCopy(mf_vel, 0, 0, 3, ng_to_copy, ng_to_copy);

          } else {
            replicate_data(*(m_leveldata[lev]->vel_g), mf_vel);
          }
        }

        {
          // Read pressure gradients
          auto prefix = amrex::MultiFabFileFullPrefix(lev, restart_file,
                                                      level_prefix, "gpx");

          MultiFab mf_gp(The_Pinned_Arena());
          VisMF::Read(mf_gp, prefix);

          if (Nrep == IntVect::TheUnitVector())
          {
            // Simply copy mf_vel into vel_g, mf_gp into gp
            const int ng_to_copy = 0;

            m_leveldata[lev]->gp->ParallelCopy(mf_gp, 0, 0, 3, ng_to_copy, ng_to_copy);

          } else {
            replicate_data(*(m_leveldata[lev]->gp), mf_gp);
          }
        }

        // Read scalar variables
        m_rw->ResetIOChkData();

        for (int i = 0; i < m_rw->chkScalarVars.size(); i++ )
        {
          if (m_rw->chkscaVarsName[i] == "level_sets") {

            amrex::Print() << "  Skipping " << m_rw->chkscaVarsName[i] << std::endl;
            continue;

          } else {
            amrex::Print() << "  Loading " << m_rw->chkscaVarsName[i] << std::endl;

            auto prefix = amrex::MultiFabFileFullPrefix(lev, restart_file,
                                                        level_prefix, 
                                                        m_rw->chkscaVarsName[i]);

            MultiFab mf(The_Pinned_Arena());
            VisMF::Read(mf, prefix);

            if (Nrep == IntVect::TheUnitVector()) {

              // Copy from the mf we used to read in to the mf we will use going forward
              const int ng_to_copy = 0;

              (*(m_rw->chkScalarVars[i][lev])).ParallelCopy(mf, 0, 0, 1, ng_to_copy, ng_to_copy);

            } else {

              replicate_data(*(m_rw->chkScalarVars[i][lev]), mf);
            }
          }
        }

        if (fluid.solve_enthalpy())
        {
          const auto& fluid_parms = fluid.parameters();

          for (int i = 0; i < m_rw->chkTVars.size(); i++ )
          {
            if (restart_from_cold_flow && m_rw->chkscaVarsName[i] == "T_g")
            {
              amrex::Print() << "  Setting T_g to T_g0 = " << fluid.cold_flow_temperature() << std::endl;
              m_leveldata[lev]->T_g->setVal(fluid.cold_flow_temperature());
              continue;

            } else if (restart_from_cold_flow && m_rw->chkscaVarsName[i] == "h_g") {

              const Real h_g0 = fluid_parms.calc_h_g<RunOn::Host>(fluid.cold_flow_temperature());

              amrex::Print() << "  Setting h_g to h_g(T_g0) = " << h_g0 << std::endl;

              const Real cp_g0 = fluid_parms.calc_cp_g<RunOn::Host>(fluid.cold_flow_temperature());
              m_leveldata[lev]->h_g->setVal(fluid.cold_flow_temperature()*cp_g0);
              continue;
            }

            amrex::Print() << "  Loading " << m_rw->chkTVarsName[i] << std::endl;

            auto prefix = amrex::MultiFabFileFullPrefix(lev, restart_file,
                                                        level_prefix,
                                                        m_rw->chkTVarsName[i]);

            MultiFab mf(The_Pinned_Arena());
            VisMF::Read(mf, prefix);

            if (Nrep == IntVect::TheUnitVector()) {

              // Copy from the mf we used to read in to the mf we will use
              // going forward
              const int ng_to_copy = 0;

              (*(m_rw->chkTVars[i][lev])).ParallelCopy(mf, 0, 0, 1, ng_to_copy, ng_to_copy);

            } else {

              replicate_data(*(m_rw->chkTVars[i][lev]), mf);
            }
          }
        }

        if (fluid.solve_species())
        {
          for (int i = 0; i < m_rw->chkSpeciesVars.size(); i++ )
          {
            amrex::Print() << "  Loading " << m_rw->chkSpeciesVarsName[i] << std::endl;

            auto prefix = amrex::MultiFabFileFullPrefix(lev, restart_file,
                                                        level_prefix,
                                                        m_rw->chkSpeciesVarsName[i]);

            MultiFab mf(The_Pinned_Arena());
            VisMF::Read(mf, prefix);

            if (Nrep == IntVect::TheUnitVector()) {

              // Copy from the mf we used to read in to the mf we will use going forward
              const int ng_to_copy = 0;

              (*(m_rw->chkSpeciesVars[i][lev])).ParallelCopy(mf, 0, 0, fluid.nspecies(),
                                                               ng_to_copy, ng_to_copy);

            } else {

              replicate_data(*(m_rw->chkSpeciesVars[i][lev]), mf);
            }
          }
        }
      }

      amrex::Print() << "  Finished reading fluid data" << std::endl;
    }

    // Make sure that the particle BoxArray is the same as the mesh data -- we can
    //      create a dual grid decomposition in the regrid operation
    if (m_dem.solve() || m_pic.solve())
    {
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
          pc->SetParticleBoxArray       (lev, grids[lev]);
          pc->SetParticleDistributionMap(lev,  dmap[lev]);
          pc->SetParticleGeometry       (lev,  geom[lev]);
          pc->setParticleFluidGridMap   (Vector<int>());
        }
        pc->Redistribute();

        int lev = 0;
        if (Nrep != IntVect::TheUnitVector())
          pc->Replicate(Nrep, geom[lev], dmap[lev], grids[lev]);

        // load level set parameters from checkpoint
        if (levelset_restart) {
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
            if (ls_ref != levelset_refinement)
                amrex::Print() << "     * Overwrote levelset_refinement = " << levelset_refinement
                               << " -> " << ls_ref << std::endl;
            if (ls_pad != levelset_pad)
                amrex::Print() << "     * Overwrote levelset_pad = " << levelset_pad
                               << " -> " << ls_pad << std::endl;
            if (eb_ref != levelset_eb_refinement)
                amrex::Print() << "     * Overwrote levelset_eb_refinement = " << levelset_eb_refinement
                               << " -> " << eb_ref << std::endl;
            if (eb_pad != levelset_eb_pad)
                amrex::Print() << "     * Overwrote levelset_eb_pad = " << levelset_eb_pad
                               << " -> " << eb_pad << std::endl;
        }
    }

  }

  {
    // make fluid and particle factories
    // This uses particles' ba and dm so it needs to be after setting up pc.
    make_eb_factories();

    if (m_dem.solve() || m_pic.solve())
      fill_eb_levelsets();


    if (fluid.solve())
    {
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
          m_leveldata[lev]->ep_g->FillBoundary(geom[lev].periodicity());

          m_leveldata[lev]->ro_g->FillBoundary(geom[lev].periodicity());
          m_leveldata[lev]->ro_go->FillBoundary(geom[lev].periodicity());

          if (fluid.solve_enthalpy()) {
            m_leveldata[lev]->T_g->FillBoundary(geom[lev].periodicity());
            m_leveldata[lev]->T_go->FillBoundary(geom[lev].periodicity());
            m_leveldata[lev]->h_g->FillBoundary(geom[lev].periodicity());
            m_leveldata[lev]->h_go->FillBoundary(geom[lev].periodicity());
          }

          // Fill the bc's just in case
          m_leveldata[lev]->vel_g->FillBoundary(geom[lev].periodicity());
          m_leveldata[lev]->vel_go->FillBoundary(geom[lev].periodicity());

          m_leveldata[lev]->gp->FillBoundary(geom[lev].periodicity());

          // Fill the bc's just in case
          if (fluid.solve_species()) {
            m_leveldata[lev]->X_gk->FillBoundary(geom[lev].periodicity());
          }
        }
    }

    // setup the multifabs to track cost and display the rank
    if (load_balance_type == "KnapSack" || load_balance_type == "SFC" ||
        load_balance_type == "Greedy")
    {
      if (m_dem.solve() || m_pic.solve()) {
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

          // Initialize particles cost
          if (load_balance_type == "KnapSack" && knapsack_weight_type == "NumParticles") {
            for (MFIXParticleContainer::MFIXParIter pti(*pc,lev); pti.isValid(); ++pti)
              pc->UpdateCost(particle_cost[lev], pti, knapsack_weight_type, 0.);
          } else {
            particle_cost[lev]->setVal(0.);
          }

          const Real proc = static_cast<Real>(ParallelDescriptor::MyProc());
          particle_proc[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                            pc->ParticleDistributionMap(lev), 1, 0);
          particle_proc[lev]->setVal(proc);

        }
      }

      if (fluid.solve()) {
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

  }

  amrex::Print() << "  Done with mfix::Restart " << std::endl;
}


void mfix::
restart_pc_from_new_chkpt (const std::string& restart_file,
                           const std::string& name)
{
  const int nspecies_s = solids.nspecies();
  const int solve_enthalpy = solids.solve_enthalpy();
  const int solve_reactions = reactions.solve();
  const int nreactions = reactions.nreactions();
  const int solve_pic = m_pic.solve();
  const int solve_cg = m_dem.cg_dem();

  // New SoArealData
  struct new_SoArealData { enum {radius = 0, mass, velx, vely, velz, omegax,
                                 omegay, omegaz, dragcoeff, dragx, dragy, dragz,
                                 count}; };

  // New SoAintData
  struct new_SoAintData { enum {phase = 0, state,
#if MFIX_POLYDISPERSE
                                ptype,
#endif
                                count}; };

  //
  constexpr int NArrayReal = new_SoArealData::count;
  constexpr int NArrayInt = new_SoAintData::count;

  using PC = amrex::NeighborParticleContainer<0,0,NArrayReal,NArrayInt>;
  using PCIter = amrex::ParIter<0,0,NArrayReal,NArrayInt>;

  PC* new_pc = new PC(GetParGDB(), 1);

  const int new_idx_X_sn = 0;
  const int new_idx_cp_s = new_idx_X_sn + nspecies_s;
  const int new_idx_temperature = new_idx_cp_s + 1*solve_enthalpy;
  const int new_idx_convection = new_idx_temperature + 1*solve_enthalpy;
  const int new_idx_vel_txfr = new_idx_convection + 1*solve_enthalpy;
  const int new_idx_h_txfr = new_idx_vel_txfr + 3*solve_reactions;
  const int new_idx_mass_txfr = new_idx_h_txfr + 1*solve_reactions;
  const int new_idx_ep_s = new_idx_mass_txfr + nspecies_s*(nreactions>0);

  const int new_idx_statwt = new_idx_ep_s +
    1*((solve_pic > 0) || (m_run_type == RunType::PIC2DEM));

  const int new_count = new_idx_statwt +
    1*((solve_pic > 0) || (solve_cg > 0) || (m_run_type == RunType::PIC2DEM));

  // Check if this is needed
  for (int n(0); n < new_count; ++n) {
    new_pc->AddRealComp(true);

    if (n == new_idx_temperature ||
        n == new_idx_convection)
      new_pc->setRealCommComp((new_SoArealData::count+2)+n, true);
    else
      new_pc->setRealCommComp((new_SoArealData::count+2)+n, false);
  }

  new_pc->Restart(restart_file, name);

  // Check if this is needed
  new_pc->Redistribute();

  const int idx_X_sn = pc->m_runtimeRealData.X_sn;
  const int idx_vel_txfr = pc->m_runtimeRealData.vel_txfr;
  const int idx_h_txfr = pc->m_runtimeRealData.h_txfr;
  const int idx_mass_txfr = pc->m_runtimeRealData.mass_txfr;

  for (int lev(0); lev < nlev; ++lev) {
    for (PCIter pti(*new_pc,lev); pti.isValid(); ++pti) {

      const int np = pti.numParticles();

      MFIXParticleContainer::PairIndex index(pti.index(), pti.LocalTileIndex());
      auto& new_particles = new_pc->GetParticles(lev);
      auto& new_ptile = new_particles[index];
      auto& new_aos = new_ptile.GetArrayOfStructs();
      auto new_pstruct = new_aos().dataPtr();
      auto& new_soa = new_ptile.GetStructOfArrays();
      auto new_realarray = new_soa.realarray();
      auto new_intarray = new_soa.intarray();
      auto new_tile_data = new_ptile.getParticleTileData();

      auto& particletile = pc->DefineAndReturnParticleTile(lev, pti);
      particletile.resize(np);
      auto& particles  = pc->GetParticles(lev);
      auto& ptile = particles[index];
      auto& aos = particletile.GetArrayOfStructs();
      auto pstruct = aos().dataPtr();
      auto& soa = particletile.GetStructOfArrays();
      auto realarray = soa.realarray();
      auto intarray = soa.intarray();
      auto tile_data = ptile.getParticleTileData();

      const int run_type = m_run_type;

      ParallelFor(np, [pstruct,new_pstruct,new_realarray,realarray,
          new_intarray,intarray,new_tile_data,tile_data,new_idx_X_sn,
          new_idx_vel_txfr,new_idx_h_txfr,new_idx_mass_txfr,idx_X_sn,new_idx_cp_s,
          new_idx_temperature,new_idx_convection,idx_vel_txfr,idx_h_txfr,
          idx_mass_txfr,new_idx_ep_s,new_idx_statwt,nspecies_s,solve_enthalpy,
          solve_reactions,solve_pic,solve_cg,run_type]
        AMREX_GPU_DEVICE (int i) noexcept
      {
        auto& part = pstruct[i];
        auto& new_part = new_pstruct[i];

        part.pos(0) = new_part.pos(0);
        part.pos(1) = new_part.pos(1);
        part.pos(2) = new_part.pos(2);

        part.id() = new_part.id();
        part.cpu() = new_part.cpu();

        intarray[SoAintData::phase][i] = new_intarray[new_SoAintData::phase][i];
        intarray[SoAintData::state][i] = new_intarray[new_SoAintData::state][i];
#if MFIX_POLYDISPERSE
        intarray[SoAintData::ptype][i] = new_intarray[new_SoAintData::ptype][i];
#endif

        const Real radius = new_realarray[new_SoArealData::radius][i];
        const Real mass = new_realarray[new_SoArealData::mass][i];

        realarray[SoArealData::radius][i] = new_realarray[new_SoArealData::radius][i];
        realarray[SoArealData::volume][i] = (4.0/3.0)*M_PI*(radius*radius*radius);
        realarray[SoArealData::mass][i] = new_realarray[new_SoArealData::mass][i];
        realarray[SoArealData::density][i] = mass / realarray[SoArealData::volume][i];
        realarray[SoArealData::oneOverI][i] = 2.5/(mass*(radius*radius));
        realarray[SoArealData::velx][i] = new_realarray[new_SoArealData::velx][i];
        realarray[SoArealData::vely][i] = new_realarray[new_SoArealData::vely][i];
        realarray[SoArealData::velz][i] = new_realarray[new_SoArealData::velz][i];
        realarray[SoArealData::omegax][i] = new_realarray[new_SoArealData::omegax][i];
        realarray[SoArealData::omegay][i] = new_realarray[new_SoArealData::omegay][i];
        realarray[SoArealData::omegaz][i] = new_realarray[new_SoArealData::omegaz][i];
        realarray[SoArealData::statwt][i] = 1.;
        realarray[SoArealData::dragcoeff][i] = new_realarray[new_SoArealData::dragcoeff][i];
        realarray[SoArealData::dragx][i] = new_realarray[new_SoArealData::dragx][i];
        realarray[SoArealData::dragy][i] = new_realarray[new_SoArealData::dragy][i];
        realarray[SoArealData::dragz][i] = new_realarray[new_SoArealData::dragz][i];

        for (int n_s(0); n_s < nspecies_s; ++n_s)
          tile_data.m_runtime_rdata[idx_X_sn+n_s][i] = new_tile_data.m_runtime_rdata[new_idx_X_sn+n_s][i];

        if (solve_enthalpy) {
          realarray[SoArealData::cp_s][i] = new_tile_data.m_runtime_rdata[new_idx_cp_s][i];
          realarray[SoArealData::temperature][i] = new_tile_data.m_runtime_rdata[new_idx_temperature][i];
          realarray[SoArealData::convection][i] = new_tile_data.m_runtime_rdata[new_idx_convection][i];
        } else {
          realarray[SoArealData::cp_s][i] = 0.;
          realarray[SoArealData::temperature][i] = 0.;
          realarray[SoArealData::convection][i] = 0.;
        }

        if (solve_reactions) {
          tile_data.m_runtime_rdata[idx_vel_txfr+0][i] = new_tile_data.m_runtime_rdata[new_idx_vel_txfr+0][i];
          tile_data.m_runtime_rdata[idx_vel_txfr+1][i] = new_tile_data.m_runtime_rdata[new_idx_vel_txfr+1][i];
          tile_data.m_runtime_rdata[idx_vel_txfr+2][i] = new_tile_data.m_runtime_rdata[new_idx_vel_txfr+2][i];

          tile_data.m_runtime_rdata[idx_h_txfr][i] = new_tile_data.m_runtime_rdata[new_idx_h_txfr][i];

          for (int n_s(0); n_s < nspecies_s; ++n_s)
            tile_data.m_runtime_rdata[idx_mass_txfr+n_s][i] = new_tile_data.m_runtime_rdata[new_idx_mass_txfr+n_s][i];
        }

        if (solve_pic || (run_type == RunType::PIC2DEM)) {
          realarray[SoArealData::oneOverI][i] = new_tile_data.m_runtime_rdata[new_idx_ep_s][i];
        }

        if (solve_pic || solve_cg || (run_type == RunType::PIC2DEM))
          realarray[SoArealData::statwt][i] = new_tile_data.m_runtime_rdata[new_idx_statwt][i];
      });
    }
  }

  delete new_pc;
}
