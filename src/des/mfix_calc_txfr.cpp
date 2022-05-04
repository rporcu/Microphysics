#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_shepard_K.H>
#include <mfix_filcc.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>

#include <mfix_mf_helpers.H>
#include <mfix_dem_parms.H>
#include <mfix_des_rrates_K.H>
#include <mfix_algorithm.H>

void
mfix::mfix_calc_txfr_fluid (Vector< MultiFab* > const& txfr_out,
                            Vector< MultiFab* > const& chem_txfr_out,
                            Vector< MultiFab* > const& ep_g_in,
                            Vector< MultiFab* > const& ro_g_in,
                            Vector< MultiFab* > const& vel_g_in,
                            Vector< MultiFab* > const& T_g_in,
                            Vector< MultiFab* > const& X_gk_in,
                            Vector< MultiFab* > const& pressure_g_in,
                            const Real time)
{
  BL_PROFILE("mfix::mfix_calc_txfr_fluid()");

  const Real strttime = ParallelDescriptor::second();

  mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in);

  // ******************************************************************************
  // Now use the transfer coeffs of individual particles to create the
  // interphase transfer terms on the fluid
  // ******************************************************************************
  for (int lev = 0; lev < nlev; lev++) {
    txfr_out[lev]->setVal(0);

    if (reactions.solve)
      chem_txfr_out[lev]->setVal(0);
  }

  if (nlev > 2)
    amrex::Abort("For right now"
        " MFIXParticleContainer::TrilinearDepositionFluidDragForce can only"
        " handle up to 2 levels");

  Vector< MultiFab* > txfr_ptr(nlev, nullptr);
  Vector< MultiFab* > chem_txfr_ptr(nlev, nullptr);

  for (int lev = 0; lev < nlev; lev++) {

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    if (lev == 0 && OnSameGrids) {

      // If we are already working with the internal mf defined on the
      // particle_box_array, then we just work with this.
      txfr_ptr[lev] = txfr_out[lev];
      chem_txfr_ptr[lev] = chem_txfr_out[lev];

    } else if (lev == 0 && (!OnSameGrids)) {

      // If beta_mf is not defined on the particle_box_array, then we need
      // to make a temporary here and copy into beta_mf at the end.
      txfr_ptr[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                   pc->ParticleDistributionMap(lev),
                                   txfr_out[lev]->nComp(),
                                   txfr_out[lev]->nGrow());

      if (reactions.solve) {
        chem_txfr_ptr[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                          pc->ParticleDistributionMap(lev),
                                          chem_txfr_out[lev]->nComp(),
                                          chem_txfr_out[lev]->nGrow());
      } 

    } else {
      // If lev > 0 we make a temporary at the coarse resolution
      BoxArray ba_crse(amrex::coarsen(pc->ParticleBoxArray(lev), this->m_gdb->refRatio(0)));
      txfr_ptr[lev] = new MultiFab(ba_crse, pc->ParticleDistributionMap(lev),
                                   txfr_out[lev]->nComp(), 1);
      chem_txfr_ptr[lev] = new MultiFab(ba_crse, pc->ParticleDistributionMap(lev),
                                        chem_txfr_out[lev]->nComp(), 1);
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (txfr_ptr[lev]->nGrow() < 1)
      amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

    if (reactions.solve) {
      if (chem_txfr_ptr[lev]->nGrow() < 1)
        amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");
    }

    txfr_ptr[lev]->setVal(0.0, 0, txfr_out[lev]->nComp(), txfr_ptr[lev]->nGrow());

    if (reactions.solve)
      chem_txfr_ptr[lev]->setVal(0.0, 0, chem_txfr_out[lev]->nComp(), chem_txfr_ptr[lev]->nGrow());
  }

  const Geometry& gm = Geom(0);
  Vector< const FabArray<EBCellFlagFab>* > flags(nlev, nullptr);
  Vector< const MultiFab* > volfrac(nlev, nullptr);
  Vector< EBFArrayBoxFactory* > crse_factory(nlev, nullptr);

  Vector< MultiFab* > tmp_eps(nlev);

  for (int lev = 0; lev < nlev; lev++) {

    tmp_eps[lev] = (MFHelpers::createFrom(*txfr_ptr[lev], 0.0)).release();

    // Use level 0 to define the EB factory. If we are not on level 0
    // then create a copy of the coarse factory to use.

    if (lev == 0) {
      flags[lev]   = &(particle_ebfactory[lev]->getMultiEBCellFlagFab());
      volfrac[lev] = &(particle_ebfactory[lev]->getVolFrac());

    } else {

      Vector<int> ngrow = {1,1,1};

      crse_factory[lev] = (makeEBFabFactory(gm, txfr_ptr[lev]->boxArray(),
                                            txfr_ptr[lev]->DistributionMap(),
                                            ngrow, EBSupport::volume)).release();

      flags[lev]   = &(crse_factory[lev]->getMultiEBCellFlagFab());
      volfrac[lev] = &(crse_factory[lev]->getVolFrac());
    }

    // Deposit the interphase transfer forces to the grid
    // Drag force: (beta and beta*particle_vel)
    // Heat transfer: gamma and gamma*particle temperature
    pc->InterphaseTxfrDeposition(lev, *tmp_eps[lev], *txfr_ptr[lev],
                                 *chem_txfr_ptr[lev], volfrac[lev], flags[lev]);
  }

  {

    // The deposition occurred on level 0, thus the next few operations
    // only need to be carried out on level 0.
    int lev(0);

    // Move any volume deposited outside the domain back into the domain
    // when BC is either a pressure inlet or mass inflow.
    mfix_deposition_bcs(lev, *txfr_ptr[lev]);

    if (reactions.solve)
      mfix_deposition_bcs(lev, *chem_txfr_ptr[lev]);

    // Sum grid boundaries to capture any material that was deposited into
    // your grid from an adjacent grid.
    txfr_ptr[lev]->SumBoundary(gm.periodicity());
    txfr_ptr[lev]->setBndry(0.0);

    if (reactions.solve) {
      chem_txfr_ptr[lev]->SumBoundary(gm.periodicity());
      chem_txfr_ptr[lev]->setBndry(0.0);
    }

    // Sum grid boundaries then fill with correct ghost values.
    tmp_eps[lev]->SumBoundary(gm.periodicity());
    tmp_eps[lev]->FillBoundary(gm.periodicity());

    // Move excessive solids volume from small cells to neighboring cells.
    // Note that we don't change tmp_eps but use the redistribution of
    // particle volume to determine how to redistribute the drag forces.
    mfix_redistribute_deposition(lev, *tmp_eps[lev], *txfr_ptr[lev], volfrac[lev], flags[lev],
                                 mfix::m_max_solids_volume_fraction);

    if (reactions.solve) {
      mfix_redistribute_deposition(lev, *tmp_eps[lev], *chem_txfr_ptr[lev], volfrac[lev], flags[lev],
                                   mfix::m_max_solids_volume_fraction);
    }

    // Sum the boundaries again to recapture any solids moved across
    // grid boundaries during the redistribute
    txfr_ptr[lev]->SumBoundary(gm.periodicity());
    txfr_ptr[lev]->FillBoundary(gm.periodicity());

    if (reactions.solve) {
      chem_txfr_ptr[lev]->SumBoundary(gm.periodicity());
      chem_txfr_ptr[lev]->FillBoundary(gm.periodicity());
    }
  }

  for (int lev(0); lev < nlev; ++lev) {
    if (lev != 0 && crse_factory[lev] != nullptr)
      delete crse_factory[lev];
  }

  // This might not need to exist on all levels. Maybe only level 0.
  for (int lev(0); lev < nlev; ++lev)
    delete tmp_eps[lev];

  int  src_nghost = 1;
  int dest_nghost = 0;
  int ng_to_copy = amrex::min(src_nghost, dest_nghost);

  for (int lev = 1; lev < nlev; lev++) {
    txfr_ptr[0]->ParallelCopy(*txfr_ptr[lev], 0, 0, txfr_ptr[0]->nComp(), ng_to_copy,
        ng_to_copy, gm.periodicity(), FabArrayBase::ADD);

    if (reactions.solve) {
      chem_txfr_ptr[0]->ParallelCopy(*chem_txfr_ptr[lev], 0, 0, chem_txfr_ptr[0]->nComp(), ng_to_copy,
          ng_to_copy, gm.periodicity(), FabArrayBase::ADD);
    }
  }

  if (nlev > 1) {

    // IntVect ref_ratio(this->m_gdb->refRatio(0));

    // Now interpolate from the coarse grid to define the fine grid ep-g
    Interpolater* mapper = &cell_cons_interp;
    int lo_bc[3] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
    int hi_bc[3] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
    Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

    BndryFuncArray bfunc(mfix_aux::filcc);

    for (int lev = 1; lev < nlev; lev++) {

      PhysBCFunct<BndryFuncArray> cphysbc(Geom(lev-1), bcs, bfunc);
      PhysBCFunct<BndryFuncArray> fphysbc(Geom(lev  ), bcs, bfunc);

      txfr_out[lev]->setVal(0);
      amrex::InterpFromCoarseLevel(*txfr_out[lev], time, *txfr_ptr[lev-1],
                                   0, 0, 1, Geom(lev-1), Geom(lev),
                                   cphysbc, 0, fphysbc, 0,
                                   ref_ratio[0], mapper,
                                   bcs, 0);

      if (reactions.solve) {
        chem_txfr_out[lev]->setVal(0);
        amrex::InterpFromCoarseLevel(*chem_txfr_out[lev], time, *chem_txfr_ptr[lev-1],
                                     0, 0, 1, Geom(lev-1), Geom(lev),
                                     cphysbc, 0, fphysbc, 0,
                                     ref_ratio[0], mapper,
                                     bcs, 0);
      }
    }
  }

  // If mf_to_be_filled is not defined on the particle_box_array, then we need
  // to copy here from txfr_ptr into mf_to_be_filled. I believe that we don't
  // need any information in ghost cells so we don't copy those.

  if (txfr_ptr[0] != txfr_out[0]) {
    txfr_out[0]->ParallelCopy(*txfr_ptr[0], 0, 0, txfr_out[0]->nComp());
  }

  if (reactions.solve) {
    if (chem_txfr_ptr[0] != chem_txfr_out[0]) {
      chem_txfr_out[0]->ParallelCopy(*chem_txfr_ptr[0], 0, 0, chem_txfr_out[0]->nComp());
    }
  }

  for (int lev = 0; lev < nlev; lev++) {
    if (txfr_ptr[lev] != txfr_out[lev]) {
      delete txfr_ptr[lev];
      delete chem_txfr_ptr[lev];
    }
  }

  if (m_verbose > 1) {
    Real stoptime = ParallelDescriptor::second() - strttime;

    ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "MFIXParticleContainer::TrilinearDepositionFluidDragForce time: " << stoptime << '\n';
  }

  if(mfix::m_deposition_diffusion_coeff > 0.) {
    // Apply mean field diffusion to drag force
    diffusion_op->diffuse_drag(txfr_out, mfix::m_deposition_diffusion_coeff);
  }
}


void
mfix::mfix_calc_txfr_particle (Real time,
                               Vector< MultiFab* > const& ep_g_in,
                               Vector< MultiFab* > const& ro_g_in,
                               Vector< MultiFab* > const& vel_g_in,
                               Vector< MultiFab* > const& T_g_in,
                               Vector< MultiFab* > const& X_gk_in,
                               Vector< MultiFab* > const& pressure_g_in,
                               Vector< MultiFab* > const& gp_in)
{
  if (m_reaction_rates_type == ReactionRatesType::RRatesUser) {
    mfix_calc_txfr_particle(time, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in,
                            pressure_g_in, gp_in, HeterogeneousRatesUser());
  } else {
    amrex::Abort("Invalid Reaction Rates Type.");
  }
}


template <typename F1>
void
mfix::mfix_calc_txfr_particle (Real time,
                               Vector< MultiFab* > const& ep_g_in,
                               Vector< MultiFab* > const& ro_g_in,
                               Vector< MultiFab* > const& vel_g_in,
                               Vector< MultiFab* > const& T_g_in,
                               Vector< MultiFab* > const& X_gk_in,
                               Vector< MultiFab* > const& pressure_g_in,
                               Vector< MultiFab* > const& gp_in,
                               F1 HeterogeneousRRates)
{
  using PairIndex = MFIXParticleContainer::PairIndex;
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  BL_PROFILE("mfix::mfix_calc_txfr_particle()");

  //***************************************************************************
  // Data for chemical reactions
  //***************************************************************************
  // Solid species data
  const int nspecies_s = solids.nspecies;

  // Fluid species data
  const int nspecies_g = fluid.nspecies;

  // Reactions data
  const int nreactions = reactions.nreactions;

  // Particles SoA starting indexes for mass fractions and rate of formations
  const int idx_X_sn   = (pc->m_runtimeRealData).X_sn;
  const int idx_mass_txfr = (pc->m_runtimeRealData).mass_txfr;
  const int idx_vel_txfr = (pc->m_runtimeRealData).vel_txfr;
  const int idx_h_txfr = (pc->m_runtimeRealData).h_txfr;

  auto& fluid_parms = *fluid.parameters;
  auto& solids_parms = *solids.parameters;
  auto& reactions_parms = *reactions.parameters;

  //***************************************************************************
  // 
  //***************************************************************************

  // Extrapolate velocity Dirichlet bc's to ghost cells
  int extrap_dir_bcs = 1;

  mfix_set_velocity_bcs(time, vel_g_in, extrap_dir_bcs);

  if (fluid.solve_enthalpy) {
    mfix_set_temperature_bcs(time, T_g_in);
  }

  if (reactions.solve) {

    const int dir_bc_in = 2;
    mfix_set_epg_bcs(ep_g_in, dir_bc_in);

    mfix_set_density_bcs(time, ro_g_in);

    mfix_set_species_bcs(time, X_gk_in);
  }

  for (int lev = 0; lev < nlev; lev++) {

    Box domain(geom[lev].Domain());
    MultiFab gp_tmp;

    gp_tmp.define(grids[lev], dmap[lev], 3, 1, MFInfo(), *ebfactory[lev]);

    MultiFab::Copy(gp_tmp, *gp_in[lev], 0, 0, 3, 1);
    gp_tmp.FillBoundary(geom[lev].periodicity());

    //
    // NOTE -- it is essential that we call set_gradp_bcs after calling FillBoundary
    //         because the set_gradp_bcs call hopefully sets the ghost cells exterior
    //         to the domain from ghost cells interior to the domain
    //

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(gp_tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      set_gradp_bcs(bx, lev, gp_tmp[mfi], domain);
    }

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    // Pointer to Multifab for interpolation
    MultiFab* interp_ptr;

    // This is just a sanity check to make sure we're not using covered values
    // We can remove these lines once we're confident in the algorithm
    EB_set_covered(*vel_g_in[lev], 0, 3, 1, covered_val);
    EB_set_covered(gp_tmp, 0, 3, 1, covered_val);

    if (fluid.solve_enthalpy)
      EB_set_covered(*T_g_in[lev], 0, 1, 1, covered_val);

    if (reactions.solve) {
      EB_set_covered(*ep_g_in[lev], 0, 1, 1, covered_val);
      EB_set_covered(*ro_g_in[lev], 0, 1, 1, covered_val);
      EB_set_covered(*X_gk_in[lev], 0, fluid.nspecies, 1, covered_val);
      EB_set_covered(*pressure_g_in[lev], 0, 1, 1, covered_val);
    }

    const int interp_ng = 1;    // Only one layer needed for interpolation
    const int interp_comp = 6 +  // 3 vel_g + 3 gp
                            1*int(fluid.solve_enthalpy) +  // 1 T_g
                            (3+fluid.nspecies)*int(reactions.solve); // ep_g + ro_g + X_gk + pressure_g

    MultiFab pressure_cc(ep_g_in[lev]->boxArray(), dmap[lev], 1, interp_ng);
    pressure_cc.setVal(0.);

    if (reactions.solve) {
      MultiFab pressure_nd(pressure_g_in[lev]->boxArray(), dmap[lev], 1, interp_ng);
      pressure_nd.setVal(0.);
      EB_set_covered(pressure_nd, 0, 1, 1, covered_val);

      MultiFab::Copy(pressure_nd, *m_leveldata[lev]->thermodynamic_p_g, 0, 0, 1, interp_ng);
      MultiFab::Add (pressure_nd, *m_leveldata[lev]->p0_g, 0, 0, 1, interp_ng);

      pressure_cc.setVal(0.);
      EB_set_covered(pressure_cc, 0, 1, 1, covered_val);

      amrex::average_node_to_cellcenter(pressure_cc, 0, pressure_nd, 0, 1);
    }

    if (OnSameGrids)
    {
      // Store gas velocity and pressure gradient for interpolation
      interp_ptr = new MultiFab(grids[lev], dmap[lev], interp_comp, interp_ng, MFInfo(), *ebfactory[lev]);

      int components_count(0);

      // Copy fluid velocity
      MultiFab::Copy(*interp_ptr, *vel_g_in[lev], 0, components_count, 3, interp_ng);
      components_count += 3;

      // Copy pressure gradient
      MultiFab::Copy(*interp_ptr, gp_tmp, 0, components_count, 3, interp_ng);
      components_count += 3;

      // Copy fluid temperature
      if (fluid.solve_enthalpy) {
        MultiFab::Copy(*interp_ptr, *T_g_in[lev], 0, components_count, 1, interp_ng);
        components_count += 1;
      }

      if (reactions.solve) {
        // Copy volume fraction
        MultiFab::Copy(*interp_ptr, *ep_g_in[lev],  0, components_count, 1, interp_ng);
        components_count += 1;

        // Copy fluid density
        MultiFab::Copy(*interp_ptr, *ro_g_in[lev],  0, components_count, 1, interp_ng);
        components_count += 1;

        // Copy species mass fractions
        MultiFab::Copy(*interp_ptr, *X_gk_in[lev],  0, components_count, fluid.nspecies, interp_ng);
        components_count += fluid.nspecies;

        // Copy thermodynamic pressure
        MultiFab::Copy(*interp_ptr, pressure_cc,  0, components_count, 1, interp_ng);
        components_count += 1;
      }

      AMREX_ALWAYS_ASSERT(interp_comp == components_count);

    } else {

      const BoxArray&            pba = pc->ParticleBoxArray(lev);
      const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(pba, pdm, interp_comp, interp_ng, MFInfo(), *particle_ebfactory[lev]);

      int components_count(0);

      // Copy fluid velocity
      interp_ptr->ParallelCopy(*vel_g_in[lev], 0, components_count, 3, interp_ng, interp_ng);
      components_count += 3;

      // Copy pressure gradient
      interp_ptr->ParallelCopy(gp_tmp, 0, components_count, 3, interp_ng, interp_ng);
      components_count += 3;

      // Copy fluid temperature
      if (fluid.solve_enthalpy) {
        interp_ptr->ParallelCopy(*T_g_in[lev], 0, components_count, 1, interp_ng, interp_ng);
        components_count += 1;
      }

      if (reactions.solve) {
        // Copy volume fraction
        interp_ptr->ParallelCopy(*ep_g_in[lev],  0, components_count, 1, interp_ng, interp_ng);
        components_count += 1;

        // Copy fluid density
        interp_ptr->ParallelCopy(*ro_g_in[lev],  0, components_count, 1, interp_ng, interp_ng);
        components_count += 1;

        // Copy species mass fractions
        interp_ptr->ParallelCopy(*X_gk_in[lev],  0, components_count, fluid.nspecies, interp_ng, interp_ng);
        components_count += fluid.nspecies;

        // Copy thermodynamic pressure
        interp_ptr->ParallelCopy(pressure_cc,  0, components_count, 1, interp_ng, interp_ng);
        components_count += 1;
      }

      AMREX_ALWAYS_ASSERT(interp_comp == components_count);
    }

    // FillBoundary on interpolation MultiFab
    interp_ptr->FillBoundary(geom[lev].periodicity());

    BL_PROFILE_VAR("particle_deposition", particle_deposition);

    {
      const auto dx_array  = geom[lev].CellSizeArray();
      const auto dxi_array = geom[lev].InvCellSizeArray();
      const auto plo_array = geom[lev].ProbLoArray();

      const RealVect  dx( dx_array[0],  dx_array[1],  dx_array[2]);
      const RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
      const RealVect plo(plo_array[0], plo_array[1], plo_array[2]);

      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(interp_ptr->Factory());

      const auto cellcent = &(factory.getCentroid());
      const auto bndrycent = &(factory.getBndryCent());
      const auto areafrac = factory.getAreaFrac();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto& ptile = pc->GetParticles(lev)[index];

        //Access to added variables
        auto ptile_data = ptile.getParticleTileData();

        auto& particles = pti.GetArrayOfStructs();
        MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        auto& soa = pti.GetStructOfArrays();
        auto p_realarray = soa.realarray();

        const int np = particles.size();

        Box bx = pti.tilebox ();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  interp_fab = static_cast<EBFArrayBox const&>((*interp_ptr)[pti]);
        const EBCellFlagFab&  flags = interp_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered) {

          const auto& interp_array = interp_ptr->array(pti);

          const auto& flags_array = flags.array();

          const int grown_bx_is_regular = (flags.getType(amrex::grow(bx,1)) == FabType::regular);

          const Array4<const Real> empty_array;

          // Cell centroids
          const auto& ccent_fab = grown_bx_is_regular ? empty_array : cellcent->const_array(pti);
          // Centroid of EB
          const auto& bcent_fab = grown_bx_is_regular ? empty_array : bndrycent->const_array(pti);
          // Area fractions
          const auto& apx_fab = grown_bx_is_regular ? empty_array : areafrac[0]->const_array(pti);
          const auto& apy_fab = grown_bx_is_regular ? empty_array : areafrac[1]->const_array(pti);
          const auto& apz_fab = grown_bx_is_regular ? empty_array : areafrac[2]->const_array(pti);

          const int solve_enthalpy = fluid.solve_enthalpy;
          const int solve_reactions = reactions.solve;

          // We need this until we remove static attribute from mfix::gp0;
          const RealVect gp0_dev(gp0);

          const Real pmult = DEM::solve ? 1.0 : 0.0;

          amrex::ParallelFor(np,
              [pstruct,p_realarray,interp_array,gp0_dev,plo,dxi,pmult,dx,
               solve_enthalpy,ccent_fab,bcent_fab,apx_fab,apy_fab,apz_fab,
               flags_array,grown_bx_is_regular,interp_comp,solve_reactions,
               ptile_data,nspecies_s,nspecies_g,idx_X_sn,idx_mass_txfr,idx_vel_txfr,
               idx_h_txfr,HeterogeneousRRates,reactions_parms,fluid_parms,
               solids_parms,nreactions]
            AMREX_GPU_DEVICE (int p_id) noexcept
          {
            MFIXParticleContainer::ParticleType& particle = pstruct[p_id];

            // Local array storing interpolated values
            GpuArray<Real,10+Species::NMAX> interp_loc;
            interp_loc.fill(0.);

            if (grown_bx_is_regular) {

              trilinear_interp(particle.pos(), &interp_loc[0],
                               interp_array, plo, dxi, interp_comp);

            } else { // FAB not all regular

              // Cell containing particle centroid
              const int iloc = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
              const int jloc = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
              const int kloc = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

              // The particle is in a covered cell.
              if (flags_array(iloc,jloc,kloc).isCovered())
              {
                p_realarray[SoArealData::dragx][p_id] = 0.0;
                p_realarray[SoArealData::dragy][p_id] = 0.0;
                p_realarray[SoArealData::dragz][p_id] = 0.0;

                if (solve_enthalpy) {
                  p_realarray[SoArealData::convection][p_id] = 0.0;
                }

                // chemical reaction txfr variables
                if (solve_reactions) {
                  for (int n_g(0); n_g < nspecies_g; n_g++)
                    ptile_data.m_runtime_rdata[idx_mass_txfr+n_g][p_id] = 0.;

                  ptile_data.m_runtime_rdata[idx_vel_txfr+0][p_id] = 0.;
                  ptile_data.m_runtime_rdata[idx_vel_txfr+1][p_id] = 0.;
                  ptile_data.m_runtime_rdata[idx_vel_txfr+2][p_id] = 0.;

                  // Write the result in the enthalpy transfer space
                  ptile_data.m_runtime_rdata[idx_h_txfr][p_id] = 0.;
                }

                return;

              // Cut or regular cell and none of the cells in the stencil is covered
              // (Note we can't assume regular cell has no covered cells in the stencil
              //      because of the diagonal case)
              } else {

                // Upper cell in trilinear stencil
                const int i = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0] + 0.5));
                const int j = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1] + 0.5));
                const int k = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2] + 0.5));

                // All cells in the stencil are regular. Use
                // traditional trilinear interpolation
                if (flags_array(i-1,j-1,k-1).isRegular() &&
                    flags_array(i  ,j-1,k-1).isRegular() &&
                    flags_array(i-1,j  ,k-1).isRegular() &&
                    flags_array(i  ,j  ,k-1).isRegular() &&
                    flags_array(i-1,j-1,k  ).isRegular() &&
                    flags_array(i  ,j-1,k  ).isRegular() &&
                    flags_array(i-1,j  ,k  ).isRegular() &&
                    flags_array(i  ,j  ,k  ).isRegular()) {

                  trilinear_interp(particle.pos(), &interp_loc[0],
                                   interp_array, plo, dxi, interp_comp);

                // At least one of the cells in the stencil is cut or covered
                } else {
#if 0
                  // TODO: This was initially split for variables that may have known
                  // EB values (e.g., no-slip velocity). However, the results changed
                  // more than expected so now EB values are not used.
                  {
                    const int srccomp = 0;
                    const int dstcomp = 0;
                    const int numcomp = 3;

                    shepard_interp_eb(particle.pos(), iloc, jloc, kloc, dx, dxi, plo,
                                      flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                      interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);
                  }
                  {
                    const int srccomp = 3;
                    const int dstcomp = 3;
                    const int numcomp = interp_comp-3;

                    shepard_interp(particle.pos(), iloc, jloc, kloc, dx, dxi, plo,
                                   flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                   interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);
                  }
#else
                  const int srccomp = 0;
                  const int dstcomp = 0;
                  const int numcomp = interp_comp;

                  shepard_interp(particle.pos(), iloc, jloc, kloc, dx, dxi, plo,
                                 flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                 interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);
#endif
                } // Cut cell
              } // Not covered
            } // if box not all regular

            Real pbeta = p_realarray[SoArealData::dragcoeff][p_id];
            Real pvol = p_realarray[SoArealData::volume][p_id];

            // Particle drag calculation.  We multiply the particle velocity
            // by "pmult" so that DEM uses the slip velocity. For PIC we
            // only want the fluid velocity as it uses a pseudo implicit
            // slip velocity for parcels.

            RealVect vel_g(interp_loc[0], interp_loc[1], interp_loc[2]);
            RealVect gp_g(interp_loc[3], interp_loc[4], interp_loc[5]);
            
            Real T_g(0);
            if (solve_enthalpy) {
              T_g = interp_loc[6];
            }

            p_realarray[SoArealData::dragx][p_id] =
              pbeta * (vel_g[0] - pmult*p_realarray[SoArealData::velx][p_id]) -
              (gp_g[0] + gp0_dev[0]) * pvol;

            p_realarray[SoArealData::dragy][p_id] =
              pbeta * (vel_g[1] - pmult*p_realarray[SoArealData::vely][p_id]) -
              (gp_g[1] + gp0_dev[1]) * pvol;

            p_realarray[SoArealData::dragz][p_id] =
              pbeta * (vel_g[2] - pmult*p_realarray[SoArealData::velz][p_id]) -
              (gp_g[2] + gp0_dev[2]) * pvol;

            if(solve_enthalpy) {
              // gamma == (heat transfer coeff) * (particle surface area)
              Real pgamma = p_realarray[SoArealData::convection][p_id];

              p_realarray[SoArealData::convection][p_id] =
                pgamma * (T_g - p_realarray[SoArealData::temperature][p_id]);
            }


            if (solve_reactions) {
              // Extract species mass fractions
              GpuArray<Real,Species::NMAX> X_sn;
              X_sn.fill(0.);

              for (int n_s(0); n_s < nspecies_s; n_s++) {
                X_sn[n_s] = ptile_data.m_runtime_rdata[idx_X_sn + n_s][p_id];
              }

              // Extract interpolated thermodynamic pressure
              const Real ep_g = interp_loc[7];
              const Real ro_g = interp_loc[8];
              Real* X_gk      = &interp_loc[9];
              const Real p_g  = interp_loc[9+nspecies_g];

              const Real ep_s = 1. - ep_g;
              const Real ro_p = p_realarray[SoArealData::density][p_id];

              GpuArray<Real,Reactions::NMAX> R_q_heterogeneous;
              R_q_heterogeneous.fill(0.);

              const RealVect pvel(p_realarray[SoArealData::velx][p_id],
                                  p_realarray[SoArealData::vely][p_id],
                                  p_realarray[SoArealData::velz][p_id]);

              const Real T_p = p_realarray[SoArealData::temperature][p_id];
              const Real DP  = 2. * p_realarray[SoArealData::radius][p_id];

              HeterogeneousRRates.template operator()<run_on>(R_q_heterogeneous.data(),
                                                              reactions_parms,
                                                              solids_parms, X_sn.data(),
                                                              ro_p, ep_s, T_p,
                                                              pvel, fluid_parms,
                                                              X_gk, ro_g, ep_g, T_g,
                                                              vel_g, DP, p_g);

              //***************************************************************
              // Loop over solids species for computing solids txfr rates
              //***************************************************************
              Real G_m_p_heterogeneous(0.);

              for (int n_s(0); n_s < nspecies_s; n_s++) {
                Real G_m_sn_heterogeneous(0.);

                // Loop over reactions to compute each contribution
                for (int q(0); q < nreactions; q++) {

                  const Real stoich_coeff = solids_parms.get_stoich_coeff<run_on>(n_s, q);

                  // Compute solids species n_s transfer rate for reaction q
                  const Real MW_sn = solids_parms.get_MW_sn<run_on>(n_s);
                  Real G_m_sn_q = stoich_coeff * MW_sn * R_q_heterogeneous[q];

                  G_m_sn_heterogeneous += G_m_sn_q;
                }

                ptile_data.m_runtime_rdata[idx_mass_txfr+n_s][p_id] = G_m_sn_heterogeneous;

                G_m_p_heterogeneous += G_m_sn_heterogeneous;
              }

              //***************************************************************
              // Loop over fluid species for computing fluid txfr rates
              //***************************************************************
              Real G_H_g_heterogeneous(0.);

              for (int n_g(0); n_g < nspecies_g; n_g++) {

                // Loop over reactions to compute each contribution
                for (int q(0); q < nreactions; q++) {

                  Real stoich_coeff = fluid_parms.get_stoich_coeff<run_on>(n_g, q);

                  // Compute fluid species n_g transfer rate for reaction q
                  const Real MW_gk = fluid_parms.get_MW_gk<run_on>(n_g);
                  Real G_m_gk_q = stoich_coeff * MW_gk * R_q_heterogeneous[q];

                  // Contribution to the particle
                  const Real h_gk_T_p = fluid_parms.calc_h_gk<run_on>(T_p, n_g);
                  const Real h_gk_T_g = fluid_parms.calc_h_gk<run_on>(T_g, n_g);

                  const Real G_H_pk_q = h_gk_T_g * amrex::min(0., G_m_gk_q);
                  const Real G_H_gk_q = h_gk_T_p * amrex::max(0., G_m_gk_q);

                  G_H_g_heterogeneous += (G_H_pk_q + G_H_gk_q);
                }
              }

              //***************************************************************
              //
              //***************************************************************
              const Real coeff = amrex::max(0., G_m_p_heterogeneous);
              ptile_data.m_runtime_rdata[idx_vel_txfr+0][p_id] = coeff*vel_g[0];
              ptile_data.m_runtime_rdata[idx_vel_txfr+1][p_id] = coeff*vel_g[1];
              ptile_data.m_runtime_rdata[idx_vel_txfr+2][p_id] = coeff*vel_g[2];

              // Write the result in the enthalpy transfer space
              ptile_data.m_runtime_rdata[idx_h_txfr][p_id] = -G_H_g_heterogeneous;
            }

          }); // particle loop
        } // FAB not covered
      } // pti
    } // omp region

    BL_PROFILE_VAR_STOP(particle_deposition);

    delete interp_ptr;

  } // lev

  // Reset velocity Dirichlet bc's to face values
  extrap_dir_bcs = 0;
  mfix_set_velocity_bcs(time, vel_g_in, extrap_dir_bcs);
}
