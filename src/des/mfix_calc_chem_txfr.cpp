#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_K.H>
#include <mfix_filcc.H>

#include <mfix_mf_helpers.H>
#include <mfix_dem_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_algorithm.H>
#include <mfix_leveldata.H>
#include <mfix_des_heterogeneous_rates_K.H>
#include <mfix_deposition_K.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>


void 
mfix::mfix_calc_chem_txfr (const Vector< MultiFab* >& chem_txfr,
                           const Vector< MultiFab* >& ep_g_in,
                           const Vector< MultiFab* >& ro_g_in,
                           const Vector< MultiFab* >& vel_g_in,
                           const Vector< MultiFab* >& T_g_in,
                           const Vector< MultiFab* >& X_gk_in,
                           const Real time)
{
  if (m_deposition_scheme == DepositionScheme::trilinear) {
    mfix_calc_chem_txfr(chem_txfr, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in,
        time, TrilinearDeposition());
  } else if (m_deposition_scheme == DepositionScheme::square_dpvm) {
    mfix_calc_chem_txfr(chem_txfr, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in,
        time, TrilinearDPVMSquareDeposition());
  } else if (m_deposition_scheme == DepositionScheme::true_dpvm) {
    mfix_calc_chem_txfr(chem_txfr, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in,
        time, TrueDPVMDeposition());
  } else if (m_deposition_scheme == DepositionScheme::centroid) {
    mfix_calc_chem_txfr(chem_txfr, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in,
        time, CentroidDeposition());
  } else {
    amrex::Abort("Don't know this deposition_scheme!");
  }
}

template <typename F1>
void 
mfix::mfix_calc_chem_txfr (const Vector< MultiFab* >& chem_txfr,
                           const Vector< MultiFab* >& ep_g_in,
                           const Vector< MultiFab* >& ro_g_in,
                           const Vector< MultiFab* >& vel_g_in,
                           const Vector< MultiFab* >& T_g_in,
                           const Vector< MultiFab* >& X_gk_in,
                           const Real time,
                           F1 WeightFunc)
{
  if (m_reaction_rates_type == ReactionRatesType::RRatesUser) {
    mfix_calc_chem_txfr(chem_txfr, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in,
        time, WeightFunc, ComputeRRateUser());
  } else {
    amrex::Abort("Invalid Reaction Rates Type.");
  }
}

template <typename F1, typename F2>
void 
mfix::mfix_calc_chem_txfr (const Vector< MultiFab* >& /*chem_txfr*/,
                           const Vector< MultiFab* >& ep_g_in,
                           const Vector< MultiFab* >& ro_g_in,
                           const Vector< MultiFab* >& vel_g_in,
                           const Vector< MultiFab* >& T_g_in,
                           const Vector< MultiFab* >& X_gk_in,
                           const Real time,
                           F1 WeightFunc,
                           F2 RRatesFunc)
{
  using PairIndex = MFIXParticleContainer::PairIndex;
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  const Real strttime = ParallelDescriptor::second();

  BL_PROFILE("mfix::mfix_calc_chem_txfr()");

  // Solid species data
  const int nspecies_s = solids.nspecies;

  Gpu::DeviceVector< int > d_species_id_s(nspecies_s);
  Gpu::DeviceVector< Real > d_MW_sn(nspecies_s);
  Gpu::copyAsync(Gpu::hostToDevice, solids.species_id.begin(), solids.species_id.end(), d_species_id_s.begin());
  Gpu::copyAsync(Gpu::hostToDevice, solids.MW_sn0.begin(), solids.MW_sn0.end(), d_MW_sn.begin());
  int* p_species_id_s = d_species_id_s.data();
  Real* p_MW_sn = d_MW_sn.data();

  // Fluid species data
  const int nspecies_g = fluid.nspecies;

  Gpu::DeviceVector< int > d_species_id_g(nspecies_g);
  Gpu::DeviceVector< Real > d_MW_gk(nspecies_g);
  Gpu::copyAsync(Gpu::hostToDevice, fluid.species_id.begin(), fluid.species_id.end(), d_species_id_g.begin());
  Gpu::copyAsync(Gpu::hostToDevice, fluid.MW_gk0.begin(), fluid.MW_gk0.end(), d_MW_gk.begin());
  int* p_species_id_g = d_species_id_g.data();
  Real* p_MW_gk = d_MW_gk.data();

  // Fluid enthalpy data
  Gpu::DeviceVector< Real > d_H_fk0(nspecies_g);
  Gpu::copyAsync(Gpu::hostToDevice, fluid.H_fk0.begin(), fluid.H_fk0.end(), d_H_fk0.begin());
  Real* p_H_fk0 = d_H_fk0.data();

  // Reactions data
  const int nreactions = REACTIONS::nreactions;

  Gpu::DeviceVector< int > d_types(nreactions);
  Gpu::HostVector  < int > h_types(nreactions);
  Gpu::DeviceVector< const int* > d_phases(nreactions);
  Gpu::HostVector  < const int* > h_phases(nreactions);
  Gpu::DeviceVector< int > d_nphases(nreactions);
  Gpu::HostVector  < int > h_nphases(nreactions);

  for (int q(0); q < nreactions; q++) {
    h_types[q] = m_chemical_reactions[q]->m_reaction_type;
    h_phases[q] = m_chemical_reactions[q]->m_phases.data();
    h_nphases[q] = m_chemical_reactions[q]->m_phases.size();
  }

  Gpu::copyAsync(Gpu::hostToDevice, h_types.begin(), h_types.end(), d_types.begin());
  Gpu::copyAsync(Gpu::hostToDevice, h_phases.begin(), h_phases.end(), d_phases.begin());
  Gpu::copyAsync(Gpu::hostToDevice, h_nphases.begin(), h_nphases.end(), d_nphases.begin());

  int* p_types = d_types.data();
  const int** p_phases = d_phases.data();
  int* p_nphases = d_nphases.data();

  Gpu::DeviceVector< int > d_nreactants(nreactions);
  Gpu::HostVector  < int > h_nreactants(nreactions);
  Gpu::DeviceVector< const int* > d_reactants_id(nreactions);
  Gpu::HostVector  < const int* > h_reactants_id(nreactions);
  Gpu::DeviceVector< const Real* > d_reactants_coeffs(nreactions);
  Gpu::HostVector  < const Real* > h_reactants_coeffs(nreactions);
  Gpu::DeviceVector< const int* > d_reactants_phases(nreactions);
  Gpu::HostVector  < const int* > h_reactants_phases(nreactions);

  for (int q(0); q < nreactions; q++) {
    h_nreactants[q] = m_chemical_reactions[q]->m_reactants.size();
    h_reactants_id[q] = m_chemical_reactions[q]->m_reactants_id.data();
    h_reactants_coeffs[q] = m_chemical_reactions[q]->m_reactants_coeffs.data();
    h_reactants_phases[q] = m_chemical_reactions[q]->m_reactants_phases.data();
  }

  Gpu::copyAsync(Gpu::hostToDevice, h_nreactants.begin(), h_nreactants.end(), d_nreactants.begin());
  Gpu::copyAsync(Gpu::hostToDevice, h_reactants_id.begin(), h_reactants_id.end(), d_reactants_id.begin());
  Gpu::copyAsync(Gpu::hostToDevice, h_reactants_coeffs.begin(), h_reactants_coeffs.end(), d_reactants_coeffs.begin());
  Gpu::copyAsync(Gpu::hostToDevice, h_reactants_phases.begin(), h_reactants_phases.end(), d_reactants_phases.begin());

  int* p_nreactants = d_nreactants.data();
  const int** p_reactants_id = d_reactants_id.data();
  const Real** p_reactants_coeffs = d_reactants_coeffs.data();
  const int** p_reactants_phases = d_reactants_phases.data();

  Gpu::DeviceVector< int > d_nproducts(nreactions);
  Gpu::HostVector  < int > h_nproducts(nreactions);
  Gpu::DeviceVector< const int* > d_products_id(nreactions);
  Gpu::HostVector  < const int* > h_products_id(nreactions);
  Gpu::DeviceVector< const Real* > d_products_coeffs(nreactions);
  Gpu::HostVector  < const Real* > h_products_coeffs(nreactions);
  Gpu::DeviceVector< const int* > d_products_phases(nreactions);
  Gpu::HostVector  < const int* > h_products_phases(nreactions);

  for (int q(0); q < nreactions; q++) {
    h_nproducts[q] = m_chemical_reactions[q]->m_products.size();
    h_products_id[q] = m_chemical_reactions[q]->m_products_id.data();
    h_products_coeffs[q] = m_chemical_reactions[q]->m_products_coeffs.data();
    h_products_phases[q] = m_chemical_reactions[q]->m_products_phases.data();
  }

  Gpu::copyAsync(Gpu::hostToDevice, h_nproducts.begin(), h_nproducts.end(), d_nproducts.begin());
  Gpu::copyAsync(Gpu::hostToDevice, h_products_id.begin(), h_products_id.end(), d_products_id.begin());
  Gpu::copyAsync(Gpu::hostToDevice, h_products_coeffs.begin(), h_products_coeffs.end(), d_products_coeffs.begin());
  Gpu::copyAsync(Gpu::hostToDevice, h_products_phases.begin(), h_products_phases.end(), d_products_phases.begin());

  int* p_nproducts = d_nproducts.data();
  const int** p_products_id = d_products_id.data();
  const Real** p_products_coeffs = d_products_coeffs.data();
  const int** p_products_phases = d_products_phases.data();

  // Solid phase integer ID
  const int Solid = CHEMICALPHASE::Solid;
  const int Fluid = CHEMICALPHASE::Fluid;
  const int Heterogeneous = REACTIONTYPE::Heterogeneous;
  const int InvalidIdx = -1; //TODO define this somewhere else

  // Particles SoA starting indexes for mass fractions and rate of formations
  const int idx_X_sn       = (pc->m_runtimeRealData).X_sn;
  const int idx_ro_sn_txfr = (pc->m_runtimeRealData).ro_sn_txfr;
  const int idx_vel_s_txfr = (pc->m_runtimeRealData).vel_s_txfr;
  const int idx_h_s_txfr   = (pc->m_runtimeRealData).h_s_txfr;

  ChemTransfer chem_txfr_idxs(fluid.nspecies, REACTIONS::nreactions);
  const int idx_ro_gk_txfr = chem_txfr_idxs.ro_gk_txfr;
  const int idx_vel_g_txfr = chem_txfr_idxs.vel_g_txfr;
  const int idx_h_g_txfr   = chem_txfr_idxs.h_g_txfr;

  Gpu::DeviceVector< Real > d_H_fn0(nspecies_s);
  Gpu::copyAsync(Gpu::hostToDevice, solids.H_fn0.begin(), solids.H_fn0.end(), d_H_fn0.begin());
  Real* p_H_fn0 = d_H_fn0.data();

  // ************************************************************************
  // Setup data structures for PC deposition
  // ************************************************************************
  for (int lev = 0; lev < nlev; lev++) {
    m_leveldata[lev]->chem_txfr->setVal(0);
  }

  if (nlev > 2)
    amrex::Abort("For right now"
        " MFIXParticleContainer::TrilinearDepositionFluidRRates can only"
        " handle up to 2 levels");

  Vector< MultiFab* > chem_txfr_ptr(nlev, nullptr);

  for (int lev = 0; lev < nlev; lev++)
  {
    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    if (lev == 0 && OnSameGrids)
    {
      // If we are already working with the internal mf defined on the
      // particle_box_array, then we just work with this.
      chem_txfr_ptr[lev] = m_leveldata[lev]->chem_txfr;
    }
    else if (lev == 0 && (! OnSameGrids))
    {
      // If beta_mf is not defined on the particle_box_array, then we need
      // to make a temporary here and copy into beta_mf at the end.
      chem_txfr_ptr[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                        pc->ParticleDistributionMap(lev),
                                        m_leveldata[lev]->chem_txfr->nComp(),
                                        m_leveldata[lev]->chem_txfr->nGrow());
    }
    else
    {
      // If lev > 0 we make a temporary at the coarse resolution
      BoxArray ba_crse(amrex::coarsen(pc->ParticleBoxArray(lev),
            this->m_gdb->refRatio(0)));

      chem_txfr_ptr[lev] = new MultiFab(ba_crse, pc->ParticleDistributionMap(lev),
                                        m_leveldata[lev]->chem_txfr->nComp(), 1);
    }

    // We must have ghost cells for each FAB so that a particle in one grid can
    // spread its effect to an adjacent grid by first putting the value into
    // ghost cells of its own grid.  The mf->sumBoundary call then adds the
    // value from one grid's ghost cell to another grid's valid region.
    if (chem_txfr_ptr[lev]->nGrow() < 1)
      amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

    chem_txfr_ptr[lev]->setVal(0.0, 0, m_leveldata[lev]->chem_txfr->nComp(),
        chem_txfr_ptr[lev]->nGrow());
  }

  const Geometry& gm = Geom(0);
  const FabArray<EBCellFlagFab>* deposition_flags = nullptr;
  const MultiFab* deposition_volfrac = nullptr;

  for (int lev = 0; lev < nlev; lev++)
  {
    // Use level 0 to define the EB factory. If we are not on level 0
    // then create a copy of the coarse factory to use.

    if (lev == 0)
    {
      deposition_flags   = &(particle_ebfactory[lev]->getMultiEBCellFlagFab());
      deposition_volfrac = &(particle_ebfactory[lev]->getVolFrac());
    }
    else
    {
      Vector<int> ngrow = {1,1,1};
      EBFArrayBoxFactory* crse_factory;

      crse_factory = (makeEBFabFactory(gm, chem_txfr_ptr[lev]->boxArray(),
                                       chem_txfr_ptr[lev]->DistributionMap(),
                                       ngrow, EBSupport::volume)).release();

      deposition_flags   = &(crse_factory->getMultiEBCellFlagFab());
      deposition_volfrac = &(crse_factory->getVolFrac());

      delete crse_factory;
    }
  }

  // **************************************************************************
  // Compute particles densities transfer rates
  // **************************************************************************
  // Extrapolate velocity Dirichlet bc's to ghost cells
  int extrap_dir_bcs = 1;

  mfix_set_velocity_bcs(time, get_vel_g(), extrap_dir_bcs);

  // We copy the value inside the domain to the outside to avoid
  // unphysical volume fractions.
  const int dir_bc_in = 2;
  mfix_set_epg_bcs(ep_g_in, dir_bc_in);

  // Set boundary conditions just in case
  mfix_set_density_bcs(time, ro_g_in);
  mfix_set_temperature_bcs(time, T_g_in);
  mfix_set_species_bcs(time, X_gk_in);

  for (int lev = 0; lev < nlev; lev++)
  {
    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    MultiFab* interp_ptr;

    EB_set_covered(*vel_g_in[lev], 0, 3, 1, covered_val);
    EB_set_covered(*ep_g_in[lev], 0, 1, 1, covered_val);
    EB_set_covered(*ro_g_in[lev], 0, 1, 1, covered_val);
    EB_set_covered(*T_g_in[lev], 0, 1, 1, covered_val);
    EB_set_covered(*X_gk_in[lev], 0, fluid.nspecies, 1, covered_val);

    const int interp_ng = 1;  // Only one layer needed for interpolation
    const int interp_comp = 6+nspecies_g; // 3 vel_g + ep_g + ro_g + T_g + X_gk

    if (OnSameGrids)
    {
      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(grids[lev], dmap[lev], interp_comp, interp_ng,
          MFInfo(), *ebfactory[lev]);

      // Copy velocity
      MultiFab::Copy(*interp_ptr, *vel_g_in[lev], 0, 0, 3, interp_ng);

      // Copy volume fraction
      MultiFab::Copy(*interp_ptr, *ep_g_in[lev], 0, 3, 1, interp_ng);

      // Copy density
      MultiFab::Copy(*interp_ptr, *ro_g_in[lev], 0, 4, 1, interp_ng);

      // Copy temperature
      MultiFab::Copy(*interp_ptr, *T_g_in[lev], 0, 5, 1, interp_ng);

      // Copy X_gk
      MultiFab::Copy(*interp_ptr, *X_gk_in[lev], 0, 6, nspecies_g, interp_ng);
    }
    else
    {
      const BoxArray&            pba = pc->ParticleBoxArray(lev);
      const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

      EBFArrayBoxFactory ebfactory_loc(*eb_levels[lev], geom[lev], pba, pdm,
                                       {nghost_eb_basic(), nghost_eb_volume(),
                                        nghost_eb_full()}, EBSupport::full);

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(pba, pdm, interp_comp, interp_ng, MFInfo(),
          ebfactory_loc);

      // Copy velocity
      interp_ptr->ParallelCopy(*vel_g_in[lev], 0, 0, 3, interp_ng, interp_ng);

      // Copy volume fraction
      interp_ptr->ParallelCopy(*ep_g_in[lev], 0, 3, 1, interp_ng, interp_ng);

      // Copy density
      interp_ptr->ParallelCopy(*ro_g_in[lev], 0, 4, 1, interp_ng, interp_ng);

      // Copy temperature
      interp_ptr->ParallelCopy(*T_g_in[lev], 0, 5, 1, interp_ng, interp_ng);

      // Copy X_gk
      interp_ptr->ParallelCopy(*X_gk_in[lev], 0, 6, fluid.nspecies, interp_ng, interp_ng);
    }

    // FillBoundary on interpolation MultiFab
    interp_ptr->FillBoundary(geom[lev].periodicity());

    // Do the interpolate
    {
      const auto dxi_array = geom[lev].InvCellSizeArray();
      const auto dx_array  = geom[lev].CellSizeArray();
      const auto plo_array = geom[lev].ProbLoArray();

      const amrex::RealVect  dx( dx_array[0],  dx_array[1],  dx_array[2]);
      const amrex::RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
      const amrex::RealVect plo(plo_array[0], plo_array[1], plo_array[2]);

      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(interp_ptr->Factory());

      const auto cellcent = &(factory.getCentroid());
      const auto bndrycent = &(factory.getBndryCent());
      const auto areafrac = factory.getAreaFrac();

      FArrayBox local_fab_chem_txfr;

      auto& fluid_parms = *fluid.parameters;
      auto& solids_parms = *solids.parameters;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
      {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto& ptile = pc->GetParticles(lev)[index];

        //Access to added variables
        auto ptile_data = ptile.getParticleTileData();

        auto& particles = pti.GetArrayOfStructs();
        MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        auto& soa = pti.GetStructOfArrays();
        auto p_realarray = soa.realarray();

        const int np = particles.size();

        // *******************************************************************
        // Reaction rates for each particle
        // *******************************************************************
        Gpu::DeviceVector<Real> G_sk_pg(np*nspecies_g);
        Real* G_sk_pg_ptr = G_sk_pg.dataPtr();

        Gpu::DeviceVector<Real> G_h_pg(np);
        Real* G_h_pg_ptr = G_h_pg.dataPtr();

        Box bx = pti.tilebox();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox& interp_fab = static_cast<EBFArrayBox const&>((*interp_ptr)[pti]);

        const EBCellFlagFab& flags = interp_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered)
        {
          const auto& interp_array = interp_ptr->const_array(pti);

          const auto& flags_array = flags.const_array();

          if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
          {
            amrex::ParallelFor(np,
              [np,pstruct,p_realarray,interp_array,RRatesFunc,G_sk_pg_ptr,
               G_h_pg_ptr,plo,dxi,ptile_data,nspecies_g,nspecies_s,nreactions,
               interp_comp,Solid,p_MW_sn,idx_X_sn,idx_ro_sn_txfr,idx_vel_s_txfr,
               idx_h_s_txfr,p_MW_gk,p_species_id_s,p_species_id_g,p_reactants_id,
               p_reactants_coeffs,p_reactants_phases,p_products_id,p_products_coeffs,
               p_products_phases,p_nreactants,p_nproducts,InvalidIdx,p_phases,
               p_nphases,p_types,Heterogeneous,p_H_fk0,p_H_fn0,fluid_parms,solids_parms]
              AMREX_GPU_DEVICE (int p_id) noexcept
            {
              auto& particle = pstruct[p_id];

              GpuArray<Real,SPECIES::NMAX> X_sn;

              for (int n_s(0); n_s < nspecies_s; n_s++) {
                const int idx = idx_X_sn + n_s;
                X_sn[n_s] = ptile_data.m_runtime_rdata[idx][p_id];
              }

              const Real T_p = p_realarray[SoArealData::temperature][p_id];

              GpuArray<Real,6+SPECIES::NMAX> interp_loc; // vel_g, ep_g, ro_g, T_g, X_gk
              interp_loc.fill(0.);

              int ip(-1); int jp(-1); int kp(-1);
              GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

              trilinear_interp(particle.pos(), ip, jp, kp, weights,
                               interp_loc.data(), interp_array, plo, dxi,
                               interp_comp);

              RealVect vel_g(interp_loc[0], interp_loc[1], interp_loc[2]);
              const Real ep_g = interp_loc[3];
              const Real ro_g = interp_loc[4];
              const Real T_g  = interp_loc[5];

              //GpuArray<Real,SPECIES::NMAX> X_gk;
              Real* X_gk = &interp_loc[6];

              //for (int n_g(0); n_g < nspecies_g; n_g++) {
              //  X_gk[n_g] = interp_loc[6+n_g];
              //}

              const Real ep_s = 1. - ep_g;
              const Real ro_s = p_realarray[SoArealData::density][p_id];

              GpuArray<Real,REACTIONS::NMAX> R_q;
              R_q.fill(0.);

              RRatesFunc(R_q.data(), nreactions, p_nreactants, p_nproducts,
                  p_reactants_id, p_reactants_coeffs, p_reactants_phases,
                  p_products_id,  p_products_coeffs,  p_products_phases,
                  p_species_id_s, X_sn.data(), p_MW_sn, nspecies_s, ro_s, ep_s,
                  p_species_id_g, X_gk,        p_MW_gk, nspecies_g, ro_g, ep_g);

              // Total transfer rates
              Real G_s_pg(0);
              Real G_s_gp(0);

              Real G_h_pg_loc(0);
              Real G_h_gp(0);


              //***************************************************************
              // Initialize to zero
              //***************************************************************
              for (int n_s(0); n_s < nspecies_s; n_s++) {
                // Initially set species n_s density transfer rate to zero
                ptile_data.m_runtime_rdata[idx_ro_sn_txfr + n_s][p_id] = 0;
              }

              ptile_data.m_runtime_rdata[idx_vel_s_txfr + 0][p_id] = 0.;
              ptile_data.m_runtime_rdata[idx_vel_s_txfr + 1][p_id] = 0.;
              ptile_data.m_runtime_rdata[idx_vel_s_txfr + 2][p_id] = 0.;

              // Initially set particle energy transfer rate to zero
              ptile_data.m_runtime_rdata[idx_h_s_txfr][p_id] = 0.;

              //***************************************************************
              // Loop over particle's species for computing particle txfr rates
              //***************************************************************
              for (int n_s(0); n_s < nspecies_s; n_s++) {
                Real G_sn_gp(0);

                // Get the ID of the current species n_s
                const int current_species_id = p_species_id_s[n_s];

                // Loop over reactions to compute each contribution
                for (int q(0); q < nreactions; q++) {
                  // Do something only if reaction is heterogeneous and contains
                  // a solid compound
                  if (p_types[q] == Heterogeneous &&
                      MFIXfind(p_phases[q], p_nphases[q], Solid) != InvalidIdx) {

                    Real stoc_coeff(0);

                    // Add reactant contribution (if any)
                    {
                      const int pos = MFIXfind(p_reactants_id[q], p_nreactants[q], current_species_id);

                      if (pos != InvalidIdx) {
                        if (p_reactants_phases[q][pos] == Solid)
                          stoc_coeff += p_reactants_coeffs[q][pos];
                      }
                    }

                    // Add products contribution (if any)
                    {
                      const int pos = MFIXfind(p_products_id[q], p_nproducts[q], current_species_id);

                      if (pos != InvalidIdx) {
                        if (p_products_phases[q][pos] == Solid)
                          stoc_coeff += p_products_coeffs[q][pos];
                      }
                    }

                    
                    // Compute particle's species n_s transfer rate for reaction q
                    Real G_sn_gp_q = stoc_coeff * p_MW_sn[n_s] * R_q[q];

                    G_sn_gp += G_sn_gp_q;
                  }
                }

                G_s_gp += G_sn_gp;

                const Real h_sn_T_p = p_H_fn0[n_s] + solids_parms.calc_h_sn<RunOn::Gpu>(T_p,n_s);
                G_h_gp += h_sn_T_p * G_sn_gp;

                // Update global variable
                ptile_data.m_runtime_rdata[idx_ro_sn_txfr + n_s][p_id] = G_sn_gp;
              }


              //***************************************************************
              // Initialize to zero
              //***************************************************************
              for (int n_g(0); n_g < nspecies_g; n_g++) {
                // Initially set fluid species n_g density transfer rate to zero
                G_sk_pg_ptr[n_g*np + p_id] = 0;
              }

              // Initially set fluid energy transfer rate to zero
              G_h_pg_ptr[p_id] = 0;

              //***************************************************************
              // Loop over fluid species for computing fluid txfr rates
              //***************************************************************
              for (int n_g(0); n_g < nspecies_g; n_g++) {
                Real G_sk_pg_loc(0.);

                // Get the ID of the current species n_g
                const int current_species_id = p_species_id_g[n_g];

                // Loop over reactions to compute each contribution
                for (int q(0); q < nreactions; q++) {
                  // Do something only if reaction is heterogeneous and contains
                  // a solid compound
                  if (p_types[q] == Heterogeneous &&
                      MFIXfind(p_phases[q], p_nphases[q], Fluid) != InvalidIdx) {

                    Real stoc_coeff(0);

                    // Add reactant contribution (if any)
                    {
                      const int pos = MFIXfind(p_reactants_id[q], p_nreactants[q], current_species_id);

                      if (pos != InvalidIdx) {
                        if (p_reactants_phases[q][pos] == Fluid)
                          stoc_coeff += p_reactants_coeffs[q][pos];
                      }
                    }

                    // Add products contribution (if any)
                    {
                      const int pos = MFIXfind(p_products_id[q], p_nproducts[q], current_species_id);

                      if (pos != InvalidIdx) {
                        if (p_products_phases[q][pos] == Fluid)
                          stoc_coeff += p_products_coeffs[q][pos];
                      }
                    }

                    // Compute particle's species n_s transfer rate for reaction q
                    Real G_sk_pg_q = stoc_coeff * p_MW_gk[n_g] * R_q[q];

                    G_sk_pg_loc += G_sk_pg_q;

                    // Contribution to the particle
                    const Real h_gk_T_p = p_H_fk0[n_g] + fluid_parms.calc_h_gk<RunOn::Gpu>(T_p,n_g);
                    const Real h_gk_T_g = p_H_fk0[n_g] + fluid_parms.calc_h_gk<RunOn::Gpu>(T_g,n_g);

                    G_h_gp += h_gk_T_p * G_sk_pg_q;
                    G_h_gp += amrex::min(0., G_sk_pg_q) * (h_gk_T_g - h_gk_T_p);

                    // Contribution to fluid energy transfer
                    G_h_pg_loc += amrex::max(0., G_sk_pg_q) * (h_gk_T_g - h_gk_T_p);
                  }
                }

                G_s_pg += G_sk_pg_loc;

                // Update global variable
                G_sk_pg_ptr[n_g*np + p_id] = G_sk_pg_loc;
              }

              // Check that (total_ro_g_txfr + total_ro_p_txfr) = 0
              //BL_ASSERT(std::abs(G_s_gp + G_s_pg) < 1.e-15);

              //***************************************************************
              // Update particle linear momentum and energy transfer
              //***************************************************************
              const Real coeff = amrex::max(0., G_s_gp);

              // NOTE: total_ro_p_txfr is computed from interpolated quantities
              // and vel_g is also interpolated to particle position
              ptile_data.m_runtime_rdata[idx_vel_s_txfr+0][p_id] = coeff*vel_g[0];
              ptile_data.m_runtime_rdata[idx_vel_s_txfr+1][p_id] = coeff*vel_g[1];
              ptile_data.m_runtime_rdata[idx_vel_s_txfr+2][p_id] = coeff*vel_g[2];

              // Write the result in the enthalpy transfer space
              ptile_data.m_runtime_rdata[idx_h_s_txfr][p_id] = G_h_gp;
              G_h_pg_ptr[p_id] = G_h_pg_loc;
            });
          } else { // FAB not all regular

            // Cell centroids
            const auto& ccent_fab = cellcent->array(pti);
            // Centroid of EB
            const auto& bcent_fab = bndrycent->array(pti);
            // Area fractions
            const auto& apx_fab = areafrac[0]->array(pti);
            const auto& apy_fab = areafrac[1]->array(pti);
            const auto& apz_fab = areafrac[2]->array(pti);

            amrex::ParallelFor(np,
              [np,pstruct,p_realarray,interp_array,RRatesFunc,plo,dx,G_h_pg_ptr,
               G_sk_pg_ptr,dxi,flags_array,ccent_fab,bcent_fab,apx_fab,apy_fab,apz_fab,
               nspecies_s,ptile_data,interp_comp,nspecies_g,nreactions,p_species_id_s,
               p_species_id_g,p_types,p_phases,p_nphases,p_products_id,p_products_coeffs,
               p_products_phases,p_reactants_id,p_reactants_coeffs,p_reactants_phases,
               InvalidIdx,p_MW_sn,idx_X_sn,idx_ro_sn_txfr,idx_vel_s_txfr,idx_h_s_txfr,p_MW_gk,Solid,
               p_nreactants,p_nproducts,Heterogeneous,p_H_fk0,p_H_fn0,fluid_parms,solids_parms]
              AMREX_GPU_DEVICE (int p_id) noexcept
            {
              auto& particle = pstruct[p_id];

              // Cell containing particle centroid
              int ip = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
              int jp = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
              int kp = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

              // No density transfer rate for particles in covered cells.
              if (flags_array(ip,jp,kp).isCovered()) {
                // Particle mass calculation
                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  // Update specie mass formation rate
                  ptile_data.m_runtime_rdata[idx_ro_sn_txfr + n_s][p_id] = 0;
                }

                ptile_data.m_runtime_rdata[idx_vel_s_txfr + 0][p_id] = 0;
                ptile_data.m_runtime_rdata[idx_vel_s_txfr + 1][p_id] = 0;
                ptile_data.m_runtime_rdata[idx_vel_s_txfr + 2][p_id] = 0;

                ptile_data.m_runtime_rdata[idx_h_s_txfr][p_id] = 0;
              } else {
                // Cut or regular cell and none of the cells in the stencil is
                // covered (Note we can't assume regular cell has no covered
                // cells in the stencil because of the diagonal case)

                GpuArray<Real,SPECIES::NMAX> X_sn;

                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  const int idx = idx_X_sn + n_s;
                  X_sn[n_s] = ptile_data.m_runtime_rdata[idx][p_id];
                }

                const Real T_p = p_realarray[SoArealData::temperature][p_id];

                GpuArray<Real,6+SPECIES::NMAX> interp_loc; // vel_g, ep_g, ro_g, T_g, X_gk
                interp_loc.fill(0.);

                GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

                // Upper cell in trilinear stencil
                int i = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0] + 0.5));
                int j = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1] + 0.5));
                int k = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2] + 0.5));

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

                  trilinear_interp(particle.pos(), ip, jp, kp, weights,
                                   interp_loc.data(), interp_array, plo, dxi,
                                   interp_comp);
                } else { // At least one of the cells in the stencil is cut or covered

                  // All the quantities to interpolate are scalar quantities
                  // No velocity to be interpolated
                  const int scomp = 3;

                  fe_interp(particle.pos(), ip, jp, kp, weights, dx, dxi, plo,
                            flags_array, ccent_fab, bcent_fab, apx_fab,
                            apy_fab, apz_fab, interp_array, interp_loc.data(),
                            interp_comp, scomp);
                } // Cut cell

                RealVect vel_g(interp_loc[0], interp_loc[1], interp_loc[2]);
                const Real ep_g = interp_loc[3];
                const Real ro_g = interp_loc[4];
                const Real T_g  = interp_loc[5];

                //GpuArray<Real,SPECIES::NMAX> X_gk;
                Real* X_gk = &interp_loc[6];

                //for (int n_g(0); n_g < nspecies_g; ++n_g) {
                //  X_gk[n_g] = interp_loc[6+n_g];
                //}

                const Real ep_s = 1. - ep_g;
                const Real ro_s = p_realarray[SoArealData::density][p_id];

                GpuArray<Real,REACTIONS::NMAX> R_q;
                R_q.fill(0.);

                RRatesFunc(R_q.data(), nreactions, p_nreactants, p_nproducts,
                    p_reactants_id, p_reactants_coeffs, p_reactants_phases,
                    p_products_id,  p_products_coeffs,  p_products_phases,
                    p_species_id_s, X_sn.data(), p_MW_sn, nspecies_s, ro_s, ep_s,
                    p_species_id_g, X_gk,        p_MW_gk, nspecies_g, ro_g, ep_g);

                // Total transfer rates
                Real G_s_pg(0);
                Real G_s_gp(0);

                Real G_h_pg_loc(0);
                Real G_h_gp(0);


                //***************************************************************
                // Initialize to zero
                //***************************************************************
                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  // Initially set species n_s density transfer rate to zero
                  ptile_data.m_runtime_rdata[idx_ro_sn_txfr + n_s][p_id] = 0;
                }

                ptile_data.m_runtime_rdata[idx_vel_s_txfr + 0][p_id] = 0.;
                ptile_data.m_runtime_rdata[idx_vel_s_txfr + 1][p_id] = 0.;
                ptile_data.m_runtime_rdata[idx_vel_s_txfr + 2][p_id] = 0.;

                // Initially set particle energy transfer rate to zero
                ptile_data.m_runtime_rdata[idx_h_s_txfr][p_id] = 0.;

                //***************************************************************
                // Loop over particle's species for computing particle txfr rates
                //***************************************************************
                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  Real G_sn_gp(0);

                  // Get the ID of the current species n_s
                  const int current_species_id = p_species_id_s[n_s];

                  // Loop over reactions to compute each contribution
                  for (int q(0); q < nreactions; q++) {
                    // Do something only if reaction is heterogeneous and contains
                    // a solid compound
                    if (p_types[q] == Heterogeneous &&
                        MFIXfind(p_phases[q], p_nphases[q], Solid) != InvalidIdx) {

                      Real stoc_coeff(0);

                      // Add reactant contribution (if any)
                      {
                        const int pos = MFIXfind(p_reactants_id[q], p_nreactants[q], current_species_id);

                        if (pos != InvalidIdx) {
                          if (p_reactants_phases[q][pos] == Solid)
                            stoc_coeff += p_reactants_coeffs[q][pos];
                        }
                      }

                      // Add products contribution (if any)
                      {
                        const int pos = MFIXfind(p_products_id[q], p_nproducts[q], current_species_id);

                        if (pos != InvalidIdx) {
                          if (p_products_phases[q][pos] == Solid)
                            stoc_coeff += p_products_coeffs[q][pos];
                        }
                      }

                      
                      // Compute particle's species n_s transfer rate for reaction q
                      Real G_sn_gp_q = stoc_coeff * p_MW_sn[n_s] * R_q[q];

                      G_sn_gp += G_sn_gp_q;
                    }
                  }

                  G_s_gp += G_sn_gp;

                  const Real h_sn_T_p = p_H_fn0[n_s] + solids_parms.calc_h_sn<RunOn::Gpu>(T_p,n_s);
                  G_h_gp += h_sn_T_p * G_sn_gp;

                  // Update global variable
                  ptile_data.m_runtime_rdata[idx_ro_sn_txfr + n_s][p_id] = G_sn_gp;
                }


                //***************************************************************
                // Initialize to zero
                //***************************************************************
                for (int n_g(0); n_g < nspecies_g; n_g++) {
                  // Initially set fluid species n_g density transfer rate to zero
                  G_sk_pg_ptr[n_g*np + p_id] = 0;
                }

                // Initially set fluid energy transfer rate to zero
                G_h_pg_ptr[p_id] = 0;

                //***************************************************************
                // Loop over fluid species for computing fluid txfr rates
                //***************************************************************
                for (int n_g(0); n_g < nspecies_g; n_g++) {
                  Real G_sk_pg_loc(0.);

                  // Get the ID of the current species n_g
                  const int current_species_id = p_species_id_g[n_g];

                  // Loop over reactions to compute each contribution
                  for (int q(0); q < nreactions; q++) {
                    // Do something only if reaction is heterogeneous and contains
                    // a solid compound
                    if (p_types[q] == Heterogeneous &&
                        MFIXfind(p_phases[q], p_nphases[q], Fluid) != InvalidIdx) {

                      Real stoc_coeff(0);

                      // Add reactant contribution (if any)
                      {
                        const int pos = MFIXfind(p_reactants_id[q], p_nreactants[q], current_species_id);

                        if (pos != InvalidIdx) {
                          if (p_reactants_phases[q][pos] == Fluid)
                            stoc_coeff += p_reactants_coeffs[q][pos];
                        }
                      }

                      // Add products contribution (if any)
                      {
                        const int pos = MFIXfind(p_products_id[q], p_nproducts[q], current_species_id);

                        if (pos != InvalidIdx) {
                          if (p_products_phases[q][pos] == Fluid)
                            stoc_coeff += p_products_coeffs[q][pos];
                        }
                      }

                      // Compute particle's species n_s transfer rate for reaction q
                      Real G_sk_pg_q = stoc_coeff * p_MW_gk[n_g] * R_q[q];

                      G_sk_pg_loc += G_sk_pg_q;

                      // Contribution to the particle
                      const Real h_gk_T_p = p_H_fk0[n_g] + fluid_parms.calc_h_gk<RunOn::Gpu>(T_p,n_g);
                      const Real h_gk_T_g = p_H_fk0[n_g] + fluid_parms.calc_h_gk<RunOn::Gpu>(T_g,n_g);

                      G_h_gp += h_gk_T_p * G_sk_pg_q;
                      G_h_gp += amrex::min(0., G_sk_pg_q) * (h_gk_T_g - h_gk_T_p);

                      // Contribution to fluid energy transfer
                      G_h_pg_loc += amrex::max(0., G_sk_pg_q) * (h_gk_T_g - h_gk_T_p);
                    }
                  }

                  G_s_pg += G_sk_pg_loc;

                  // Update global variable
                  G_sk_pg_ptr[n_g*np + p_id] = G_sk_pg_loc;
                }

                // Check that (total_ro_g_txfr + total_ro_p_txfr) = 0
                //Print() << "G_s_gp = " << G_s_gp << "\n";
                //Print() << "G_s_pg = " << G_s_pg << "\n";
                //BL_ASSERT(std::abs(G_s_gp + G_s_pg) < 1.e-15);

                //***************************************************************
                // Update particle linear momentum and energy transfer
                //***************************************************************
                const Real coeff = amrex::max(0., G_s_gp);

                // NOTE: total_ro_p_txfr is computed from interpolated quantities
                // and vel_g is also interpolated to particle position
                ptile_data.m_runtime_rdata[idx_vel_s_txfr+0][p_id] = coeff*vel_g[0];
                ptile_data.m_runtime_rdata[idx_vel_s_txfr+1][p_id] = coeff*vel_g[1];
                ptile_data.m_runtime_rdata[idx_vel_s_txfr+2][p_id] = coeff*vel_g[2];

                // Write the result in the enthalpy transfer space
                ptile_data.m_runtime_rdata[idx_h_s_txfr][p_id] = G_h_gp;
                G_h_pg_ptr[p_id] = G_h_pg_loc;
              } // Not covered
            }); // p_id

          } // type of FAB
        } // if entire FAB not covered

        // Synchronization is needed because the reaction rates need to be fully
        // updated before we start the deposition algorithm
        Gpu::synchronize();

        // ******************************************************************
        // Deposit the interphase chemical transfer quantities to the grid
        // Density transfer:
        // Momentum transfer:
        // Energy transfer
        // ******************************************************************
        FArrayBox& fab_chem_txfr = (*chem_txfr_ptr[lev])[pti];

        const Box& box = pti.tilebox(); // I need a box without ghosts

        if ((*deposition_flags)[pti].getType(box) != FabType::covered) {
          auto arr_chem_txfr   = fab_chem_txfr.array();
          const auto& flagsarr = (*deposition_flags)[pti].array();
          const auto& vfrac    = (*deposition_volfrac)[pti].array();

          const amrex::Real deposition_scale_factor = m_deposition_scale_factor;

#ifdef _OPENMP
          const int ncomp = chem_txfr_ptr[lev]->nComp();
          Box tile_box = box;

          if (Gpu::notInLaunchRegion()) {
            tile_box.grow(chem_txfr_ptr[lev]->nGrow());

            local_fab_chem_txfr.resize(tile_box, ncomp);
            local_fab_chem_txfr.setVal<RunOn::Host>(0.0);
            arr_chem_txfr = local_fab_chem_txfr.array();
          }
#endif

          const auto reg_cell_vol = dx[0]*dx[1]*dx[2];

          const long nrp = pti.numParticles();

          amrex::ParallelFor(nrp,
            [nrp,pstruct,p_realarray,ptile_data,plo_array,dx_array,dxi_array,
             vfrac,flagsarr,deposition_scale_factor,reg_cell_vol,WeightFunc,
             arr_chem_txfr,idx_ro_sn_txfr,idx_ro_gk_txfr,idx_vel_g_txfr,
             idx_h_g_txfr,G_sk_pg_ptr,G_h_pg_ptr,nspecies_g]
            AMREX_GPU_DEVICE (int p_id) noexcept
          {
            const auto& p = pstruct[p_id];

            GpuArray<Real,SPECIES::NMAX> G_sk_pg_loc;
            G_sk_pg_loc.fill(0.);

            Real G_s_pg(0.);

            for (int n_g(0); n_g < nspecies_g; ++n_g) {
              G_sk_pg_loc[n_g] = G_sk_pg_ptr[n_g*nrp + p_id];
              G_s_pg += G_sk_pg_loc[n_g];
            }

            const Real G_h_pg_loc = G_h_pg_ptr[p_id];

            const Real coeff = amrex::max(0., G_s_pg);

            int i(0);
            int j(0);
            int k(0);

            GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

            WeightFunc(plo_array, dx_array, dxi_array, flagsarr, p.pos(),
              p_realarray[SoArealData::radius][p_id], i, j, k, weights,
              deposition_scale_factor);

            const RealVect vel_p(p_realarray[SoArealData::velx][p_id],
                                 p_realarray[SoArealData::vely][p_id],
                                 p_realarray[SoArealData::velz][p_id]);

            // Deposition
            for (int ii = -1; ii <= 0; ++ii) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int kk = -1; kk <= 0; ++kk) {
                  if (! flagsarr(i+ii,j+jj,k+kk).isCovered()) {

                    amrex::Real weight_vol = weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk);

                    for (int n_g(0); n_g < nspecies_g; n_g++) {
                      // Deposition of Rrates
                      Gpu::Atomic::Add(&arr_chem_txfr(i+ii,j+jj,k+kk,idx_ro_gk_txfr+n_g), weight_vol*G_sk_pg_loc[n_g]);
                    }

                    if (coeff > 0.) {
                      // Deposition of linear momentum
                      Gpu::Atomic::Add(&arr_chem_txfr(i+ii,j+jj,k+kk,idx_vel_g_txfr+0), weight_vol*coeff*vel_p[0]);
                      Gpu::Atomic::Add(&arr_chem_txfr(i+ii,j+jj,k+kk,idx_vel_g_txfr+1), weight_vol*coeff*vel_p[1]);
                      Gpu::Atomic::Add(&arr_chem_txfr(i+ii,j+jj,k+kk,idx_vel_g_txfr+2), weight_vol*coeff*vel_p[2]);
                    }

                    Gpu::Atomic::Add(&arr_chem_txfr(i+ii,j+jj,k+kk,idx_h_g_txfr), weight_vol*G_h_pg_loc);
                  }
                }
              }
            }
          });

#ifdef _OPENMP
          if (Gpu::notInLaunchRegion()) {
            fab_chem_txfr.atomicAdd<RunOn::Host>(local_fab_chem_txfr,
                tile_box, tile_box, 0, 0, ncomp);
          }
#endif
        }

      } // pti
    } // GPU region

    delete interp_ptr;
  } // lev

  // Reset the volume fractions back to the correct values at
  // inflow faces.
  const int dir_bc_out = 1;
  mfix_set_epg_bcs(ep_g_in, dir_bc_out);
  // End compute particles density transfer rates

  // **************************************************************************
  // Deposit particles density transfer rates into fluid and compute fluid
  // density transfer rates
  // **************************************************************************
  {
    // The deposition occurred on level 0, thus the next few operations
    // only need to be carried out on level 0.
    int lev(0);

    // Move any volume deposited outside the domain back into the domain
    // when BC is either a pressure inlet or mass inflow.
    mfix_deposition_bcs(lev, *chem_txfr_ptr[lev]);

    // Sum grid boundaries to capture any material that was deposited into
    // your grid from an adjacent grid.
    chem_txfr_ptr[lev]->SumBoundary(gm.periodicity());
    chem_txfr_ptr[lev]->FillBoundary(gm.periodicity());
  }

  int  src_nghost = 1;
  int dest_nghost = 0;
  int ng_to_copy = amrex::min(src_nghost, dest_nghost);

  for (int lev = 1; lev < nlev; lev++) {
    chem_txfr_ptr[0]->ParallelCopy(*chem_txfr_ptr[lev], 0, 0, chem_txfr_ptr[0]->nComp(),
        ng_to_copy, ng_to_copy, gm.periodicity(), FabArrayBase::ADD);
  }

  if (nlev > 1)
  {
    // IntVect ref_ratio(this->m_gdb->refRatio(0));

    // Now interpolate from the coarse grid to define the fine grid ep-g
    Interpolater* mapper = &cell_cons_interp;
    int lo_bc[3] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
    int hi_bc[3] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
    Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

    BndryFuncArray bfunc(mfix_aux::filcc);

    for (int lev = 1; lev < nlev; lev++)
    {
      PhysBCFunct<BndryFuncArray> cphysbc(Geom(lev-1), bcs, bfunc);
      PhysBCFunct<BndryFuncArray> fphysbc(Geom(lev  ), bcs, bfunc);

      m_leveldata[lev]->chem_txfr->setVal(0);

      amrex::InterpFromCoarseLevel(*m_leveldata[lev]->chem_txfr, time,
                                   *chem_txfr_ptr[lev-1],
                                   0, 0, 1, Geom(lev-1), Geom(lev),
                                   cphysbc, 0, fphysbc, 0,
                                   ref_ratio[0], mapper,
                                   bcs, 0);
    }
  }

  // If mf_to_be_filled is not defined on the particle_box_array, then we need
  // to copy here from txfr_ptr into mf_to_be_filled. I believe that we don't
  // need any information in ghost cells so we don't copy those.

  if (chem_txfr_ptr[0] != m_leveldata[0]->chem_txfr) {
    m_leveldata[0]->chem_txfr->ParallelCopy(*chem_txfr_ptr[0], 0, 0,
        m_leveldata[0]->chem_txfr->nComp());
  }

  for (int lev = 0; lev < nlev; lev++) {
    if (chem_txfr_ptr[lev] != m_leveldata[lev]->chem_txfr)
      delete chem_txfr_ptr[lev];
  }

  if (m_verbose > 1) {
    Real stoptime = ParallelDescriptor::second() - strttime;

    ParallelDescriptor::ReduceRealMax(stoptime,
                                      ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "MFIXParticleContainer::TrilinearChemDepositionFluid"
      " time: " << stoptime << '\n';
  }

  // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
  for (int lev = 0; lev < nlev; lev++)
    m_leveldata[lev]->chem_txfr->FillBoundary(geom[lev].periodicity());
  // End compute fluid density transfer rate from particles transfer rates
  // deposition
}
