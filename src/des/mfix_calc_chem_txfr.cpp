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
#include <mfix_reactions_rates_K.H>
#include <mfix_deposition_K.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>


void 
mfix::mfix_calc_chem_txfr (const Vector< MultiFab* >& chem_txfr,
                           const Vector< MultiFab* >& ep_g_in,
                           const Vector< MultiFab* >& ro_g_in,
                           const Vector< MultiFab* >& vel_g_in,
                           const Vector< MultiFab* >& p_g_in,
                           const Vector< MultiFab* >& T_g_in,
                           const Vector< MultiFab* >& X_gk_in,
                           const Real time)
{
  if (m_deposition_scheme == DepositionScheme::trilinear) {
    mfix_calc_chem_txfr(chem_txfr, ep_g_in, ro_g_in, vel_g_in, p_g_in, T_g_in, X_gk_in,
        time, TrilinearDeposition());
  } else if (m_deposition_scheme == DepositionScheme::square_dpvm) {
    mfix_calc_chem_txfr(chem_txfr, ep_g_in, ro_g_in, vel_g_in, p_g_in, T_g_in, X_gk_in,
        time, TrilinearDPVMSquareDeposition());
  } else if (m_deposition_scheme == DepositionScheme::true_dpvm) {
    mfix_calc_chem_txfr(chem_txfr, ep_g_in, ro_g_in, vel_g_in, p_g_in, T_g_in, X_gk_in,
        time, TrueDPVMDeposition());
  } else if (m_deposition_scheme == DepositionScheme::centroid) {
    mfix_calc_chem_txfr(chem_txfr, ep_g_in, ro_g_in, vel_g_in, p_g_in, T_g_in, X_gk_in,
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
                           const Vector< MultiFab* >& p_g_in,
                           const Vector< MultiFab* >& T_g_in,
                           const Vector< MultiFab* >& X_gk_in,
                           const Real time,
                           F1 WeightFunc)
{
  if (m_reaction_rates_type == ReactionRatesType::RRatesUser) {
    mfix_calc_chem_txfr(chem_txfr, ep_g_in, ro_g_in, vel_g_in, p_g_in, T_g_in, X_gk_in,
        time, WeightFunc, HeterogeneousRatesUser(), HomogeneousRatesUser());
  } else {
    amrex::Abort("Invalid Reaction Rates Type.");
  }
}

template <typename F1, typename F2, typename F3>
void 
mfix::mfix_calc_chem_txfr (const Vector< MultiFab* >& chem_txfr,
                           const Vector< MultiFab* >& ep_g_in,
                           const Vector< MultiFab* >& ro_g_in,
                           const Vector< MultiFab* >& vel_g_in,
                           const Vector< MultiFab* >& p_g_in,
                           const Vector< MultiFab* >& T_g_in,
                           const Vector< MultiFab* >& X_gk_in,
                           const Real time,
                           F1 WeightFunc,
                           F2 HeterogeneousRatesFunc,
                           F3 HomogeneousRatesFunc)
{
  using PairIndex = MFIXParticleContainer::PairIndex;
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

  const Real strttime = ParallelDescriptor::second();

  BL_PROFILE("mfix::mfix_calc_chem_txfr()");

  // Solid species data
  const int nspecies_s = solids.nspecies;

  int* p_species_id_s = solids.d_species_id.data();
  Real* p_MW_sn = solids.d_MW_sn0.data();

  // Fluid species data
  const int nspecies_g = fluid.nspecies;

  int* p_species_id_g = fluid.d_species_id.data();
  Real* p_MW_gk = fluid.d_MW_gk0.data();

  // Fluid enthalpy data
  const Real T_ref = fluid.T_ref;

  // Reactions data
  const int nreactions = reactions.nreactions;

  int* p_types = reactions.d_types.data();
  const int** p_phases = reactions.d_phases.data();
  int* p_nphases = reactions.d_nphases.data();

  int* p_nreactants = reactions.d_nreactants.data();
  const int** p_reactants_id = reactions.d_reactants_id.data();
  const Real** p_reactants_coeffs = reactions.d_reactants_coeffs.data();
  const int** p_reactants_phases = reactions.d_reactants_phases.data();

  int* p_nproducts = reactions.d_nproducts.data();
  const int** p_products_id = reactions.d_products_id.data();
  const Real** p_products_coeffs = reactions.d_products_coeffs.data();
  const int** p_products_phases = reactions.d_products_phases.data();

  // Solid phase integer ID
  constexpr int Solid = ChemicalReaction::CHEMICALPHASE::Solid;
  constexpr int Fluid = ChemicalReaction::CHEMICALPHASE::Fluid;
  constexpr int Heterogeneous = ChemicalReaction::REACTIONTYPE::Heterogeneous;
  constexpr int Homogeneous   = ChemicalReaction::REACTIONTYPE::Homogeneous;
  const int InvalidIdx = -1; //TODO define this somewhere else

  // Particles SoA starting indexes for mass fractions and rate of formations
  const int idx_X_sn       = (pc->m_runtimeRealData).X_sn;
  const int idx_mass_sn_txfr = (pc->m_runtimeRealData).mass_sn_txfr;
  const int idx_vel_s_txfr = (pc->m_runtimeRealData).vel_s_txfr;
  const int idx_h_s_txfr   = (pc->m_runtimeRealData).h_s_txfr;

  ChemTransfer chem_txfr_idxs(fluid.nspecies, reactions.nreactions);
  const int idx_ro_gk_txfr = chem_txfr_idxs.ro_gk_txfr;
  const int idx_vel_g_txfr = chem_txfr_idxs.vel_g_txfr;
  const int idx_h_g_txfr   = chem_txfr_idxs.h_g_txfr;

  auto& fluid_parms = *fluid.parameters;
  auto& solids_parms = *solids.parameters;
  auto& reactions_parms = *reactions.parameters;

  const Geometry& gm = Geom(0);

  // **************************************************************************
  // Compute particles densities transfer rates
  // **************************************************************************
  // Extrapolate velocity Dirichlet bc's to ghost cells
  int extrap_dir_bcs = 1;

  mfix_set_velocity_bcs(time, vel_g_in, extrap_dir_bcs);

  // We copy the value inside the domain to the outside to avoid
  // unphysical volume fractions.
  const int dir_bc_in = 2;
  mfix_set_epg_bcs(ep_g_in, dir_bc_in);

  // Set boundary conditions just in case
  mfix_set_density_bcs(time, ro_g_in);
  mfix_set_temperature_bcs(time, T_g_in);
  mfix_set_species_bcs(time, X_gk_in);

  const auto cell_dx = gm.CellSizeArray();
  const auto reg_cell_vol = cell_dx[0]*cell_dx[1]*cell_dx[2];

  Vector< MultiFab* > tmp_eps(nlev, nullptr);
  Vector< MultiFab* > chem_txfr_ptr(nlev, nullptr);

  Vector<std::map<PairIndex, Gpu::DeviceVector<Real>>> G_ro_gk_deposition(nlev);
  Vector<std::map<PairIndex, Gpu::DeviceVector<Real>>> G_h_g_deposition(nlev);

  for (int lev = 0; lev < nlev && (DEM::solve || PIC::solve); lev++) {

    MultiFab* interp_ptr;

    // Set to zero the chem_txfr multifab
    chem_txfr[lev]->setVal(0);

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    EB_set_covered(*ep_g_in[lev], 0, 1, 1, covered_val);
    EB_set_covered(*ro_g_in[lev], 0, 1, 1, covered_val);
    EB_set_covered(*vel_g_in[lev], 0, 3, 1, covered_val);
    EB_set_covered(*T_g_in[lev], 0, 1, 1, covered_val);
    EB_set_covered(*X_gk_in[lev], 0, fluid.nspecies, 1, covered_val);
    EB_set_covered(*chem_txfr[lev], 0, chem_txfr_idxs.count, 1, covered_val);

    const int interp_ng = 1;  // Only one layer needed for interpolation
    const int interp_comp = 7 + nspecies_g; // 3 vel_g + ep_g + ro_g + p_g + T_g + X_g

    MultiFab p_nd(p_g_in[lev]->boxArray(), dmap[lev], 1, interp_ng);
    p_nd.setVal(0.);

    MultiFab::Copy(p_nd, *m_leveldata[lev]->p_g, 0, 0, 1, interp_ng);
    MultiFab::Add (p_nd, *m_leveldata[lev]->p0_g, 0, 0, 1, interp_ng);

    MultiFab p_cc(ep_g_in[lev]->boxArray(), dmap[lev], 1, interp_ng);
    amrex::average_node_to_cellcenter(p_cc, 0, p_nd, 0, 1);

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

      // Copy perturbational pressure
      MultiFab::Copy(*interp_ptr, p_cc, 0, 5, 1, interp_ng);

      // Copy temperature
      MultiFab::Copy(*interp_ptr, *T_g_in[lev], 0, 6, 1, interp_ng);

      // Copy X_gk
      MultiFab::Copy(*interp_ptr, *X_gk_in[lev], 0, 7, nspecies_g, interp_ng);

      // Set chem_txfr_ptr
      chem_txfr_ptr[lev] = chem_txfr[lev];
    }
    else
    {
      const BoxArray&            pba = pc->ParticleBoxArray(lev);
      const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(pba, pdm, interp_comp, interp_ng, MFInfo(),
                                *particle_ebfactory[lev]);

      // Copy velocity
      interp_ptr->ParallelCopy(*vel_g_in[lev], 0, 0, 3, interp_ng, interp_ng);

      // Copy volume fraction
      interp_ptr->ParallelCopy(*ep_g_in[lev], 0, 3, 1, interp_ng, interp_ng);

      // Copy density
      interp_ptr->ParallelCopy(*ro_g_in[lev], 0, 4, 1, interp_ng, interp_ng);

      // Copy pressure
      interp_ptr->ParallelCopy(p_cc, 0, 5, 1, interp_ng, interp_ng);

      // Copy temperature
      interp_ptr->ParallelCopy(*T_g_in[lev], 0, 6, 1, interp_ng, interp_ng);

      // Copy X_gk
      interp_ptr->ParallelCopy(*X_gk_in[lev], 0, 7, fluid.nspecies, interp_ng, interp_ng);

      // Set chem_txfr_ptr
      chem_txfr_ptr[lev] = new MultiFab(pba, pdm, interp_comp, interp_ng, MFInfo(),
                                        *particle_ebfactory[lev]);
    }

    // FillBoundary on interpolation MultiFab
    interp_ptr->FillBoundary(geom[lev].periodicity());

    tmp_eps[lev] = (MFHelpers::createFrom(*interp_ptr, 0.0)).release();

    // We must have ghost cells for each FAB so that a particle in one grid can
    // spread its effect to an adjacent grid by first putting the value into
    // ghost cells of its own grid.  The mf->sumBoundary call then adds the
    // value from one grid's ghost cell to another grid's valid region.
    if (chem_txfr_ptr[lev]->nGrow() < 1)
      amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

    // Do the interpolate
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
    const auto& volfrac = factory.getVolFrac();

    // *************************************************************************
    // Homogeneous fluid reactions
    // *************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*interp_ptr,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      // This is to check efficiently if this tile contains any eb stuff
      const EBFArrayBox& flags_fab = static_cast<EBFArrayBox const&>((*interp_ptr)[mfi]);
      const EBCellFlagFab& flags = flags_fab.getEBCellFlagFab();

      if (flags.getType(amrex::grow(bx,0)) != FabType::covered) {
        
        // Get arrays
        const auto& ep_g_array       = interp_ptr->const_array(mfi, 3);
        const auto& ro_g_array       = interp_ptr->const_array(mfi, 4);
        const auto& T_g_array        = interp_ptr->const_array(mfi, 6);
        const auto& X_gk_array       = interp_ptr->const_array(mfi, 7);
        const auto& ro_gk_txfr_array = chem_txfr_ptr[lev]->array(mfi, idx_ro_gk_txfr);
        const auto& h_g_txfr_array   = chem_txfr_ptr[lev]->array(mfi, idx_h_g_txfr);

        auto const& volfrac_arr = volfrac.const_array(mfi);

        amrex::ParallelFor(bx,
          [ro_gk_txfr_array,h_g_txfr_array,HomogeneousRatesFunc,nspecies_g,
           nreactions,p_MW_gk,p_species_id_g,p_reactants_id,p_reactants_coeffs,
           p_reactants_phases,p_products_id,p_products_coeffs,p_products_phases,
           p_nreactants,p_nproducts,InvalidIdx,p_phases,p_nphases,p_types,
           Homogeneous,T_ref,fluid_parms,reactions_parms,reg_cell_vol,
           volfrac_arr,run_on_device,ep_g_array,ro_g_array,T_g_array,X_gk_array]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const Real ep_g  = ep_g_array(i,j,k);
          const Real ro_g  = ro_g_array(i,j,k);
          const Real T_g   = T_g_array(i,j,k);
          const Real vfrac = volfrac_arr(i,j,k);

          const Real fluid_vol = ep_g * vfrac * reg_cell_vol;

          GpuArray<Real,SPECIES::NMAX> X_gk;
          X_gk.fill(0.);

          for (int n_g(0); n_g < nspecies_g; ++n_g) {
            X_gk[n_g] = X_gk_array(i,j,k,n_g);
          }

          // Structure for storing reaction rates in the stencil
          GpuArray<Real,Reactions::NMAX> R_q_homogeneous;
          R_q_homogeneous.fill(0.);

          if (run_on_device) {
            HomogeneousRatesFunc.template operator()<RunOn::Device>(R_q_homogeneous.data(),
                                                                    reactions_parms,
                                                                    fluid_parms, X_gk.data(),
                                                                    ro_g, ep_g);
          } else {
            HomogeneousRatesFunc.template operator()<RunOn::Host>(R_q_homogeneous.data(),
                                                                  reactions_parms,
                                                                  fluid_parms, X_gk.data(),
                                                                  ro_g, ep_g);
          }

          // Total transfer rates
          Real G_m_g_homogeneous(0.);
          Real G_H_g_homogeneous(0.);

          //***************************************************************
          // Loop over fluid species for computing fluid txfr rates
          //***************************************************************
          for (int n_g(0); n_g < nspecies_g; n_g++) {
            Real G_m_gk_homogeneous(0.);

            // Get the ID of the current species n_g
            const int current_species_id = p_species_id_g[n_g];

            // Loop over reactions to compute each contribution
            for (int q(0); q < nreactions; q++) {
              // Do something only if reaction is heterogeneous and contains
              // a solid compound
              if (p_types[q] == Homogeneous &&
                  MFIXfind(p_phases[q], p_nphases[q], Fluid) != InvalidIdx) {

                Real stoc_coeff(0);

                // Add reactant contribution (if any)
                {
                  const int pos = MFIXfind(p_reactants_id[q], p_nreactants[q], current_species_id);

                  if (pos != InvalidIdx) {
                    BL_ASSERT(p_reactants_phases[q][pos] == Fluid);
                    stoc_coeff += p_reactants_coeffs[q][pos];
                  }
                }

                // Add products contribution (if any)
                {
                  const int pos = MFIXfind(p_products_id[q], p_nproducts[q], current_species_id);

                  if (pos != InvalidIdx) {
                    BL_ASSERT(p_products_phases[q][pos] == Fluid);
                    stoc_coeff += p_products_coeffs[q][pos];
                  }
                }

                // Compute particle's species n_s transfer rate for reaction q
                // Note R[q] is in mol/s
                Real G_m_gk_q = stoc_coeff * p_MW_gk[n_g] * R_q_homogeneous[q];

                G_m_gk_homogeneous += G_m_gk_q;

                // Contribution to the particle
                const Real h_gk_T_g = run_on_device ?
                  fluid_parms.calc_h_gk<RunOn::Device>(T_g,n_g) :
                  fluid_parms.calc_h_gk<RunOn::Host>(T_g,n_g);

                G_H_g_homogeneous += h_gk_T_g * G_m_gk_q;
              }
            }

            G_m_g_homogeneous += G_m_gk_homogeneous;
            ro_gk_txfr_array(i,j,k,n_g) += G_m_gk_homogeneous / fluid_vol;
          }

          // Check that (total_ro_g_txfr_homogeneous) = 0
          //BL_ASSERT(std::abs(G_ro_g_homogeneous) < 1.e-15);

          h_g_txfr_array(i,j,k) += G_H_g_homogeneous / fluid_vol;
        });
      } // if entire FAB not covered
    } // mfi

    // *************************************************************************
    // Heterogeneous fluid reactions
    // *************************************************************************
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

      // *******************************************************************
      // Reaction rates for each particle
      // *******************************************************************
      Box bx = pti.tilebox();

      // This is to check efficiently if this tile contains any eb stuff
      const EBFArrayBox& interp_fab = static_cast<EBFArrayBox const&>((*interp_ptr)[pti]);

      const EBCellFlagFab& flags = interp_fab.getEBCellFlagFab();

      if (flags.getType(amrex::grow(bx,0)) != FabType::covered) {

        G_ro_gk_deposition[lev][index] = Gpu::DeviceVector<Real>(np*nspecies_g);
        Real* G_ro_gk_ptr = G_ro_gk_deposition[lev][index].dataPtr();

        G_h_g_deposition[lev][index] = Gpu::DeviceVector<Real>(np);
        Real* G_h_g_ptr = G_h_g_deposition[lev][index].dataPtr();

        // Get arrays
        const auto& ep_g_array = interp_ptr->const_array(pti, 3);
        const auto& interp_array = interp_ptr->const_array(pti);

        const auto& flags_array = flags.const_array();
        auto const& volfrac_arr = volfrac.const_array(pti);

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

        const int dem_solve = DEM::solve;

        amrex::ParallelFor(np,
          [np,pstruct,p_realarray,interp_array,HeterogeneousRatesFunc,
           HomogeneousRatesFunc,G_ro_gk_ptr,G_h_g_ptr,plo,dxi,ptile_data,
           nspecies_g,nspecies_s,nreactions,interp_comp,Solid,p_MW_sn,idx_X_sn,
           idx_mass_sn_txfr,idx_vel_s_txfr,idx_h_s_txfr,p_MW_gk,p_species_id_s,
           p_species_id_g,p_reactants_id,p_reactants_coeffs,p_reactants_phases,
           p_products_id,p_products_coeffs,p_products_phases,p_nreactants,
           p_nproducts,InvalidIdx,p_phases,p_nphases,p_types,Heterogeneous,
           Homogeneous,T_ref,fluid_parms,solids_parms,grown_bx_is_regular,
           ccent_fab,bcent_fab,apx_fab,apy_fab,apz_fab,flags_array,dx,
           reactions_parms,reg_cell_vol,volfrac_arr,dem_solve,run_on_device,
           ep_g_array]
          AMREX_GPU_DEVICE (int p_id) noexcept
        {
          auto& particle = pstruct[p_id];

          // Cell containing particle centroid
          int ip = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
          int jp = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
          int kp = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));
          
          if (!grown_bx_is_regular and flags_array(ip,jp,kp).isCovered()) {

            // Particle mass calculation
            for (int n_s(0); n_s < nspecies_s; n_s++) {
              // Update specie mass formation rate
              ptile_data.m_runtime_rdata[idx_mass_sn_txfr + n_s][p_id] = 0;
            }

            if (dem_solve) {
              ptile_data.m_runtime_rdata[idx_vel_s_txfr + 0][p_id] = 0;
              ptile_data.m_runtime_rdata[idx_vel_s_txfr + 1][p_id] = 0;
              ptile_data.m_runtime_rdata[idx_vel_s_txfr + 2][p_id] = 0;
            }

            ptile_data.m_runtime_rdata[idx_h_s_txfr][p_id] = 0;

          } else {

            GpuArray<Real,SPECIES::NMAX> X_sn;

            for (int n_s(0); n_s < nspecies_s; n_s++) {
              const int idx = idx_X_sn + n_s;
              X_sn[n_s] = ptile_data.m_runtime_rdata[idx][p_id];
            }

            const Real T_p = p_realarray[SoArealData::temperature][p_id];

            GpuArray<Real,7+SPECIES::NMAX> interp_loc; // vel_g, ep_g, ro_g, p_g, T_g, X_gk
            interp_loc.fill(0.);

            GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

            if (grown_bx_is_regular) {

              trilinear_interp(particle.pos(), ip, jp, kp, weights,
                               interp_loc.data(), interp_array, plo, dxi,
                               interp_comp);

            } else if (!flags_array(ip,jp,kp).isCovered()) {

              // Cut or regular cell and none of the cells in the stencil is
              // covered (Note we can't assume regular cell has no covered
              // cells in the stencil because of the diagonal case)

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
            }

            // Common part from here
            const RealVect vel_g(interp_loc[0], interp_loc[1], interp_loc[2]);
            const Real ep_g = interp_loc[3];
            const Real ro_g = interp_loc[4];
            const Real p_g  = interp_loc[5] + 101325; // TODO TODO fix this
            const Real T_g  = interp_loc[6];

            Real* X_gk = &interp_loc[7];

            const Real ep_s = 1. - ep_g;
            const Real ro_s = p_realarray[SoArealData::density][p_id];

            // NOTE: this quanity is not interpolated to particle position since
            // it is used to divide the fluid chem_txfr quantities to obtain
            // volume divided mass and enthalpy transfers
            const Real fluid_vol = ep_g_array(ip,jp,kp) * volfrac_arr(ip,jp,kp) * reg_cell_vol;

            GpuArray<Real,Reactions::NMAX> R_q_heterogeneous;
            R_q_heterogeneous.fill(0.);

            const Real T_s = p_realarray[SoArealData::temperature][p_id];
            const Real DP  = 2. * p_realarray[SoArealData::radius][p_id];
            const RealVect vel_s(p_realarray[SoArealData::velx][p_id],
                                 p_realarray[SoArealData::vely][p_id],
                                 p_realarray[SoArealData::velz][p_id]);

            if (run_on_device) {
              HeterogeneousRatesFunc.template operator()<RunOn::Device>(R_q_heterogeneous.data(),
                                                                        reactions_parms,
                                                                        solids_parms, X_sn.data(),
                                                                        ro_s, ep_s, T_s,
                                                                        vel_s, fluid_parms,
                                                                        X_gk, ro_g, ep_g, T_g,
                                                                        vel_g, DP, p_g);
            } else {
              HeterogeneousRatesFunc.template operator()<RunOn::Host>(R_q_heterogeneous.data(),
                                                                      reactions_parms,
                                                                      solids_parms, X_sn.data(),
                                                                      ro_s, ep_s, T_s,
                                                                      vel_s, fluid_parms,
                                                                      X_gk, ro_g, ep_g, T_g,
                                                                      vel_g, DP, p_g);
            }

            GpuArray<Real,Reactions::NMAX> R_q_homogeneous;
            R_q_homogeneous.fill(0.);

            if (run_on_device) {
              HomogeneousRatesFunc.template operator()<RunOn::Device>(R_q_homogeneous.data(),
                                                                      reactions_parms,
                                                                      solids_parms, X_sn.data(),
                                                                      ro_s, ep_s);
            } else {
              HomogeneousRatesFunc.template operator()<RunOn::Host>(R_q_homogeneous.data(),
                                                                    reactions_parms,
                                                                    solids_parms, X_sn.data(),
                                                                    ro_s, ep_s);
            }

            // Total transfer rates
            Real G_m_g_heterogeneous(0.);
            Real G_H_g_heterogeneous(0.);

            Real G_m_p_heterogeneous(0.);
            Real G_m_p_homogeneous(0.);

            //***************************************************************
            // Initialize to zero
            //***************************************************************
            for (int n_s(0); n_s < nspecies_s; n_s++) {
              // Initially set species n_s density transfer rate to zero
              ptile_data.m_runtime_rdata[idx_mass_sn_txfr + n_s][p_id] = 0;
            }

            if (dem_solve) {
              ptile_data.m_runtime_rdata[idx_vel_s_txfr + 0][p_id] = 0.;
              ptile_data.m_runtime_rdata[idx_vel_s_txfr + 1][p_id] = 0.;
              ptile_data.m_runtime_rdata[idx_vel_s_txfr + 2][p_id] = 0.;
            }

            // Initially set particle energy transfer rate to zero
            ptile_data.m_runtime_rdata[idx_h_s_txfr][p_id] = 0.;

            //***************************************************************
            // Loop over particle's species for computing particle txfr rates
            //***************************************************************
            for (int n_s(0); n_s < nspecies_s; n_s++) {
              Real G_m_pk_heterogeneous(0.);
              Real G_m_pk_homogeneous(0.);

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
                  Real G_m_pk_q = stoc_coeff * p_MW_sn[n_s] * R_q_heterogeneous[q];

                  G_m_pk_heterogeneous += G_m_pk_q;
                }

                // Do something only if reaction is heterogeneous and contains
                // a solid compound
                if (p_types[q] == Homogeneous &&
                    MFIXfind(p_phases[q], p_nphases[q], Solid) != InvalidIdx) {

                  Real stoc_coeff(0);

                  // Add reactant contribution (if any)
                  {
                    const int pos = MFIXfind(p_reactants_id[q], p_nreactants[q], current_species_id);

                    if (pos != InvalidIdx) {
                      BL_ASSERT(p_reactants_phases[q][pos] == Solid);
                      stoc_coeff += p_reactants_coeffs[q][pos];
                    }
                  }

                  // Add products contribution (if any)
                  {
                    const int pos = MFIXfind(p_products_id[q], p_nproducts[q], current_species_id);

                    if (pos != InvalidIdx) {
                      BL_ASSERT(p_products_phases[q][pos] == Solid);
                      stoc_coeff += p_products_coeffs[q][pos];
                    }
                  }

                  // Compute particle's species n_s transfer rate for reaction q
                  Real G_m_pk_q = stoc_coeff * p_MW_sn[n_s] * R_q_homogeneous[q];

                  G_m_pk_homogeneous += G_m_pk_q;
                }
              }

              G_m_p_heterogeneous += G_m_pk_heterogeneous;
              G_m_p_homogeneous   += G_m_pk_homogeneous;

              // Update global variable
              ptile_data.m_runtime_rdata[idx_mass_sn_txfr + n_s][p_id] =
                G_m_pk_heterogeneous + G_m_pk_homogeneous;
            }

            //***************************************************************
            // Initialize to zero
            //***************************************************************
            for (int n_g(0); n_g < nspecies_g; n_g++) {
              // Initially set fluid species n_g density transfer rate to zero
              G_ro_gk_ptr[n_g*np + p_id] = 0;
            }

            // Initially set fluid energy transfer rate to zero
            G_h_g_ptr[p_id] = 0;

            //***************************************************************
            // Loop over fluid species for computing fluid txfr rates
            //***************************************************************
            for (int n_g(0); n_g < nspecies_g; n_g++) {
              Real G_m_gk_heterogeneous(0.);

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

                  // Compute fluid species n_g transfer rate for reaction q
                  Real G_m_gk_q = stoc_coeff * p_MW_gk[n_g] * R_q_heterogeneous[q];

                  G_m_gk_heterogeneous += G_m_gk_q;

                  const Real h_gk_T_p = run_on_device ?
                    fluid_parms.calc_h_gk<RunOn::Device>(T_p,n_g) :
                    fluid_parms.calc_h_gk<RunOn::Host>(T_p,n_g);

                  const Real h_gk_T_g = run_on_device ?
                    fluid_parms.calc_h_gk<RunOn::Device>(T_g,n_g) :
                    fluid_parms.calc_h_gk<RunOn::Host>(T_g,n_g);

                  const Real G_H_pk_q = h_gk_T_g * amrex::min(0., G_m_gk_q);
                  const Real G_H_gk_q = h_gk_T_p * amrex::max(0., G_m_gk_q);

                  G_H_g_heterogeneous += (G_H_pk_q + G_H_gk_q);
                }
              }

              G_m_g_heterogeneous += G_m_gk_heterogeneous;

              // Update global variable
              G_ro_gk_ptr[n_g*np + p_id] = G_m_gk_heterogeneous / fluid_vol;
            }

            // Check that (total_ro_g_txfr + total_ro_p_txfr) = 0
            BL_ASSERT(std::abs(G_m_p_heterogeneous + G_m_g_heterogeneous) < 1.e-15);
            BL_ASSERT(std::abs(G_m_p_homogeneous) < 1.e-15);

            //***************************************************************
            // Update particle linear momentum and energy transfer
            //***************************************************************
            if (dem_solve) {
              const Real coeff = amrex::max(0., G_m_p_heterogeneous);

              // NOTE: total_ro_p_txfr is computed from interpolated quantities
              // and vel_g is also interpolated to particle position
              ptile_data.m_runtime_rdata[idx_vel_s_txfr+0][p_id] = coeff*vel_g[0];
              ptile_data.m_runtime_rdata[idx_vel_s_txfr+1][p_id] = coeff*vel_g[1];
              ptile_data.m_runtime_rdata[idx_vel_s_txfr+2][p_id] = coeff*vel_g[2];
            }

            // Write the result in the enthalpy transfer space
            ptile_data.m_runtime_rdata[idx_h_s_txfr][p_id] = -G_H_g_heterogeneous;
            G_h_g_ptr[p_id] = G_H_g_heterogeneous / fluid_vol;
          }
        });

      } // if entire FAB not covered
    }

    delete interp_ptr;
  }

  Vector< const FabArray<EBCellFlagFab>* > flags(nlev, nullptr);
  Vector< const MultiFab* > volfrac(nlev, nullptr);
  Vector< EBFArrayBoxFactory* > crse_factory(nlev, nullptr);

  // Deposition algorithm
  for (int lev = 0; lev < nlev && (DEM::solve || PIC::solve); lev++) {

    // Do the interpolate
    const auto dxi_array = geom[lev].InvCellSizeArray();
    const auto dx_array  = geom[lev].CellSizeArray();
    const auto plo_array = geom[lev].ProbLoArray();

    // Use level 0 to define the EB factory. If we are not on level 0
    // then create a copy of the coarse factory to use.
    if (lev == 0) {
      flags[lev]   = &(particle_ebfactory[lev]->getMultiEBCellFlagFab());
      volfrac[lev] = &(particle_ebfactory[lev]->getVolFrac());
    } else {
      Vector<int> ngrow = {1,1,1};

      crse_factory[lev] = (makeEBFabFactory(gm, chem_txfr_ptr[lev]->boxArray(),
                                            chem_txfr_ptr[lev]->DistributionMap(),
                                            ngrow, EBSupport::volume)).release();

      flags[lev]   = &(crse_factory[lev]->getMultiEBCellFlagFab());
      volfrac[lev] = &(crse_factory[lev]->getVolFrac());
    }

    for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {

      PairIndex index(pti.index(), pti.LocalTileIndex());

      auto& particles = pti.GetArrayOfStructs();
      MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const int np = particles.size();

      // ******************************************************************
      // Deposit the interphase chemical transfer quantities to the grid
      // Density transfer:
      // Momentum transfer:
      // Energy transfer
      // ******************************************************************

      FArrayBox& eps_fab  = (*tmp_eps[lev])[pti];
      FArrayBox& chem_txfr_fab  = (*chem_txfr_ptr[lev])[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags[lev])[pti].getType(box) != FabType::covered) {

        Real* G_ro_gk_ptr = G_ro_gk_deposition[lev][index].dataPtr();
        Real* G_h_g_ptr = G_h_g_deposition[lev][index].dataPtr();

        auto ro_gk_txfr_array = chem_txfr_fab.array(idx_ro_gk_txfr);
        auto vel_g_txfr_array = chem_txfr_fab.array(idx_vel_g_txfr);
        auto h_g_txfr_array   = chem_txfr_fab.array(idx_h_g_txfr);
        auto vol_arr          = eps_fab.array();

        const auto& flagsarr = (*flags[lev])[pti].const_array();
        const auto& vfrac    = (*volfrac[lev])[pti].const_array();

        const amrex::Real deposition_scale_factor = m_deposition_scale_factor;

#ifdef _OPENMP
        FArrayBox local_fab_chem_txfr;
        FArrayBox local_fab_vol;

        Box tile_box = box;

        if (Gpu::notInLaunchRegion()) {
          tile_box.grow(chem_txfr_ptr[lev]->nGrow());

          local_fab_chem_txfr.resize(tile_box, chem_txfr_ptr[lev]->nComp());
          local_fab_chem_txfr.setVal<RunOn::Host>(0.0);

          ro_gk_txfr_array = local_fab_chem_txfr.array(idx_ro_gk_txfr);
          vel_g_txfr_array = local_fab_chem_txfr.array(idx_vel_g_txfr);
          h_g_txfr_array   = local_fab_chem_txfr.array(idx_h_g_txfr);

          local_fab_vol.resize(tile_box, tmp_eps[lev]->nComp());
          local_fab_vol.setVal<RunOn::Host>(0.0);
          vol_arr = local_fab_vol.array();
        }
#endif

        const long nrp = pti.numParticles();

        amrex::ParallelFor(nrp,
          [nrp,pstruct,p_realarray,plo_array,dx_array,dxi_array,vfrac,flagsarr,
           deposition_scale_factor,WeightFunc,ro_gk_txfr_array,vel_g_txfr_array,
           h_g_txfr_array,idx_mass_sn_txfr,G_ro_gk_ptr,G_h_g_ptr,nspecies_g,
           run_on_device,vol_arr,reg_cell_vol,local_cg_dem=DEM::cg_dem]
          AMREX_GPU_DEVICE (int p_id) noexcept
        {
          const auto& p = pstruct[p_id];

          int i(0);
          int j(0);
          int k(0);

          const Real statwt = p_realarray[SoArealData::statwt][p_id];

          GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

          WeightFunc(plo_array, dx_array, dxi_array, flagsarr, p.pos(),
            p_realarray[SoArealData::radius][p_id], i, j, k, weights,
            deposition_scale_factor);

          GpuArray<Real,SPECIES::NMAX> G_ro_gk_heterogeneous;
          G_ro_gk_heterogeneous.fill(0.);

          Real G_ro_g_heterogeneous(0.);

          for (int n_g(0); n_g < nspecies_g; ++n_g) {
            G_ro_gk_heterogeneous[n_g] = G_ro_gk_ptr[n_g*nrp + p_id];
            G_ro_g_heterogeneous += G_ro_gk_heterogeneous[n_g];
          }

          const Real G_h_g_heterogeneous = G_h_g_ptr[p_id];

          const Real coeff = amrex::max(0., G_ro_g_heterogeneous);

          Real pvol = statwt * p_realarray[SoArealData::volume][p_id] / reg_cell_vol;

          if (local_cg_dem){
            pvol = pvol / p_realarray[SoArealData::statwt][p_id];
          }

          const RealVect vel_p(p_realarray[SoArealData::velx][p_id],
                               p_realarray[SoArealData::vely][p_id],
                               p_realarray[SoArealData::velz][p_id]);

          // Deposition
          for (int ii = -1; ii <= 0; ++ii) {
            for (int jj = -1; jj <= 0; ++jj) {
              for (int kk = -1; kk <= 0; ++kk) {
                if (! flagsarr(i+ii,j+jj,k+kk).isCovered()) {

                  amrex::Real weight_vol = weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk);

                  amrex::Gpu::Atomic::Add(&vol_arr(i+ii,j+jj,k+kk), weight_vol*pvol);

                  for (int n_g(0); n_g < nspecies_g; n_g++) {
                    // Deposition of Rrates
                    Gpu::Atomic::Add(&ro_gk_txfr_array(i+ii,j+jj,k+kk,n_g),
                                     weight_vol*G_ro_gk_heterogeneous[n_g]);
                  }

                  if (coeff > 0.) {
                    // Deposition of linear momentum
                    Gpu::Atomic::Add(&vel_g_txfr_array(i+ii,j+jj,k+kk,0), weight_vol*coeff*vel_p[0]);
                    Gpu::Atomic::Add(&vel_g_txfr_array(i+ii,j+jj,k+kk,1), weight_vol*coeff*vel_p[1]);
                    Gpu::Atomic::Add(&vel_g_txfr_array(i+ii,j+jj,k+kk,2), weight_vol*coeff*vel_p[2]);
                  }

                  Gpu::Atomic::Add(&h_g_txfr_array(i+ii,j+jj,k+kk), weight_vol*G_h_g_heterogeneous);
                }
              }
            }
          }
        });

#ifdef _OPENMP
        if (Gpu::notInLaunchRegion()) {
          chem_txfr_fab.atomicAdd<RunOn::Host>(local_fab_chem_txfr, tile_box, tile_box, 0, 0, chem_txfr_ptr[lev]->nComp());
          eps_fab.atomicAdd<RunOn::Host>(local_fab_vol, tile_box, tile_box, 0, 0, tmp_eps[lev]->nComp());

        }
#endif
      }

    } // pti

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
    chem_txfr_ptr[lev]->setBndry(0.0);

    // Sum grid boundaries then fill with correct ghost values.
    tmp_eps[lev]->SumBoundary(gm.periodicity());
    tmp_eps[lev]->FillBoundary(gm.periodicity());

    // Move excessive solids volume from small cells to neighboring cells.
    // Note that we don't change tmp_eps but use the redistribution of
    // particle volume to determine how to redistribute the drag forces.
    mfix_redistribute_deposition(lev, *tmp_eps[lev], *chem_txfr_ptr[lev], volfrac[lev], flags[lev],
                                 mfix::m_max_solids_volume_fraction);
  }

  // This might not need to exist on all levels. Maybe only level 0.
  for (int lev(0); lev < nlev; ++lev) {

    if (lev != 0 && crse_factory[lev] != nullptr)
      delete crse_factory[lev];

    delete tmp_eps[lev];
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
      chem_txfr[lev]->setVal(0);
      amrex::InterpFromCoarseLevel(*chem_txfr[lev], time, *chem_txfr_ptr[lev-1],
                                   0, 0, 1, Geom(lev-1), Geom(lev),
                                   cphysbc, 0, fphysbc, 0,
                                   ref_ratio[0], mapper,
                                   bcs, 0);
    }
  }

  // If mf_to_be_filled is not defined on the particle_box_array, then we need
  // to copy here from txfr_ptr into mf_to_be_filled. I believe that we don't
  // need any information in ghost cells so we don't copy those.

  if (chem_txfr_ptr[0] != chem_txfr[0]) {
    chem_txfr[0]->ParallelCopy(*chem_txfr_ptr[0], 0, 0, chem_txfr[0]->nComp());
  }

  for (int lev = 0; lev < nlev; lev++) {
    if (chem_txfr_ptr[lev] != m_leveldata[lev]->chem_txfr)
      delete chem_txfr_ptr[lev];
  }

  if (m_verbose > 1) {
    Real stoptime = ParallelDescriptor::second() - strttime;

    ParallelDescriptor::ReduceRealMax(stoptime, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "MFIXParticleContainer::TrilinearChemDepositionFluid time: " << stoptime << '\n';
  }

  // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
  for (int lev = 0; lev < nlev; lev++)
    chem_txfr[lev]->FillBoundary(geom[lev].periodicity());
}
