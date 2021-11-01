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

  auto& reactions_parms = *reactions.parameters;

  // ************************************************************************
  // Reset chem_txfr values to 0 before updating them
  // ************************************************************************
  for (int lev = 0; lev < nlev; lev++) {
    chem_txfr[lev]->setVal(0);
  }


  Vector< MultiFab* > chem_txfr_ptr(nlev, nullptr);

  for (int lev = 0; lev < nlev; lev++)
  {
    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    if (lev == 0 && OnSameGrids)
    {
      // If we are already working with the internal mf defined on the
      // particle_box_array, then we just work with this.
      chem_txfr_ptr[lev] = chem_txfr[lev];
    }
    else if (lev == 0 && (! OnSameGrids))
    {
      // If beta_mf is not defined on the particle_box_array, then we need
      // to make a temporary here and copy into beta_mf at the end.
      chem_txfr_ptr[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                        pc->ParticleDistributionMap(lev),
                                        chem_txfr[lev]->nComp(),
                                        chem_txfr[lev]->nGrow());
    }
    else
    {
      // If lev > 0 we make a temporary at the coarse resolution
      BoxArray ba_crse(amrex::coarsen(pc->ParticleBoxArray(lev),
                                      this->m_gdb->refRatio(0)));

      chem_txfr_ptr[lev] = new MultiFab(ba_crse, pc->ParticleDistributionMap(lev),
                                        chem_txfr[lev]->nComp(), 1);
    }

    // We must have ghost cells for each FAB so that a particle in one grid can
    // spread its effect to an adjacent grid by first putting the value into
    // ghost cells of its own grid.  The mf->sumBoundary call then adds the
    // value from one grid's ghost cell to another grid's valid region.
    if (chem_txfr_ptr[lev]->nGrow() < 1)
      amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

    chem_txfr_ptr[lev]->setVal(0.0, 0, chem_txfr[lev]->nComp(), chem_txfr_ptr[lev]->nGrow());
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

#if 0
  for (int lev = 0; lev < nlev && fluid.solve; lev++) {

    auto& fluid_parms = *fluid.parameters;

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ep_g_in[lev]->Factory());
    const auto& volfrac = factory.getVolFrac();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      // This is to check efficiently if this tile contains any eb stuff
      const EBFArrayBox& flags_fab = static_cast<EBFArrayBox const&>((*ep_g_in[lev])[mfi]);
      const EBCellFlagFab& flags = flags_fab.getEBCellFlagFab();

      if (flags.getType(amrex::grow(bx,0)) != FabType::covered) {
        
        // Get arrays
        const auto& rho_gk_txfr_array = chem_txfr_ptr[lev]->array(mfi,idx_ro_gk_txfr);
        const auto& h_g_txfr_array    = chem_txfr_ptr[lev]->array(mfi,idx_h_g_txfr);
        const auto& ep_g_array        = ep_g_in[lev]->const_array(mfi);
        const auto& ro_g_array        = ro_g_in[lev]->const_array(mfi);
        const auto& T_g_array         = T_g_in[lev]->const_array(mfi);
        const auto& X_gk_array        = X_gk_in[lev]->const_array(mfi);

        auto const& volfrac_arr = volfrac.const_array(mfi);

        amrex::ParallelFor(bx,
          [rho_gk_txfr_array,h_g_txfr_array,ep_g_array,ro_g_array,T_g_array,X_gk_array,
           HomogeneousRatesFunc,nspecies_g,nreactions,p_MW_gk,p_species_id_g,
           p_reactants_id,p_reactants_coeffs,p_reactants_phases,p_products_id,
           p_products_coeffs,p_products_phases,p_nreactants,p_nproducts,InvalidIdx,
           p_phases,p_nphases,p_types,Homogeneous,T_ref,fluid_parms,
           reactions_parms,reg_cell_vol,volfrac_arr,run_on_device]
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
          Real G_rho_g_homogeneous(0.);
          Real G_h_g_homogeneous(0.);

          //***************************************************************
          // Loop over fluid species for computing fluid txfr rates
          //***************************************************************
          for (int n_g(0); n_g < nspecies_g; n_g++) {
            Real G_sk_gg_homogeneous(0.);

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
                Real G_sk_gg_q = stoc_coeff * p_MW_gk[n_g] * R_q_homogeneous[q];

                G_sk_gg_homogeneous += G_sk_gg_q;

                // Contribution to the particle
                const Real h_gk_T_g = run_on_device ?
                  fluid_parms.calc_h_gk<RunOn::Device>(T_g,n_g) :
                  fluid_parms.calc_h_gk<RunOn::Host>(T_g,n_g);

                G_h_g_homogeneous += h_gk_T_g * G_sk_gg_q;
              }
            }

            G_rho_g_homogeneous += G_sk_gg_homogeneous;
            rho_gk_txfr_array(i,j,k,n_g) += G_sk_gg_homogeneous / fluid_vol;
          }

          // Check that (total_ro_g_txfr_homogeneous) = 0
          //BL_ASSERT(std::abs(G_rho_g_homogeneous) < 1.e-15);

          h_g_txfr_array(i,j,k) += G_h_g_homogeneous / fluid_vol;
        });
      } // if entire FAB not covered
    } // mfi
  } // lev
#endif

  for (int lev = 0; lev < nlev && (DEM::solve || PIC::solve); lev++) {
    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    MultiFab* interp_ptr;

    EB_set_covered(*ep_g_in[lev], 0, 1, 1, covered_val);
    EB_set_covered(*ro_g_in[lev], 0, 1, 1, covered_val);
    EB_set_covered(*vel_g_in[lev], 0, 3, 1, covered_val);
    EB_set_covered(*T_g_in[lev], 0, 1, 1, covered_val);
    EB_set_covered(*X_gk_in[lev], 0, fluid.nspecies, 1, covered_val);

    const int interp_ng = 1;  // Only one layer needed for interpolation
    const int interp_comp = 7+nspecies_g; // 3 vel_g + ep_g + ro_g + p_g + T_g + X_gk

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

      // Copy pressure
      interp_ptr->ParallelCopy(p_cc, 0, 5, 1, interp_ng, interp_ng);

      // Copy temperature
      interp_ptr->ParallelCopy(*T_g_in[lev], 0, 6, 1, interp_ng, interp_ng);

      // Copy X_gk
      interp_ptr->ParallelCopy(*X_gk_in[lev], 0, 7, fluid.nspecies, interp_ng, interp_ng);
    }

    // FillBoundary on interpolation MultiFab
    interp_ptr->FillBoundary(geom[lev].periodicity());

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

      if (flags.getType(amrex::grow(bx,0)) != FabType::covered) {
        // Get arrays
        // const auto& ep_g_array   = ep_g_in[lev]->const_array(pti);
        const auto& interp_array = interp_ptr->const_array(pti);

        const auto& ep_g_array   = interp_ptr->const_array(pti,3);
        const auto& flags_array  = flags.const_array();
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
           HomogeneousRatesFunc,G_sk_pg_ptr,G_h_pg_ptr,plo,dxi,ptile_data,
           nspecies_g,nspecies_s,nreactions,interp_comp,Solid,p_MW_sn,idx_X_sn,
           idx_mass_sn_txfr,idx_vel_s_txfr,idx_h_s_txfr,p_MW_gk,p_species_id_s,
           p_species_id_g,p_reactants_id,p_reactants_coeffs,p_reactants_phases,
           p_products_id,p_products_coeffs,p_products_phases,p_nreactants,
           p_nproducts,InvalidIdx,p_phases,p_nphases,p_types,Heterogeneous,
           Homogeneous,T_ref,fluid_parms,solids_parms,
           grown_bx_is_regular,ccent_fab,bcent_fab,apx_fab,apy_fab,apz_fab,
           flags_array,dx,reactions_parms,ep_g_array,reg_cell_vol,
           volfrac_arr,dem_solve,run_on_device]
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
            RealVect vel_g(interp_loc[0], interp_loc[1], interp_loc[2]);
            const Real ep_g = interp_loc[3];
            const Real ro_g = interp_loc[4];
            const Real p_g  = interp_loc[5] + 101325;
            const Real T_g  = interp_loc[6];

            Real* X_gk = &interp_loc[7];

            const Real ep_s = 1. - ep_g;
            const Real ro_s = p_realarray[SoArealData::density][p_id];

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
            Real G_rho_g_heterogeneous(0.);

            Real G_mass_p_heterogeneous(0.);
            Real G_mass_p_homogeneous(0.);

            Real G_h_p_heterogeneous(0.);
            Real G_h_p_homogeneous(0.);


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
              Real G_sn_gp_heterogeneous(0.);
              Real G_sn_pp_homogeneous(0.);

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
                  Real G_sn_gp_q = stoc_coeff * p_MW_sn[n_s] * R_q_heterogeneous[q];

                  G_sn_gp_heterogeneous += G_sn_gp_q;
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
                  Real G_sn_gp_q = stoc_coeff * p_MW_sn[n_s] * R_q_homogeneous[q];

                  G_sn_pp_homogeneous += G_sn_gp_q;
                }
              }

              G_mass_p_heterogeneous += G_sn_gp_heterogeneous;
              G_mass_p_homogeneous   += G_sn_pp_homogeneous;

              //const Real h_sn_T_p = run_on_device ?
              //  solids_parms.calc_h_sn<RunOn::Device>(T_p,n_s) :
              //  solids_parms.calc_h_sn<RunOn::Host>(T_p,n_s);

              // G_h_p_heterogeneous += h_sn_T_p * G_sn_gp_heterogeneous;
              // G_h_p_homogeneous   += h_sn_T_p * G_sn_pp_homogeneous;

              // Update global variable
              ptile_data.m_runtime_rdata[idx_mass_sn_txfr + n_s][p_id] =
                G_sn_gp_heterogeneous + G_sn_pp_homogeneous;
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
              Real G_sk_pg_heterogeneous(0.);

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
                  Real G_sk_pg_q = stoc_coeff * p_MW_gk[n_g] * R_q_heterogeneous[q];

                  G_sk_pg_heterogeneous += G_sk_pg_q;

                  // Contribution to the particle
                  const Real h_gk_T_p = run_on_device ?
                    fluid_parms.calc_h_gk<RunOn::Device>(T_p,n_g) :
                    fluid_parms.calc_h_gk<RunOn::Host>(T_p,n_g);

                  const Real h_gk_T_g = run_on_device ?
                    fluid_parms.calc_h_gk<RunOn::Device>(T_g,n_g) :
                    fluid_parms.calc_h_gk<RunOn::Host>(T_g,n_g);

                  G_h_p_heterogeneous += h_gk_T_p * G_sk_pg_q;
                  G_h_p_heterogeneous += amrex::min(0., G_sk_pg_q) * (h_gk_T_g - h_gk_T_p);

                }
              }

              G_rho_g_heterogeneous += G_sk_pg_heterogeneous;

              // Update global variable
              G_sk_pg_ptr[n_g*np + p_id] = G_sk_pg_heterogeneous / fluid_vol;
            }

            // Check that (total_ro_g_txfr + total_ro_p_txfr) = 0
            //BL_ASSERT(std::abs(G_rho_p_heterogeneous + G_rho_g_heterogeneous) < 1.e-15);
            //BL_ASSERT(std::abs(G_rho_p_homogeneous) < 1.e-15);

            //***************************************************************
            // Update particle linear momentum and energy transfer
            //***************************************************************
            if (dem_solve) {
              const Real coeff = amrex::max(0., G_mass_p_heterogeneous);

              // NOTE: total_ro_p_txfr is computed from interpolated quantities
              // and vel_g is also interpolated to particle position
              ptile_data.m_runtime_rdata[idx_vel_s_txfr+0][p_id] = coeff*vel_g[0];
              ptile_data.m_runtime_rdata[idx_vel_s_txfr+1][p_id] = coeff*vel_g[1];
              ptile_data.m_runtime_rdata[idx_vel_s_txfr+2][p_id] = coeff*vel_g[2];
            }

            // Write the result in the enthalpy transfer space
            ptile_data.m_runtime_rdata[idx_h_s_txfr][p_id] = G_h_p_heterogeneous + G_h_p_homogeneous;
            G_h_pg_ptr[p_id] = -G_h_p_heterogeneous / fluid_vol;
          }
        });

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
        auto        rho_gk_txfr_array = fab_chem_txfr.array(idx_ro_gk_txfr);
        auto        vel_g_txfr_array  = fab_chem_txfr.array(idx_vel_g_txfr);
        auto        h_g_txfr_array    = fab_chem_txfr.array(idx_h_g_txfr);
        const auto& flagsarr          = (*deposition_flags)[pti].const_array();
        const auto& vfrac             = (*deposition_volfrac)[pti].const_array();

        const amrex::Real deposition_scale_factor = m_deposition_scale_factor;

#ifdef _OPENMP
        const int ncomp = chem_txfr_ptr[lev]->nComp();
        Box tile_box = box;

        if (Gpu::notInLaunchRegion()) {
          tile_box.grow(chem_txfr_ptr[lev]->nGrow());

          local_fab_chem_txfr.resize(tile_box, ncomp);
          local_fab_chem_txfr.setVal<RunOn::Host>(0.0);

          rho_gk_txfr_array = local_fab_chem_txfr.array(idx_ro_gk_txfr);
          vel_g_txfr_array  = local_fab_chem_txfr.array(idx_vel_g_txfr);
          h_g_txfr_array    = local_fab_chem_txfr.array(idx_h_g_txfr);
        }
#endif

        const long nrp = pti.numParticles();

        amrex::ParallelFor(nrp,
          [nrp,pstruct,p_realarray,plo_array,dx_array,dxi_array,
           vfrac,flagsarr,deposition_scale_factor,WeightFunc,
           rho_gk_txfr_array,vel_g_txfr_array,h_g_txfr_array,idx_mass_sn_txfr,
           G_sk_pg_ptr,G_h_pg_ptr,nspecies_g,run_on_device]
          AMREX_GPU_DEVICE (int p_id) noexcept
        {
          const auto& p = pstruct[p_id];

          int i(0);
          int j(0);
          int k(0);

          GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

          WeightFunc(plo_array, dx_array, dxi_array, flagsarr, p.pos(),
            p_realarray[SoArealData::radius][p_id], i, j, k, weights,
            deposition_scale_factor);

          GpuArray<Real,SPECIES::NMAX> G_sk_pg_heterogeneous;
          G_sk_pg_heterogeneous.fill(0.);

          Real G_rho_g_heterogeneous(0.);

          for (int n_g(0); n_g < nspecies_g; ++n_g) {
            G_sk_pg_heterogeneous[n_g] = G_sk_pg_ptr[n_g*nrp + p_id];
            G_rho_g_heterogeneous += G_sk_pg_heterogeneous[n_g];
          }

          const Real G_h_g_heterogeneous = G_h_pg_ptr[p_id];

          const Real coeff = amrex::max(0., G_rho_g_heterogeneous);

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
                    Gpu::Atomic::Add(&rho_gk_txfr_array(i+ii,j+jj,k+kk,n_g),
                                     weight_vol*G_sk_pg_heterogeneous[n_g]);
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
          fab_chem_txfr.atomicAdd<RunOn::Host>(local_fab_chem_txfr,
              tile_box, tile_box, 0, 0, ncomp);
        }
#endif
      }

    } // pti

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

      chem_txfr[lev]->setVal(0);

      amrex::InterpFromCoarseLevel(*chem_txfr[lev], time,
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

  if (chem_txfr_ptr[0] != chem_txfr[0]) {
    chem_txfr[0]->ParallelCopy(*chem_txfr_ptr[0], 0, 0, chem_txfr[0]->nComp());
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
