#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_shepard_K.H>
#include <mfix_des_drag_K.H>
#include <mfix_des_conv_coeff_K.H>
#include <mfix_des_rrates_K.H>
#include <mfix_mf_helpers.H>
#include <mfix_algorithm.H>

void mfix::mfix_calc_transfer_coeffs (Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in,
                                      Vector< MultiFab* > const& X_gk_in,
                                      Vector< MultiFab* > const& pressure_g_in)
{
  if (m_drag_type == DragType::WenYu) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                              ComputeDragWenYu(DEM::small_number, DEM::large_number, DEM::eps));
  }
  else if (m_drag_type == DragType::Gidaspow) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                              ComputeDragGidaspow(DEM::small_number,DEM::large_number,DEM::eps));
  }
  else if (m_drag_type == DragType::BVK2) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                              ComputeDragBVK2(DEM::small_number,DEM::large_number,DEM::eps));
  }
  else if (m_drag_type == DragType::UserDrag) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                              ComputeDragUser(DEM::small_number,DEM::large_number,DEM::eps));
  }
  else {
    amrex::Abort("Invalid Drag Type.");
  }
}

template <typename F1>
void mfix::mfix_calc_transfer_coeffs (Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in,
                                      Vector< MultiFab* > const& X_gk_in,
                                      Vector< MultiFab* > const& pressure_g_in,
                                      F1 DragFunc)
{
  if (advect_enthalpy)
  {
    if (m_convection_type == ConvectionType::RanzMarshall) {
      mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                                DragFunc,
                                ComputeConvRanzMarshall(DEM::small_number,DEM::large_number,DEM::eps));
    }
    else if (m_convection_type == ConvectionType::Gunn) {
      mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                                DragFunc,
                                ComputeConvGunn(DEM::small_number,DEM::large_number,DEM::eps));
    }
    else if (m_convection_type == ConvectionType::NullConvection) {
      mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                                DragFunc,
                                NullConvectionCoeff());
    }
    else {
      amrex::Abort("Invalid Convection Type.");
    }
  }
  else {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                              DragFunc,
                              NullConvectionCoeff());
  }
}

template <typename F1, typename F2>
void mfix::mfix_calc_transfer_coeffs (Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in,
                                      Vector< MultiFab* > const& X_gk_in,
                                      Vector< MultiFab* > const& pressure_g_in,
                                      F1 DragFunc,
                                      F2 ConvectionCoeff)
{
  if (m_reaction_rates_type == ReactionRatesType::RRatesUser) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                              DragFunc, ConvectionCoeff,
                              HeterogeneousRatesUser());
  } else {
    amrex::Abort("Invalid Reaction Rates Type.");
  }
}

template <typename F1, typename F2, typename F3>
void mfix::mfix_calc_transfer_coeffs (Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in,
                                      Vector< MultiFab* > const& X_gk_in,
                                      Vector< MultiFab* > const& pressure_g_in,
                                      F1 DragFunc,
                                      F2 ConvectionCoeff,
                                      F3 HeterogeneousRRates)
{
  using PairIndex = MFIXParticleContainer::PairIndex;
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  BL_PROFILE("mfix::mfix_calc_transfer_coeff()");

  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

  //***************************************************************************
  // Data for chemical reactions
  //***************************************************************************
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

  //***************************************************************************
  // 
  //***************************************************************************

  // We copy the value inside the domain to the outside to avoid
  // unphysical volume fractions.
  const int dir_bc_in = 2;
  mfix_set_epg_bcs(ep_g_in, dir_bc_in);

  for (int lev = 0; lev < nlev; lev++) {

    // This is just a sanity check to make sure we're not using covered values
    // We can remove these lines once we're confident in the algorithm
    EB_set_covered(*vel_g_in[lev], 0, 3, 1, covered_val);
    EB_set_covered(*ep_g_in[lev], 0, 1, 1, covered_val);
    EB_set_covered(*ro_g_in[lev], 0, 1, 1, covered_val);

    if (advect_enthalpy) {
      EB_set_covered(*T_g_in[lev], 0, 1, 1, covered_val);
    }

    if (solve_species) {
      EB_set_covered(*X_gk_in[lev], 0, fluid.nspecies, 1, covered_val);
    }

    if (reactions.solve) {
      EB_set_covered(*pressure_g_in[lev], 0, 1, 1, covered_val);
    }

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    MultiFab* interp_ptr;

    const int interp_ng = 1;    // Only one layer needed for interpolation
    const int interp_comp = 5 + // 3 vel_g + 1 ep_g + 1 ro_g
                            1*int(advect_enthalpy) +  // 1 T_g
                            fluid.nspecies*int(solve_species) +  // Ng X_gk
                            1*int(reactions.solve);  //  1 pressure_g

    MultiFab pressure_cc(ep_g_in[lev]->boxArray(), dmap[lev], 1, interp_ng);

    if (reactions.solve) {
      MultiFab pressure_nd(pressure_g_in[lev]->boxArray(), dmap[lev], 1, interp_ng);
      pressure_nd.setVal(0.);
      EB_set_covered(pressure_nd, 0, 1, 1, covered_val);

      MultiFab::Copy(pressure_nd, *m_leveldata[lev]->pressure_g, 0, 0, 1, interp_ng);
      MultiFab::Add (pressure_nd, *m_leveldata[lev]->p0_g, 0, 0, 1, interp_ng);

      pressure_cc.setVal(0.);
      EB_set_covered(pressure_cc, 0, 1, 1, covered_val);

      amrex::average_node_to_cellcenter(pressure_cc, 0, pressure_nd, 0, 1);
    }

    if (OnSameGrids) {

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(grids[lev], dmap[lev], interp_comp, interp_ng, MFInfo(), *ebfactory[lev]);

      int components_count(0);

      // Copy fluid velocity
      MultiFab::Copy(*interp_ptr, *vel_g_in[lev], 0, components_count, 3, interp_ng);
      components_count += 3;

      // Copy volume fraction
      MultiFab::Copy(*interp_ptr, *ep_g_in[lev],  0, components_count, 1, interp_ng);
      components_count += 1;

      // Copy fluid density
      MultiFab::Copy(*interp_ptr, *ro_g_in[lev],  0, components_count, 1, interp_ng);
      components_count += 1;

      if (advect_enthalpy) {
        // Copy fluid temperature
        MultiFab::Copy(*interp_ptr, *T_g_in[lev],  0, components_count, 1, interp_ng);
        components_count += 1;
      }

      if (solve_species) {
        // Copy volume fraction
        MultiFab::Copy(*interp_ptr, *X_gk_in[lev],  0, components_count, fluid.nspecies, interp_ng);
        components_count += fluid.nspecies;
      }

      if (reactions.solve) {
        // Copy volume fraction
        MultiFab::Copy(*interp_ptr, pressure_cc,  0, components_count, 1, interp_ng);
        components_count += 1;
      }
    }
    else
    {
      const BoxArray&            pba = pc->ParticleBoxArray(lev);
      const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

      // Temporary arrays  -- copies with no ghost cells
      //const int ng_to_copy = 0; UNUSED VARIABLE

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(pba, pdm, interp_comp, interp_ng, MFInfo(), *particle_ebfactory[lev]);

      int components_count(0);

      // Copy fluid velocity
      interp_ptr->ParallelCopy(*vel_g_in[lev], 0, components_count, 3, interp_ng, interp_ng);
      components_count += 3;

      // Copy volume fraction
      interp_ptr->ParallelCopy(*ep_g_in[lev],  0, components_count, 1, interp_ng, interp_ng);
      components_count += 1;

      // Copy fluid density
      interp_ptr->ParallelCopy(*ro_g_in[lev],  0, components_count, 1, interp_ng, interp_ng);
      components_count += 1;

      if (advect_enthalpy) {
        // Copy fluid temperature
        interp_ptr->ParallelCopy(*T_g_in[lev],  0, components_count, 1, interp_ng, interp_ng);
        components_count += 1;
      }

      if (solve_species) {
        // Copy fluid species
        interp_ptr->ParallelCopy(*X_gk_in[lev],  0, components_count, fluid.nspecies, interp_ng, interp_ng);
        components_count += fluid.nspecies;
      }

      if (reactions.solve) {
        // Copy fluid species
        interp_ptr->ParallelCopy(pressure_cc,  0, components_count, 1, interp_ng, interp_ng);
        components_count += 1;
      }
    }

    // FillBoundary on interpolation MultiFab
    interp_ptr->FillBoundary(geom[lev].periodicity());

    {
      const auto dxi_array = geom[lev].InvCellSizeArray();
      const auto dx_array  = geom[lev].CellSizeArray();
      const auto plo_array = geom[lev].ProbLoArray();

      const RealVect  dx( dx_array[0],  dx_array[1],  dx_array[2]);
      const RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
      const RealVect plo(plo_array[0], plo_array[1], plo_array[2]);

      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(interp_ptr->Factory());

      const auto cellcent = &(factory.getCentroid());
      const auto bndrycent = &(factory.getBndryCent());
      const auto areafrac = factory.getAreaFrac();
      const auto& volfrac = factory.getVolFrac();

      auto& fluid_parms = *fluid.parameters;
      auto& solids_parms = *solids.parameters;

      const int fluid_is_a_mixture = fluid.is_a_mixture;

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

        Box bx = pti.tilebox();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  interp_fab = static_cast<EBFArrayBox const&>((*interp_ptr)[pti]);
        const EBCellFlagFab&  flags = interp_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered) {

          const auto& interp_array = interp_ptr->const_array(pti);

          const auto& flags_array  = flags.array();

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

          const int adv_enthalpy = advect_enthalpy;
          const int fluid_is_a_mixture = fluid.is_a_mixture;
          const int solve_reactions = reactions.solve;

          const Real mu_g0 = fluid.mu_g0;

          const int nspecies_g = fluid.nspecies;

          auto local_cg_dem = DEM::cg_dem;

          amrex::ParallelFor(np,
              [pstruct,p_realarray,interp_array,DragFunc,ConvectionCoeff,
               HeterogeneousRRates,plo,dxi,adv_enthalpy,fluid_is_a_mixture,
               mu_g0,fluid_parms,nspecies_g,interp_comp,local_cg_dem,
               run_on_device,nspecies_s,ptile_data,nreactions,
               Solid,p_MW_sn,idx_X_sn,idx_mass_sn_txfr,idx_vel_s_txfr,
               idx_h_s_txfr,p_MW_gk,p_species_id_s,p_species_id_g,
               p_reactants_id,p_reactants_coeffs,p_reactants_phases,
               p_products_id,p_products_coeffs,p_products_phases,p_nreactants,
               p_nproducts,InvalidIdx,p_phases,p_nphases,p_types,Heterogeneous,
               T_ref,solids_parms,reactions_parms,flags_array,
               dem_solve,grown_bx_is_regular,dx,ccent_fab,bcent_fab,apx_fab,
               apy_fab,apz_fab,solve_reactions]
            AMREX_GPU_DEVICE (int p_id) noexcept
          {
            auto& particle = pstruct[p_id];

            // Cell containing particle centroid
            int ip = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
            int jp = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
            int kp = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

            const int cell_is_covered = static_cast<int>(flags_array(ip,jp,kp).isCovered());

            if (cell_is_covered) {

              p_realarray[SoArealData::dragcoeff][p_id] = 0.;

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

              GpuArray<Real,6+SPECIES::NMAX> interp_loc; // vel_g, ep_g, ro_g, T_g, X_gk
              interp_loc.fill(0.);

              GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

              if (grown_bx_is_regular) {

                trilinear_interp(particle.pos(), ip, jp, kp, weights,
                                 interp_loc.data(), interp_array, plo, dxi,
                                 interp_comp);

              } else {

                // Cut or regular cell and none of the cells in the stencil is
                // covered (Note we can't assume regular cell has no covered
                // cells in the stencil because of the diagonal case)

                // Upper cell in stencil
                int i = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0] + 0.5));
                int j = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1] + 0.5));
                int k = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2] + 0.5));

                // All cells in the stencil are regular. Use traditional
                // trilinear interpolation
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

#if 0
                  // TODO: This was initially split for variables that may have known
                  // EB values (e.g., no-slip velocity). However, the results changed
                  // more than expected so now EB values are not used.
                  {
                    const int srccomp = 0;
                    const int dstcomp = 0;
                    const int numcomp = 3;

                    shepard_interp_eb(particle.pos(), ip, jp, kp, dx, dxi, plo,
                                      flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                      interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);

                  }
                  {
                    const int srccomp = 3;
                    const int dstcomp = 3;
                    const int numcomp = interp_comp-3; // ep_g, ro_g, T_g, X_gk

                    shepard_interp(particle.pos(), ip, jp, kp, dx, dxi, plo,
                                   flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                   interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);
                  }
#else
                  const int srccomp = 0;
                  const int dstcomp = 0;
                  const int numcomp = interp_comp; // vel_g, ep_g, ro_g, T_g, X_gk

                  shepard_interp(particle.pos(), ip, jp, kp, dx, dxi, plo,
                                 flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                 interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);

#endif
                } // Cut cell
              }

              RealVect vel_g(interp_loc[0], interp_loc[1], interp_loc[2]);
              const Real ep_g = interp_loc[3];
              const Real ro_g = interp_loc[4];

              int comp_count = 5;

              Real T_g(0);
              if (adv_enthalpy) {
                T_g = interp_loc[comp_count];
                comp_count += 1;
              }

              Real* X_gk;
              if (fluid_is_a_mixture) {
                X_gk = &interp_loc[comp_count];
                comp_count += nspecies_g;
              }

              // Indices of cell where particle is located
              int iloc = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
              int jloc = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
              int kloc = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

              Real  mu_g(0);

              if (adv_enthalpy)
                mu_g = fluid_parms.calc_mu_g(T_g);
              else
                mu_g = mu_g0;

              Real rad = p_realarray[SoArealData::radius][p_id];
              Real vol = p_realarray[SoArealData::volume][p_id];

              int ip = particle.id();

              RealVect pvel(0.);
              pvel[0] = p_realarray[SoArealData::velx][p_id];
              pvel[1] = p_realarray[SoArealData::vely][p_id];
              pvel[2] = p_realarray[SoArealData::velz][p_id];

              Real rop_g = ro_g * ep_g;

              RealVect vslp(0.);
              vslp[0] = vel_g[0] - pvel[0];
              vslp[1] = vel_g[1] - pvel[1];
              vslp[2] = vel_g[2] - pvel[2];

              Real vrel = sqrt(dot_product(vslp, vslp));
              Real dp = 2.0*rad;

              if (local_cg_dem) {
                 dp = dp/std::cbrt(p_realarray[SoArealData::statwt][p_id]);
              }

              Real ep_s = 1.0 - ep_g;

              Real beta = vol*DragFunc(ep_g, mu_g, rop_g, vrel, dp, dp, ep_s,
                 vel_g[0], vel_g[1], vel_g[2], iloc, jloc, kloc, ip);

              p_realarray[SoArealData::dragcoeff][p_id] = beta;

              if(adv_enthalpy){
                Real k_g = fluid_parms.calc_k_g(T_g);

                Real cp_g(0);
                if (!fluid_is_a_mixture)
                  cp_g = run_on_device ?
                    fluid_parms.calc_cp_g<RunOn::Device>(T_g) :
                    fluid_parms.calc_cp_g<RunOn::Host>(T_g);
                else {
                  for (int n_g(0); n_g < nspecies_g; ++n_g) {
                    cp_g += run_on_device ?
                      X_gk[n_g]*fluid_parms.calc_cp_gk<RunOn::Device>(T_g,n_g) :
                      X_gk[n_g]*fluid_parms.calc_cp_gk<RunOn::Host>(T_g,n_g);
                  }
                }

                Real gamma = ConvectionCoeff(ep_g, mu_g, k_g, cp_g, rop_g, vrel,
                                             dp, iloc, jloc, kloc, ip);

                p_realarray[SoArealData::convection][p_id] = 4.0*M_PI*rad*rad*gamma;
              }

              if (solve_reactions) {
//                // Compute fluid molecular weight
//                Real MW_g(0);
//
//                if (fluid_is_a_mixture) {
//                  for (int n_g(0); n_g < nspecies_g; ++n_g) {
//                    MW_g += run_on_device ?
//                      X_gk[n_g] / fluid_parms.get_MW_gk<RunOn::Device>(n_g) :
//                      X_gk[n_g] / fluid_parms.get_MW_gk<RunOn::Host>(n_g);
//                  }
//
//                  MW_g = 1. / MW_g;
//                } else {
//                  MW_g = run_on_device ?
//                    fluid_parms.get_MW_g<RunOn::Device>() :
//                    fluid_parms.get_MW_g<RunOn::Host>();
//                }
//
//                const Real p_g = ro_g*fluid_parms.R*T_g / MW_g;
//
                // Extract interpolated thermodynamic pressure
                const Real p_g = interp_loc[comp_count];

                const Real ep_s = 1. - ep_g;
                const Real ro_p = p_realarray[SoArealData::density][p_id];

                GpuArray<Real,Reactions::NMAX> R_q_heterogeneous;
                R_q_heterogeneous.fill(0.);

                const Real T_p = p_realarray[SoArealData::temperature][p_id];
                const Real DP  = 2. * p_realarray[SoArealData::radius][p_id];
                const RealVect vel_p(p_realarray[SoArealData::velx][p_id],
                                     p_realarray[SoArealData::vely][p_id],
                                     p_realarray[SoArealData::velz][p_id]);

                if (run_on_device) {
                  HeterogeneousRRates.template operator()<RunOn::Device>(R_q_heterogeneous.data(),
                                                                         reactions_parms,
                                                                         solids_parms, X_sn.data(),
                                                                         ro_p, ep_s, T_p,
                                                                         vel_p, fluid_parms,
                                                                         X_gk, ro_g, ep_g, T_g,
                                                                         vel_g, DP, p_g);
                } else {
                  HeterogeneousRRates.template operator()<RunOn::Host>(R_q_heterogeneous.data(),
                                                                       reactions_parms,
                                                                       solids_parms, X_sn.data(),
                                                                       ro_p, ep_s, T_p,
                                                                       vel_p, fluid_parms,
                                                                       X_gk, ro_g, ep_g, T_g,
                                                                       vel_g, DP, p_g);
                }

                // Total transfer rates
                Real G_m_g_heterogeneous(0.);
                Real G_H_g_heterogeneous(0.);

                Real G_m_p_heterogeneous(0.);

                //*************************************************************
                // Initialize to zero
                //*************************************************************
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

                //*************************************************************
                // Loop over particle's species for computing particle txfr
                // rates
                //*************************************************************
                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  Real G_m_pk_heterogeneous(0.);

                  // Get the ID of the current species n_s
                  const int current_species_id = p_species_id_s[n_s];

                  // Loop over reactions to compute each contribution
                  for (int q(0); q < nreactions; q++) {
                    // Do something only if reaction is heterogeneous and
                    // contains a solid compound
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
                  }

                  G_m_p_heterogeneous += G_m_pk_heterogeneous;

                  // Update global variable
                  ptile_data.m_runtime_rdata[idx_mass_sn_txfr + n_s][p_id] = G_m_pk_heterogeneous;
                }

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

                      // Contribution to the particle
                      const Real h_gk_T_p = run_on_device ?
                        fluid_parms.calc_h_gk<RunOn::Device>(T_p, n_g, cell_is_covered) :
                        fluid_parms.calc_h_gk<RunOn::Host>(T_p, n_g, cell_is_covered);

                      const Real h_gk_T_g = run_on_device ?
                        fluid_parms.calc_h_gk<RunOn::Device>(T_g, n_g, cell_is_covered) :
                        fluid_parms.calc_h_gk<RunOn::Host>(T_g, n_g, cell_is_covered);

                      const Real G_H_pk_q = h_gk_T_g * amrex::min(0., G_m_gk_q);
                      const Real G_H_gk_q = h_gk_T_p * amrex::max(0., G_m_gk_q);

                      G_H_g_heterogeneous += (G_H_pk_q + G_H_gk_q);
                    }
                  }

                  G_m_g_heterogeneous += G_m_gk_heterogeneous;

                  // Update global variable
//                  G_ro_gk_ptr[n_g*np + ip] = G_m_gk_heterogeneous / fluid_vol;
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
//                G_h_g_ptr[p_id] = G_H_g_heterogeneous / fluid_vol;
              }
            }
          });
        } // if entire FAB not covered
      } // pti
    } // GPU region

    delete interp_ptr;
  } // lev


  // Reset the volume fractions back to the correct values at
  // inflow faces.
  const int dir_bc_out = 1;
  mfix_set_epg_bcs(ep_g_in, dir_bc_out);

}
