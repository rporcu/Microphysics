#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_K.H>
#include <mfix_filcc.H>

#include <mfix_mf_helpers.H>
#include <mfix_dem_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_algorithm.H>
#include <mfix_des_heterogeneous_rates_K.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>


void 
mfix::mfix_calc_chem_txfr (const Real time,
                           const Vector< MultiFab* >& ep_g_in,
                           const Vector< MultiFab* >& ro_g_in,
                           const Vector< MultiFab* >& X_gk_in)
{
  if (m_reaction_rates_type == ReactionRatesType::RRatesUser) {
    mfix_calc_chem_txfr(time, ep_g_in, ro_g_in, X_gk_in, ComputeRRateUser());
  }
  else {
    amrex::Abort("Invalid Reaction Rates Type.");
  }
}


template <typename F1>
void 
mfix::mfix_calc_chem_txfr (const Real time,
                           const Vector< MultiFab* >& ep_g_in,
                           const Vector< MultiFab* >& ro_g_in,
                           const Vector< MultiFab* >& X_gk_in,
                           F1 RRatesFunc)
{
  using PairIndex = MFIXParticleContainer::PairIndex;
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  const Real strttime = ParallelDescriptor::second();

  BL_PROFILE("mfix::mfix_calc_chem_txfr()");

  // Solid species data
  const int nspecies_s = SOLIDS::nspecies;

  Gpu::ManagedVector< int > mng_species_id_s(nspecies_s);
  Gpu::ManagedVector< Real > mng_MW_sn(nspecies_s);

  for (int n(0); n < nspecies_s; n++) {
    mng_species_id_s[n] = SOLIDS::species_id[n];
    mng_MW_sn[n] = SOLIDS::MW_sn0[n];
  }

  int* p_species_id_s = mng_species_id_s.data();
  Real* p_MW_sn = mng_MW_sn.data();

  // Fluid species data
  const int nspecies_g = FLUID::nspecies;

  Gpu::ManagedVector< int > mng_species_id_g(nspecies_g);
  Gpu::ManagedVector< Real > mng_MW_gk(nspecies_g);

  for (int n(0); n < nspecies_g; n++) {
    mng_species_id_g[n] = FLUID::species_id[n];
    mng_MW_gk[n] = FLUID::MW_gk0[n];
  }

  int* p_species_id_g = mng_species_id_g.data();
  Real* p_MW_gk = mng_MW_gk.data();

  // Reactions data
  const auto& chemical_reactions = REACTIONS::chemical_reactions;
  const int nreactions = REACTIONS::nreactions;

  Gpu::ManagedVector< int > mng_types(nreactions);
  Gpu::ManagedVector< const int* > mng_phases(nreactions);
  Gpu::ManagedVector< int > mng_nphases(nreactions);

  for (int q(0); q < nreactions; q++) {
    mng_types[q] = chemical_reactions[q].m_reaction_type;
    mng_phases[q] = chemical_reactions[q].m_phases.data();
    mng_nphases[q] = chemical_reactions[q].m_phases.size();
  }

  int* p_types = mng_types.data();
  const int** p_phases = mng_phases.data();
  int* p_nphases = mng_nphases.data();

  Gpu::ManagedVector< int > mng_nreactants(nreactions);
  Gpu::ManagedVector< const int* > mng_reactants_id(nreactions);
  Gpu::ManagedVector< const Real* > mng_reactants_coeffs(nreactions);
  Gpu::ManagedVector< const int* > mng_reactants_phases(nreactions);

  for (int q(0); q < nreactions; q++) {
    mng_nreactants[q] = chemical_reactions[q].m_reactants.size();
    mng_reactants_id[q] = chemical_reactions[q].m_reactants_id.data();
    mng_reactants_coeffs[q] = chemical_reactions[q].m_reactants_coeffs.data();
    mng_reactants_phases[q] = chemical_reactions[q].m_reactants_phases.data();
  }

  int* p_nreactants = mng_nreactants.data();
  const int** p_reactants_id = mng_reactants_id.data();
  const Real** p_reactants_coeffs = mng_reactants_coeffs.data();
  const int** p_reactants_phases = mng_reactants_phases.data();

  Gpu::ManagedVector< int > mng_nproducts(nreactions);
  Gpu::ManagedVector< const int* > mng_products_id(nreactions);
  Gpu::ManagedVector< const Real* > mng_products_coeffs(nreactions);
  Gpu::ManagedVector< const int* > mng_products_phases(nreactions);

  for (int q(0); q < nreactions; q++) {
    mng_nproducts[q] = chemical_reactions[q].m_products.size();
    mng_products_id[q] = chemical_reactions[q].m_products_id.data();
    mng_products_coeffs[q] = chemical_reactions[q].m_products_coeffs.data();
    mng_products_phases[q] = chemical_reactions[q].m_products_phases.data();
  }

  int* p_nproducts = mng_nproducts.data();
  const int** p_products_id = mng_products_id.data();
  const Real** p_products_coeffs = mng_products_coeffs.data();
  const int** p_products_phases = mng_products_phases.data();

  // Solid phase integer ID
  const int Solid = CHEMICALPHASE::Solid;
  const int Heterogeneous = REACTIONTYPE::Heterogeneous;
  const int InvalidIdx = -1; //TODO define this somewhere else

  // Particles indexes
  const int idx_X = speciesData::X_sn*nspecies_s;

  const int idx_G = speciesData::count*nspecies_s + reactionsData::G_sn_pg_q*nreactions;

  // START MFIX_CALC_CHEM_TRANSFER_COEFFS
  {
    // We copy the value inside the domain to the outside to avoid
    // unphysical volume fractions.
    const int dir_bc_in = 2;
    mfix_set_epg_bcs(ep_g_in, dir_bc_in);

    // Set boundary conditions just in case
    mfix_set_density_bcs(time, ro_g_in);

    // TODO: we only need to do this on X_gk
    mfix_set_species_bcs(time, X_gk_in, get_D_gk(), get_cp_gk(), get_h_gk());

    for (int lev = 0; lev < nlev; lev++)
    {
      bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) and
                           (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

      MultiFab* interp_ptr;

      EB_set_covered(*ep_g_in[0], 0, 1, 1, covered_val);
      EB_set_covered(*ro_g_in[0], 0, 1, 1, covered_val);
      EB_set_covered(*X_gk_in[0], 0, FLUID::nspecies, 1, covered_val);

      const int interp_ng = 1;    // Only one layer needed for interpolation

      // Fluid species components + 1 ep_g + 1 ro_g
      const int interp_comp = FLUID::nspecies+2;  

      if (OnSameGrids)
      {
        // Store gas velocity and volume fraction for interpolation
        interp_ptr = new MultiFab(grids[lev], dmap[lev], interp_comp, interp_ng,
            MFInfo(), *ebfactory[lev]);

        // Copy fluid species mass fractions
        interp_ptr->copy(*X_gk_in[lev], 0, 0, FLUID::nspecies, interp_ng, interp_ng);

        // Copy volume fraction
        interp_ptr->copy(*ep_g_in[lev], 0, interp_comp-2, 1, interp_ng, interp_ng);

        // Copy density
        interp_ptr->copy(*ro_g_in[lev], 0, interp_comp-1, 1, interp_ng, interp_ng);

        interp_ptr->FillBoundary(geom[lev].periodicity());
      }
      else
      {
        const BoxArray&            pba = pc->ParticleBoxArray(lev);
        const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

        EBFArrayBoxFactory ebfactory_loc(*eb_levels[lev], geom[lev], pba, pdm,
                                         {m_eb_basic_grow_cells,
                                          m_eb_volume_grow_cells,
                                          m_eb_full_grow_cells},
                                         EBSupport::full);

        // Store gas velocity and volume fraction for interpolation
        interp_ptr = new MultiFab(pba, pdm, interp_comp, interp_ng, MFInfo(),
            ebfactory_loc);

        // Copy fluid velocity
        interp_ptr->copy(*X_gk_in[lev], 0, 0, FLUID::nspecies, interp_ng, interp_ng);

        // Copy volume fraction
        interp_ptr->copy(*ep_g_in[lev], 0, interp_comp-2, 1, interp_ng, interp_ng);

        // Copy density
        interp_ptr->copy(*ro_g_in[lev], 0, interp_comp-1, 1, interp_ng, interp_ng);

        interp_ptr->FillBoundary(geom[lev].periodicity());
      }

      {
        const auto dxi_array = geom[lev].InvCellSizeArray();
        const auto dx_array  = geom[lev].CellSizeArray();
        const auto plo_array = geom[lev].ProbLoArray();

        const amrex::RealVect  dx( dx_array[0],  dx_array[1],  dx_array[2]);
        const amrex::RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
        const amrex::RealVect plo(plo_array[0], plo_array[1], plo_array[2]);

        const auto& factory =
          dynamic_cast<EBFArrayBoxFactory const&>(interp_ptr->Factory());

        const auto cellcent = &(factory.getCentroid());
        const auto bndrycent = &(factory.getBndryCent());
        const auto areafrac = factory.getAreaFrac();

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

          const int np = particles.size();

          Box bx = pti.tilebox();

          // This is to check efficiently if this tile contains any eb stuff
          const EBFArrayBox& interp_fab =
            static_cast<EBFArrayBox const&>((*interp_ptr)[pti]);

          const EBCellFlagFab& flags = interp_fab.getEBCellFlagFab();

          if (flags.getType(amrex::grow(bx,0)) != FabType::covered)
          {
            const auto& interp_array = interp_ptr->array(pti);

            const auto& flags_array = flags.array();

            if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
            {
              amrex::ParallelFor(np,
                [pstruct,interp_array,RRatesFunc,plo,dxi,ptile_data,
                 nspecies_g,nspecies_s,nreactions,interp_comp,Solid,
                 p_MW_sn,idx_X,idx_G,p_MW_gk,p_species_id_s,p_species_id_g,
                 p_reactants_id,p_reactants_coeffs,p_reactants_phases,
                 p_products_id,p_products_coeffs,p_products_phases,p_nreactants,
                 p_nproducts,InvalidIdx,p_phases,p_nphases,p_types,Heterogeneous]
              AMREX_GPU_DEVICE (int p_id) noexcept
              {
                auto& particle = pstruct[p_id];
                
                Real* X_sn = new Real [nspecies_s];

                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  const int idx = idx_X + n_s;
                  X_sn[n_s] = ptile_data.m_runtime_rdata[idx][p_id];
                }

                // Pointer to this particle's species rate of formation
                Real** G_sn_pg_q = new Real* [nspecies_s];

                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  G_sn_pg_q[n_s] = new Real [nreactions];
                  for (int q(0); q < nreactions; q++)
                    G_sn_pg_q[n_s][q] = 0;
                }

                Real* interp_loc = new Real[interp_comp];

                for(int n(0); n < interp_comp; n++)
                  interp_loc[n] = 0;

                Real* R_q = new Real[nreactions];

                for(int q(0); q < nreactions; q++)
                  R_q[q] = 0;

                trilinear_interp(particle.pos(), interp_loc,
                                 interp_array, plo, dxi, interp_comp);

                Real* X_gk = new Real [nspecies_g];

                for (int n_g(0); n_g < nspecies_g; n_g++)
                  X_gk[n_g] = interp_loc[n_g];
                
                Real ep_g = interp_loc[interp_comp-2];

                Real ro_g = interp_loc[interp_comp-1];

                Real ep_s = 1. - ep_g;

                Real ro_s = particle.rdata(realData::mass) /
                  particle.rdata(realData::volume);

                RRatesFunc(R_q, nreactions, p_nreactants, p_nproducts,
                    p_reactants_id, p_reactants_coeffs, p_reactants_phases,
                    p_products_id, p_products_coeffs, p_products_phases,
                    p_species_id_s, X_sn, p_MW_sn, nspecies_s, ro_s, ep_s,
                    p_species_id_g, X_gk, p_MW_gk, nspecies_g, ro_g, ep_g);

                for (int n_s(0); n_s < nspecies_s; n_s++)
                {
                  const int current_species_id = p_species_id_s[n_s];

                  for (int q(0); q < nreactions; q++)
                  {
                    // Do something only if reaction is heterogeneous and contains
                    // a solid compound
                    if (p_types[q] == Heterogeneous and
                        MFIXfind(p_phases[q], p_nphases[q], Solid) != InvalidIdx)
                    {
                      Real stoc_coeff(0);

                      // Add reactant contribution (if any)
                      {
                        const int pos = MFIXfind(p_reactants_id[q],
                            p_nreactants[q], current_species_id);

                        if (pos != InvalidIdx) {
                          if (p_reactants_phases[q][pos] == Solid)
                            stoc_coeff += p_reactants_coeffs[q][pos];
                        }
                      }

                      // Add products contribution (if any)
                      {
                        const int pos = MFIXfind(p_products_id[q], p_nproducts[q],
                            current_species_id);

                        if (pos != InvalidIdx) {
                          if (p_products_phases[q][pos] == Solid)
                            stoc_coeff += p_products_coeffs[q][pos];
                        }
                      }

                      G_sn_pg_q[n_s][q] = stoc_coeff * p_MW_sn[n_s] * R_q[q];
                    }
                    else {
                      G_sn_pg_q[n_s][q] = 0;
                    }
                  }
                }

                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  for (int q(0); q < nreactions; q++) {
                    const int idx = idx_G + n_s*nreactions + q;
                    ptile_data.m_runtime_rdata[idx][p_id] = G_sn_pg_q[n_s][q];
                  }
                }

                delete[] interp_loc;
                delete[] X_gk;
                delete[] X_sn;
                delete[] R_q;

                for (int n_s(0); n_s < nspecies_s; n_s++)
                  delete[] G_sn_pg_q[n_s];

                delete[] G_sn_pg_q;
              });
            }
            else // FAB not all regular
            {
              // Cell centroids
              const auto& ccent_fab = cellcent->array(pti);
              // Centroid of EB
              const auto& bcent_fab = bndrycent->array(pti);
              // Area fractions
              const auto& apx_fab = areafrac[0]->array(pti);
              const auto& apy_fab = areafrac[1]->array(pti);
              const auto& apz_fab = areafrac[2]->array(pti);

              amrex::ParallelFor(np,
                [pstruct,interp_array,RRatesFunc,plo,dx,dxi,flags_array,
                 ccent_fab,bcent_fab,apx_fab,apy_fab, apz_fab,nspecies_s,
                 ptile_data,interp_comp,nspecies_g,nreactions,p_species_id_s,
                 p_species_id_g,p_types,p_phases,p_nphases,p_products_id,
                 p_products_coeffs,p_products_phases,p_reactants_id,
                 p_reactants_coeffs,p_reactants_phases,InvalidIdx,p_MW_sn,idx_X,
                 idx_G,p_MW_gk,Solid,p_nreactants,p_nproducts,Heterogeneous]
              AMREX_GPU_DEVICE (int p_id) noexcept
              {
                auto& particle = pstruct[p_id];

                // Pointer to this particle's species rate of formation
                Real** G_sn_pg_q = new Real* [nspecies_s];
                
                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  G_sn_pg_q[n_s] = new Real [nreactions];
                  for (int q(0); q < nreactions; q++)
                    G_sn_pg_q[n_s][q] = 0.;
                }
                  
                // Cell containing particle centroid
                int ip = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
                int jp = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
                int kp = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

                // No drag force for particles in covered cells.
                if (flags_array(ip,jp,kp).isCovered())
                {
                  // Particle mass calculation
                  for (int n_s(0); n_s < nspecies_s; n_s++) {
                    for (int q(0); q < nreactions; q++) {
                      // Update specie mass formation rate
                      G_sn_pg_q[n_s][q] = 0;
                    }
                  }
                }
                else {
                  // Cut or regular cell and none of the cells in the stencil is
                  // covered (Note we can't assume regular cell has no covered
                  // cells in the stencil because of the diagonal case)

                  // Upper cell in trilinear stencil
                  int i = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0] + 0.5));
                  int j = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1] + 0.5));
                  int k = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2] + 0.5));

                  // Local array storing interpolated values
                  Real* interp_loc = new Real [interp_comp];

                  // All cells in the stencil are regular. Use
                  // traditional trilinear interpolation
                  if (flags_array(i-1,j-1,k-1).isRegular() and
                      flags_array(i  ,j-1,k-1).isRegular() and
                      flags_array(i-1,j  ,k-1).isRegular() and
                      flags_array(i  ,j  ,k-1).isRegular() and
                      flags_array(i-1,j-1,k  ).isRegular() and
                      flags_array(i  ,j-1,k  ).isRegular() and
                      flags_array(i-1,j  ,k  ).isRegular() and
                      flags_array(i  ,j  ,k  ).isRegular()) {

                    trilinear_interp(particle.pos(), interp_loc,
                                     interp_array, plo, dxi, interp_comp);
                  // At least one of the cells in the stencil is cut or covered
                  }
                  else
                  {
                    // All the quantities to interpolate are scalar quantities
                    // No velocity to be interpolated
                    const int scomp = 0;

                    fe_interp(particle.pos(), ip, jp, kp, dx, dxi, plo, flags_array,
                              ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                              interp_array, interp_loc, interp_comp, scomp);
                  } // Cut cell

                  Real* X_gk = new Real [nspecies_g];
                  Real* X_sn = new Real [nspecies_s];
                  Real* R_q = new Real[nreactions];

                  BL_ASSERT(interp_comp-2 == nspecies_g);

                  for (int n_s(0); n_s < nspecies_s; n_s++) {
                    const int idx = idx_X + n_s;
                    X_sn[n_s] = ptile_data.m_runtime_rdata[idx][p_id];
                  }
                  
                  for (int n_g(0); n_g < nspecies_g; n_g++)
                    X_gk[n_g] = interp_loc[n_g];
                  
                  Real ep_g = interp_loc[interp_comp-2];
                  Real ro_g = interp_loc[interp_comp-1];

                  Real ep_s = 1. - ep_g;

                  Real ro_s = particle.rdata(realData::mass) /
                    particle.rdata(realData::volume);

                  RRatesFunc(R_q, nreactions, p_nreactants, p_nproducts,
                      p_reactants_id, p_reactants_coeffs, p_reactants_phases,
                      p_products_id, p_products_coeffs, p_products_phases,
                      p_species_id_s, X_sn, p_MW_sn, nspecies_s, ro_s, ep_s,
                      p_species_id_g, X_gk, p_MW_gk, nspecies_g, ro_g, ep_g);

                  for (int n_s(0); n_s < nspecies_s; n_s++)
                  {
                    const int current_species_id = p_species_id_s[n_s];

                    for (int q(0); q < nreactions; q++)
                    {
                      // Do something only if reaction is heterogeneous and contains
                      // a solid compound
                      if (p_types[q] == Heterogeneous and
                          MFIXfind(p_phases[q], p_nphases[q], Solid) != InvalidIdx)
                      {
                        Real stoc_coeff(0);

                        // Add reactant contribution (if any)
                        {
                          const int pos = MFIXfind(p_reactants_id[q], p_nreactants[q],
                              current_species_id);

                          if (pos != InvalidIdx) {
                            if (p_reactants_phases[q][pos] == Solid)
                              stoc_coeff += p_reactants_coeffs[q][pos];
                          }
                        }

                        // Add products contribution (if any)
                        {
                          const int pos = MFIXfind(p_products_id[q], p_nproducts[q],
                              current_species_id);

                          if (pos != InvalidIdx) {
                            if (p_products_phases[q][pos] == Solid)
                              stoc_coeff += p_products_coeffs[q][pos];
                          }
                        }

                        G_sn_pg_q[n_s][q] = stoc_coeff * p_MW_sn[n_s] * R_q[q];
                      }
                      else {
                        G_sn_pg_q[n_s][q] = 0;
                      }
                    }
                  }

                  for (int n_s(0); n_s < nspecies_s; n_s++) {
                    for (int q(0); q < nreactions; q++) {
                      const int idx = idx_G + n_s*nreactions + q;

                      ptile_data.m_runtime_rdata[idx][p_id] = G_sn_pg_q[n_s][q];
                    }
                  }

                  delete[] interp_loc;
                  delete[] X_gk;
                  delete[] X_sn;
                  delete[] R_q;

                } // Not covered

                for (int n_s(0); n_s < nspecies_s; n_s++)
                  delete[] G_sn_pg_q[n_s];

                delete[] G_sn_pg_q;

              }); // p_id

            } // type of FAB
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
  // END MFIX_CALC_CHEM_TRANSFER_COEFFS

  {
    // ***************************************************************************
    // Now use the chem transfer coeffs of individual particles to create the
    // interphase transfer terms on the fluid
    // ***************************************************************************
    for (int lev = 0; lev < nlev; lev++)
      m_leveldata[lev]->ro_gk_txfr->setVal(0);

    if (nlev > 2)
      amrex::Abort("For right now"
          " MFIXParticleContainer::TrilinearDepositionFluidRRates can only"
          " handle up to 2 levels");

    Vector< MultiFab* > ro_gk_txfr_ptr(nlev, nullptr);

    for (int lev = 0; lev < nlev; lev++)
    {
      bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) and
                           (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

      if (lev == 0 and OnSameGrids)
      {
        // If we are already working with the internal mf defined on the
        // particle_box_array, then we just work with this.
        ro_gk_txfr_ptr[lev] = m_leveldata[lev]->ro_gk_txfr;
      }
      else if (lev == 0 and (not OnSameGrids))
      {
        // If beta_mf is not defined on the particle_box_array, then we need
        // to make a temporary here and copy into beta_mf at the end.
        ro_gk_txfr_ptr[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                           pc->ParticleDistributionMap(lev),
                                           m_leveldata[lev]->ro_gk_txfr->nComp(),
                                           m_leveldata[lev]->ro_gk_txfr->nGrow());
      }
      else
      {
        // If lev > 0 we make a temporary at the coarse resolution
        BoxArray ba_crse(amrex::coarsen(pc->ParticleBoxArray(lev),
              this->m_gdb->refRatio(0)));

        ro_gk_txfr_ptr[lev] = new MultiFab(ba_crse, pc->ParticleDistributionMap(lev),
                                           m_leveldata[lev]->ro_gk_txfr->nComp(), 1);
      }

      // We must have ghost cells for each FAB so that a particle in one grid can
      // spread its effect to an adjacent grid by first putting the value into
      // ghost cells of its own grid.  The mf->sumBoundary call then adds the
      // value from one grid's ghost cell to another grid's valid region.
      if (ro_gk_txfr_ptr[lev]->nGrow() < 1)
        amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

      ro_gk_txfr_ptr[lev]->setVal(0.0, 0, m_leveldata[lev]->ro_gk_txfr->nComp(),
          ro_gk_txfr_ptr[lev]->nGrow());
    }

    const Geometry& gm = Geom(0);
    const FabArray<EBCellFlagFab>* flags = nullptr;
    const MultiFab* volfrac = nullptr;

    Vector< MultiFab* > tmp_eps(nlev);

    for (int lev = 0; lev < nlev; lev++)
    {
      tmp_eps[lev] = (MFHelpers::createFrom(*ro_gk_txfr_ptr[lev], 0.0)).release();

      // Use level 0 to define the EB factory. If we are not on level 0
      // then create a copy of the coarse factory to use.

      if (lev == 0)
      {
        flags   = &(particle_ebfactory[lev]->getMultiEBCellFlagFab());
        volfrac = &(particle_ebfactory[lev]->getVolFrac());
      }
      else
      {
        Vector<int> ngrow = {1,1,1};
        EBFArrayBoxFactory* crse_factory;

        crse_factory = (makeEBFabFactory(gm, ro_gk_txfr_ptr[lev]->boxArray(),
                                         ro_gk_txfr_ptr[lev]->DistributionMap(),
                                         ngrow, EBSupport::volume)).release();

        flags   = &(crse_factory->getMultiEBCellFlagFab());
        volfrac = &(crse_factory->getVolFrac());

        delete crse_factory;
      }

      // Deposit the interphase transfer forces to the grid
      // Drag force: (beta and beta*particle_vel)
      // Heat transfer: gamma and gamma*particle temperature
      pc->InterphaseChemDeposition(lev, *tmp_eps[lev], *ro_gk_txfr_ptr[lev],
          volfrac, flags);
    }

    {
      // The deposition occurred on level 0, thus the next few operations
      // only need to be carried out on level 0.
      int lev(0);

      // Move any volume deposited outside the domain back into the domain
      // when BC is either a pressure inlet or mass inflow.
      mfix_deposition_bcs(lev, *ro_gk_txfr_ptr[lev]);

      // Sum grid boundaries to capture any material that was deposited into
      // your grid from an adjacent grid.
      ro_gk_txfr_ptr[lev]->SumBoundary(gm.periodicity());
      ro_gk_txfr_ptr[lev]->setBndry(0.0);

      // Sum grid boundaries then fill with correct ghost values.
      tmp_eps[lev]->SumBoundary(gm.periodicity());
      tmp_eps[lev]->FillBoundary(gm.periodicity());

      // Move excessive solids volume from small cells to neighboring cells.
      // Note that we don't change tmp_eps but use the redistribution of
      // particle volume to determine how to redistribute the drag forces.
      mfix_redistribute_deposition(lev, *tmp_eps[lev], *ro_gk_txfr_ptr[lev],
                                   volfrac, flags,
                                   mfix::m_max_solids_volume_fraction);

      // Sum the boundaries again to recapture any solids moved across
      // grid boundaries during the redistribute
      ro_gk_txfr_ptr[lev]->SumBoundary(gm.periodicity());
      ro_gk_txfr_ptr[lev]->FillBoundary(gm.periodicity());
    }

    // This might not need to exist on all levels. Maybe only level 0.
    for (int lev(0); lev < nlev; ++lev)
      delete tmp_eps[lev];

    int  src_nghost = 1;
    int dest_nghost = 0;
    int ng_to_copy = amrex::min(src_nghost, dest_nghost);

    for (int lev = 1; lev < nlev; lev++) {
      ro_gk_txfr_ptr[0]->copy(*ro_gk_txfr_ptr[lev], 0, 0, ro_gk_txfr_ptr[0]->nComp(),
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

        m_leveldata[lev]->ro_gk_txfr->setVal(0);

        amrex::InterpFromCoarseLevel(*m_leveldata[lev]->ro_gk_txfr, time,
                                     *ro_gk_txfr_ptr[lev-1],
                                     0, 0, 1, Geom(lev-1), Geom(lev),
                                     cphysbc, 0, fphysbc, 0,
                                     ref_ratio[0], mapper,
                                     bcs, 0);
      }
    }

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from txfr_ptr into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.

    if (ro_gk_txfr_ptr[0] != m_leveldata[0]->ro_gk_txfr) {
      m_leveldata[0]->ro_gk_txfr->copy(*ro_gk_txfr_ptr[0], 0, 0,
          m_leveldata[0]->ro_gk_txfr->nComp());
    }

    for (int lev = 0; lev < nlev; lev++) {
      if (ro_gk_txfr_ptr[lev] != m_leveldata[lev]->ro_gk_txfr)
        delete ro_gk_txfr_ptr[lev];
    }

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,
                                        ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::TrilinearDepositionFluidRrates"
        " time: " << stoptime << '\n';
    }

    // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
    for (int lev = 0; lev < nlev; lev++)
      m_leveldata[lev]->ro_gk_txfr->FillBoundary(geom[lev].periodicity());
  }
}
