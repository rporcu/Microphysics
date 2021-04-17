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
mfix::mfix_calc_chem_txfr (const Vector< MultiFab* >& chem_txfr,
                           const Vector< MultiFab* >& ep_g_in,
                           const Vector< MultiFab* >& ro_g_in,
                           const Vector< MultiFab* >& X_gk_in,
                           const Vector< MultiFab* >& D_gk_in,
                           const Vector< MultiFab* >& cp_gk_in,
                           const Vector< MultiFab* >& h_gk_in,
                           const Real time)
{
  if (m_reaction_rates_type == ReactionRatesType::RRatesUser) {
    mfix_calc_chem_txfr(chem_txfr, ep_g_in, ro_g_in, X_gk_in, D_gk_in, cp_gk_in,
                        h_gk_in, time, ComputeRRateUser());
  }
  else {
    amrex::Abort("Invalid Reaction Rates Type.");
  }
}


template <typename F1>
void 
mfix::mfix_calc_chem_txfr (const Vector< MultiFab* >& chem_txfr,
                           const Vector< MultiFab* >& ep_g_in,
                           const Vector< MultiFab* >& ro_g_in,
                           const Vector< MultiFab* >& X_gk_in,
                           const Vector< MultiFab* >& D_gk_in,
                           const Vector< MultiFab* >& cp_gk_in,
                           const Vector< MultiFab* >& h_gk_in,
                           const Real time,
                           F1 RRatesFunc)
{
  using PairIndex = MFIXParticleContainer::PairIndex;
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  const Real strttime = ParallelDescriptor::second();

  BL_PROFILE("mfix::mfix_calc_chem_txfr()");

  // Solid species data
  const int nspecies_s = SOLIDS::nspecies;

  Gpu::DeviceVector< int > d_species_id_s(nspecies_s);
  Gpu::DeviceVector< Real > d_MW_sn(nspecies_s);
  Gpu::copyAsync(Gpu::hostToDevice, SOLIDS::species_id.begin(), SOLIDS::species_id.end(), d_species_id_s.begin());
  Gpu::copyAsync(Gpu::hostToDevice, SOLIDS::MW_sn0.begin(), SOLIDS::MW_sn0.end(), d_MW_sn.begin());
  int* p_species_id_s = d_species_id_s.data();
  Real* p_MW_sn = d_MW_sn.data();

  // Fluid species data
  const int nspecies_g = FLUID::nspecies;

  Gpu::DeviceVector< int > d_species_id_g(nspecies_g);
  Gpu::DeviceVector< Real > d_MW_gk(nspecies_g);
  Gpu::copyAsync(Gpu::hostToDevice, FLUID::species_id.begin(), FLUID::species_id.end(), d_species_id_g.begin());
  Gpu::copyAsync(Gpu::hostToDevice, FLUID::MW_gk0.begin(), FLUID::MW_gk0.end(), d_MW_gk.begin());
  int* p_species_id_g = d_species_id_g.data();
  Real* p_MW_gk = d_MW_gk.data();

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
  const int Heterogeneous = REACTIONTYPE::Heterogeneous;
  const int InvalidIdx = -1; //TODO define this somewhere else

  // Particles indexes
  const int idx_X = SoAspeciesData::X_sn*nspecies_s;

  const int idx_G = SoAspeciesData::count*nspecies_s + SoAreactionsData::G_sn_pg_q*nreactions;

  Gpu::synchronize();

  // START MFIX_CALC_CHEM_TRANSFER_COEFFS
  {
    // We copy the value inside the domain to the outside to avoid
    // unphysical volume fractions.
    const int dir_bc_in = 2;
    mfix_set_epg_bcs(ep_g_in, dir_bc_in);

    // Set boundary conditions just in case
    mfix_set_density_bcs(time, ro_g_in);

    // TODO: we only need to do this on X_gk
    mfix_set_species_bcs(time, X_gk_in, D_gk_in, cp_gk_in, h_gk_in);

    for (int lev = 0; lev < nlev; lev++)
    {
      bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
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
        MultiFab::Copy(*interp_ptr, *X_gk_in[lev], 0, 0, FLUID::nspecies, interp_ng);

        // Copy volume fraction
        MultiFab::Copy(*interp_ptr, *ep_g_in[lev], 0, interp_comp-2, 1, interp_ng);

        // Copy density
        MultiFab::Copy(*interp_ptr, *ro_g_in[lev], 0, interp_comp-1, 1, interp_ng);

        interp_ptr->FillBoundary(geom[lev].periodicity());
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

        const RealVect  dx( dx_array[0],  dx_array[1],  dx_array[2]);
        const RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
        const RealVect plo(plo_array[0], plo_array[1], plo_array[2]);

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

          auto& soa = pti.GetStructOfArrays();
          auto p_realarray = soa.realarray();

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
                [pstruct,p_realarray,interp_array,RRatesFunc,plo,dxi,ptile_data,
                 nspecies_g,nspecies_s,nreactions,interp_comp,Solid,
                 p_MW_sn,idx_X,idx_G,p_MW_gk,p_species_id_s,p_species_id_g,
                 p_reactants_id,p_reactants_coeffs,p_reactants_phases,
                 p_products_id,p_products_coeffs,p_products_phases,p_nreactants,
                 p_nproducts,InvalidIdx,p_phases,p_nphases,p_types,Heterogeneous]
              AMREX_GPU_DEVICE (int p_id) noexcept
              {
                auto& particle = pstruct[p_id];
                
                GpuArray<Real,SPECIES::NMAX> X_sn;

                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  const int idx = idx_X + n_s;
                  X_sn[n_s] = ptile_data.m_runtime_rdata[idx][p_id];
                }

                // Pointer to this particle's species rate of formation
                GpuArray<GpuArray<Real,REACTIONS::NMAX>,SPECIES::NMAX> G_sn_pg_q;

                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  for (int q(0); q < nreactions; q++)
                    G_sn_pg_q[n_s][q] = 0;
                }

                GpuArray<Real,SPECIES::NMAX+2> interp_loc;

                for(int n(0); n < interp_comp; n++)
                  interp_loc[n] = 0;

                GpuArray<Real,REACTIONS::NMAX> R_q;

                for(int q(0); q < nreactions; q++)
                  R_q[q] = 0;

                trilinear_interp(particle.pos(), interp_loc.data(),
                                 interp_array, plo, dxi, interp_comp);

                GpuArray<Real,SPECIES::NMAX> X_gk;

                for (int n_g(0); n_g < nspecies_g; n_g++)
                  X_gk[n_g] = interp_loc[n_g];
                
                Real ep_g = interp_loc[interp_comp-2];

                Real ro_g = interp_loc[interp_comp-1];

                Real ep_s = 1. - ep_g;

                Real ro_s = p_realarray[SoArealData::mass][p_id] /
                  p_realarray[SoArealData::volume][p_id];

                RRatesFunc(R_q.data(), nreactions, p_nreactants, p_nproducts,
                    p_reactants_id, p_reactants_coeffs, p_reactants_phases,
                    p_products_id, p_products_coeffs, p_products_phases,
                    p_species_id_s, X_sn.data(), p_MW_sn, nspecies_s, ro_s, ep_s,
                    p_species_id_g, X_gk.data(), p_MW_gk, nspecies_g, ro_g, ep_g);

                for (int n_s(0); n_s < nspecies_s; n_s++)
                {
                  const int current_species_id = p_species_id_s[n_s];

                  for (int q(0); q < nreactions; q++)
                  {
                    // Do something only if reaction is heterogeneous and contains
                    // a solid compound
                    if (p_types[q] == Heterogeneous &&
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
                [pstruct,p_realarray,interp_array,RRatesFunc,plo,dx,dxi,flags_array,
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
                GpuArray<GpuArray<Real,REACTIONS::NMAX>,SPECIES::NMAX> G_sn_pg_q;
                
                for (int n_s(0); n_s < nspecies_s; n_s++) {
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
                  GpuArray<Real,SPECIES::NMAX+2> interp_loc;

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

                    trilinear_interp(particle.pos(), interp_loc.data(),
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
                              interp_array, interp_loc.data(), interp_comp, scomp);
                  } // Cut cell

                  GpuArray<Real,SPECIES::NMAX> X_gk;
                  GpuArray<Real,SPECIES::NMAX> X_sn;
                  GpuArray<Real,REACTIONS::NMAX> R_q;

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

                  Real ro_s = p_realarray[SoArealData::mass][p_id] /
                    p_realarray[SoArealData::volume][p_id];

                  RRatesFunc(R_q.data(), nreactions, p_nreactants, p_nproducts,
                      p_reactants_id, p_reactants_coeffs, p_reactants_phases,
                      p_products_id, p_products_coeffs, p_products_phases,
                      p_species_id_s, X_sn.data(), p_MW_sn, nspecies_s, ro_s, ep_s,
                      p_species_id_g, X_gk.data(), p_MW_gk, nspecies_g, ro_g, ep_g);

                  for (int n_s(0); n_s < nspecies_s; n_s++)
                  {
                    const int current_species_id = p_species_id_s[n_s];

                    for (int q(0); q < nreactions; q++)
                    {
                      // Do something only if reaction is heterogeneous and contains
                      // a solid compound
                      if (p_types[q] == Heterogeneous &&
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
                } // Not covered
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
      chem_txfr[lev]->setVal(0);

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
        chem_txfr_ptr[lev] = chem_txfr[lev];
      }
      else if (lev == 0 && (!OnSameGrids))
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
    const FabArray<EBCellFlagFab>* flags = nullptr;
    const MultiFab* volfrac = nullptr;

    Vector< MultiFab* > tmp_eps(nlev);

    for (int lev = 0; lev < nlev; lev++)
    {
      tmp_eps[lev] = (MFHelpers::createFrom(*chem_txfr_ptr[lev], 0.0)).release();

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

        crse_factory = (makeEBFabFactory(gm, chem_txfr_ptr[lev]->boxArray(),
                                         chem_txfr_ptr[lev]->DistributionMap(),
                                         ngrow, EBSupport::volume)).release();

        flags   = &(crse_factory->getMultiEBCellFlagFab());
        volfrac = &(crse_factory->getVolFrac());

        delete crse_factory;
      }

      // Deposit the interphase transfer forces to the grid
      // Drag force: (beta and beta*particle_vel)
      // Heat transfer: gamma and gamma*particle temperature
      pc->InterphaseChemDeposition(lev, *tmp_eps[lev], *chem_txfr_ptr[lev],
          volfrac, flags, m_chemical_reactions);
    }

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
      mfix_redistribute_deposition(lev, *tmp_eps[lev], *chem_txfr_ptr[lev],
                                   volfrac, flags,
                                   mfix::m_max_solids_volume_fraction);

      // Sum the boundaries again to recapture any solids moved across
      // grid boundaries during the redistribute
      chem_txfr_ptr[lev]->SumBoundary(gm.periodicity());
      chem_txfr_ptr[lev]->FillBoundary(gm.periodicity());
    }

    // This might not need to exist on all levels. Maybe only level 0.
    for (int lev(0); lev < nlev; ++lev)
      delete tmp_eps[lev];

    int  src_nghost = 1;
    int dest_nghost = 0;
    int ng_to_copy = amrex::min(src_nghost, dest_nghost);

    for (int lev = 1; lev < nlev; lev++) {
      chem_txfr_ptr[0]->copy(*chem_txfr_ptr[lev], 0, 0, chem_txfr_ptr[0]->nComp(),
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
      chem_txfr[0]->copy(*chem_txfr_ptr[0], 0, 0, chem_txfr[0]->nComp());
    }

    for (int lev = 0; lev < nlev; lev++) {
      if (chem_txfr_ptr[lev] != chem_txfr[lev])
        delete chem_txfr_ptr[lev];
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
      chem_txfr[lev]->FillBoundary(geom[lev].periodicity());
  }
}
