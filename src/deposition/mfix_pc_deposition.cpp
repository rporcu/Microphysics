#include <AMReX.H>
#include <AMReX_Particles.H>
#include <mfix_pc.H>

#include <mfix_deposition_K.H>
#include <mfix_dem_parms.H>
#include <mfix_species_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_algorithm.H>

using namespace amrex;

void MFIXParticleContainer::
SolidsVolumeDeposition (int lev,
                  MultiFab & mf_to_be_filled,
                  const MultiFab * volfrac,
                  const amrex::FabArray<EBCellFlagFab>* flags)
{

  if (mfix::m_deposition_scheme == DepositionScheme::trilinear) {

    SolidsVolumeDeposition(TrilinearDeposition(),
                     lev, mf_to_be_filled, volfrac, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::square_dpvm) {

    SolidsVolumeDeposition(TrilinearDPVMSquareDeposition(),
                     lev, mf_to_be_filled, volfrac, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::true_dpvm) {

    SolidsVolumeDeposition(TrueDPVMDeposition(),
                     lev, mf_to_be_filled, volfrac, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::centroid) {

    SolidsVolumeDeposition(CentroidDeposition(),
                     lev, mf_to_be_filled, volfrac, flags);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");

  }


}


template <typename F>
void MFIXParticleContainer::
SolidsVolumeDeposition (F WeightFunc, int lev,
                        MultiFab & mf_to_be_filled,
                        const MultiFab * volfrac,
                        const amrex::FabArray<EBCellFlagFab>* flags)
{
  BL_PROFILE("MFIXParticleContainer::SolidsVolumeDeposition()");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      plo = gm.ProbLoArray();
  const auto      dx  = gm.CellSizeArray();
  const auto      dxi = gm.InvCellSizeArray();

  const auto      reg_cell_vol = dx[0]*dx[1]*dx[2];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox local_fab_to_be_filled;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();
      FArrayBox& fab_to_be_filled = mf_to_be_filled[pti];

      const Box& bx  = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(bx) != FabType::covered ) {

        auto volarr = fab_to_be_filled.array();
        const auto& flagsarr = (*flags)[pti].array();
        const auto& vfrac = (*volfrac)[pti].array();

#ifdef _OPENMP
        const int ncomp = mf_to_be_filled.nComp();
        Box tile_box = bx;

        if(Gpu::notInLaunchRegion())
        {
          tile_box.grow(mf_to_be_filled.nGrow());
          local_fab_to_be_filled.resize(tile_box, ncomp);
          local_fab_to_be_filled.setVal<RunOn::Host>(0.0);
          volarr = local_fab_to_be_filled.array();
        }
#endif

        const Real deposition_scale_factor =
          mfix::m_deposition_scale_factor;

        amrex::ParallelFor(nrp,
          [pstruct,p_realarray,plo,dx,dxi,vfrac,deposition_scale_factor,volarr,
           reg_cell_vol,WeightFunc,flagsarr,local_cg_dem=DEM::cg_dem]
          AMREX_GPU_DEVICE (int ip) noexcept
          {
            const ParticleType& p = pstruct[ip];

            int i;
            int j;
            int k;

            GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

            WeightFunc(plo, dx, dxi, flagsarr, p.pos(), p_realarray[SoArealData::radius][ip], i, j, k, weights,
                deposition_scale_factor);

            Real pvol = p_realarray[SoArealData::statwt][ip] * p_realarray[SoArealData::volume][ip] / reg_cell_vol;

            if (local_cg_dem){
               pvol = pvol / p_realarray[SoArealData::statwt][ip];
            }

            for (int kk = -1; kk <= 0; ++kk) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int ii = -1; ii <= 0; ++ii) {
                  if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                    continue;

                  Real weight_vol = weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk);

                  amrex::Gpu::Atomic::Add(&volarr(i+ii,j+jj,k+kk), weight_vol*pvol);
                }
              }
            }
          });

#ifdef _OPENMP
        if(Gpu::notInLaunchRegion())
        {
          fab_to_be_filled.atomicAdd<RunOn::Host>(local_fab_to_be_filled,
              tile_box, tile_box, 0, 0, ncomp);
        }
#endif

      }
    }
  }
}



void MFIXParticleContainer::
InterphaseTxfrDeposition (int lev,
                          MultiFab & mf_tmp_eps,
                          MultiFab & txfr_mf,
                          const MultiFab * volfrac,
                          const amrex::FabArray<EBCellFlagFab>* flags,
                          const int advect_enthalpy)
{

  if (mfix::m_deposition_scheme == DepositionScheme::trilinear) {

    InterphaseTxfrDeposition(TrilinearDeposition(), lev, mf_tmp_eps, txfr_mf,
        volfrac, flags, advect_enthalpy);

  } else if (mfix::m_deposition_scheme == DepositionScheme::square_dpvm) {

    InterphaseTxfrDeposition(TrilinearDPVMSquareDeposition(), lev, mf_tmp_eps,
        txfr_mf, volfrac, flags, advect_enthalpy);

  } else if (mfix::m_deposition_scheme == DepositionScheme::true_dpvm) {

    InterphaseTxfrDeposition(TrueDPVMDeposition(), lev, mf_tmp_eps, txfr_mf,
        volfrac, flags, advect_enthalpy);

  } else if (mfix::m_deposition_scheme == DepositionScheme::centroid) {

    InterphaseTxfrDeposition(CentroidDeposition(), lev, mf_tmp_eps, txfr_mf,
        volfrac, flags, advect_enthalpy);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");

  }

}



template <typename F>
void MFIXParticleContainer::
InterphaseTxfrDeposition (F WeightFunc, int lev,
                          MultiFab & mf_tmp_eps,
                          MultiFab & txfr_mf,
                          const MultiFab * volfrac,
                          const amrex::FabArray<EBCellFlagFab>* flags,
                          const int advect_enthalpy)
{
  BL_PROFILE("MFIXParticleContainer::InterphaseTxfrDeposition()");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      plo = gm.ProbLoArray();
  const auto      dx  = gm.CellSizeArray();
  const auto      dxi = gm.InvCellSizeArray();

  const auto      reg_cell_vol = dx[0]*dx[1]*dx[2];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox local_txfr;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();

      FArrayBox& eps_fab  = mf_tmp_eps[pti];
      FArrayBox& txfr_fab = txfr_mf[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered ) {

        auto        txfr_arr = txfr_fab.array();
        const auto&   volarr = eps_fab.array();
        const auto& flagsarr = (*flags)[pti].array();
        const auto&    vfrac = (*volfrac)[pti].array();

        const Real deposition_scale_factor =
          mfix::m_deposition_scale_factor;

#ifdef _OPENMP
        const int ncomp = txfr_mf.nComp();
        Box tile_box = box;

        if (Gpu::notInLaunchRegion())
        {
          tile_box.grow(txfr_mf.nGrow());
          local_txfr.resize(tile_box, ncomp);
          local_txfr.setVal<RunOn::Host>(0.0);
          txfr_arr = local_txfr.array();
        }
#endif

        amrex::ParallelFor(nrp,
          [pstruct,p_realarray,plo,dx,dxi,vfrac,volarr,deposition_scale_factor,
           reg_cell_vol,WeightFunc,flagsarr,txfr_arr,advect_enthalpy,
           local_cg_dem=DEM::cg_dem] AMREX_GPU_DEVICE (int ip) noexcept
          {
            const ParticleType& p = pstruct[ip];

            int i;
            int j;
            int k;

            const Real statwt = p_realarray[SoArealData::statwt][ip];

            GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

            WeightFunc(plo, dx, dxi, flagsarr, p.pos(), p_realarray[SoArealData::radius][ip], i, j, k, weights,
                       deposition_scale_factor);

            Real pvol = statwt * p_realarray[SoArealData::volume][ip] / reg_cell_vol;

            Real pbeta = statwt * p_realarray[SoArealData::dragcoeff][ip] / reg_cell_vol;

            Real pgamma = advect_enthalpy ?
              statwt * p_realarray[SoArealData::convection][ip] / reg_cell_vol : 0;

            if (local_cg_dem){
               pvol = pvol / statwt;
               pbeta = pbeta / statwt;
            }

            Real pvx   = p_realarray[SoArealData::velx][ip] * pbeta;
            Real pvy   = p_realarray[SoArealData::vely][ip] * pbeta;
            Real pvz   = p_realarray[SoArealData::velz][ip] * pbeta;

            Real pTp   = advect_enthalpy ?
              p_realarray[SoArealData::temperature][ip] * pgamma : 0;

            for (int ii = -1; ii <= 0; ++ii) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int kk = -1; kk <= 0; ++kk) {
                  if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                    continue;

                  Real weight_vol = weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk);

                  amrex::Gpu::Atomic::Add(&volarr(i+ii,j+jj,k+kk), weight_vol*pvol);

                  amrex::Gpu::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,Transfer::velx), weight_vol*pvx);
                  amrex::Gpu::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,Transfer::vely), weight_vol*pvy);
                  amrex::Gpu::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,Transfer::velz), weight_vol*pvz);

                  amrex::Gpu::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,Transfer::beta), weight_vol*pbeta);

                  if (advect_enthalpy)
                  {
                    amrex::Gpu::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,Transfer::gammaTp), weight_vol*pTp);
                    amrex::Gpu::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,Transfer::gamma), weight_vol*pgamma);
                  }
                }
              }
            }
          });

#ifdef _OPENMP
        if (Gpu::notInLaunchRegion())
        {
          txfr_fab.atomicAdd<RunOn::Host>(local_txfr, tile_box, tile_box, 0, 0, ncomp);
        }
#endif

      }
    }
  }
}


void MFIXParticleContainer::
InterphaseChemDeposition (int lev,
                          MultiFab & mf_tmp_eps,
                          MultiFab & Rrates,
                          const MultiFab * volfrac,
                          const FabArray<EBCellFlagFab>* flags,
                          const Vector< REACTIONS::ChemicalReaction* >& chemical_reactions)
{
  if (mfix::m_deposition_scheme == DepositionScheme::trilinear)
  {
    InterphaseChemDeposition(TrilinearDeposition(), lev, mf_tmp_eps, Rrates,
                             volfrac, flags, chemical_reactions);
  }
  else if (mfix::m_deposition_scheme == DepositionScheme::square_dpvm)
  {
    InterphaseChemDeposition(TrilinearDPVMSquareDeposition(), lev, mf_tmp_eps,
                             Rrates, volfrac, flags, chemical_reactions);
  }
  else if (mfix::m_deposition_scheme == DepositionScheme::true_dpvm)
  {
    InterphaseChemDeposition(TrueDPVMDeposition(), lev, mf_tmp_eps, Rrates,
                             volfrac, flags, chemical_reactions);
  }
  else if (mfix::m_deposition_scheme == DepositionScheme::centroid)
  {
    InterphaseChemDeposition(CentroidDeposition(), lev, mf_tmp_eps, Rrates,
                             volfrac, flags, chemical_reactions);
  }
  else
  {
    amrex::Abort("Don't know this deposition_scheme!");
  }
}


template <typename F>
void MFIXParticleContainer::
InterphaseChemDeposition (F WeightFunc,
                          int lev,
                          MultiFab& mf_tmp_eps,
                          MultiFab& mf_G_gk_fp,
                          const MultiFab* volfrac,
                          const FabArray<EBCellFlagFab>* flags,
                          const Vector< REACTIONS::ChemicalReaction* >& chemical_reactions)
{
  BL_PROFILE("MFIXParticleContainer::InterphaseChemDeposition()");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto plo = gm.ProbLoArray();
  const auto dx  = gm.CellSizeArray();
  const auto dxi = gm.InvCellSizeArray();

  const auto reg_cell_vol = dx[0]*dx[1]*dx[2];

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
    h_types[q] = chemical_reactions[q]->m_reaction_type;
    h_phases[q] = chemical_reactions[q]->m_phases.data();
    h_nphases[q] = chemical_reactions[q]->m_phases.size();
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
    h_nreactants[q] = chemical_reactions[q]->m_reactants.size();
    h_reactants_id[q] = chemical_reactions[q]->m_reactants_id.data();
    h_reactants_coeffs[q] = chemical_reactions[q]->m_reactants_coeffs.data();
    h_reactants_phases[q] = chemical_reactions[q]->m_reactants_phases.data();
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
    h_nproducts[q] = chemical_reactions[q]->m_products.size();
    h_products_id[q] = chemical_reactions[q]->m_products_id.data();
    h_products_coeffs[q] = chemical_reactions[q]->m_products_coeffs.data();
    h_products_phases[q] = chemical_reactions[q]->m_products_phases.data();
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

  // Particles indexes
  const int idx_G = SoAspeciesData::count*nspecies_s + SoAreactionsData::G_sn_pg_q*nreactions;

  Gpu::synchronize();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox local_fab_txfr;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      //Access to added variables
      PairIndex index(pti.index(), pti.LocalTileIndex());
      auto& ptile = GetParticles(lev)[index];
      auto ptile_data = ptile.getParticleTileData();

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();
      const long nrp = pti.numParticles();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      FArrayBox& fab_G_gk_fp = mf_G_gk_fp[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered )
      {
        auto     arr_G_gk_fp = fab_G_gk_fp.array();
        const auto& flagsarr = (*flags)[pti].array();
        const auto&    vfrac = (*volfrac)[pti].array();

        const Real deposition_scale_factor =
          mfix::m_deposition_scale_factor;

#ifdef _OPENMP
        const int ncomp = mf_G_gk_fp.nComp();
        Box tile_box = box;

        if (Gpu::notInLaunchRegion())
        {
          tile_box.grow(mf_G_gk_fp.nGrow());
          local_fab_txfr.resize(tile_box, ncomp);
          local_fab_txfr.setVal<RunOn::Host>(0.0);
          arr_G_gk_fp = local_fab_txfr.array();
        }
#endif

        amrex::ParallelFor(nrp,
          [pstruct,p_realarray,plo,dx,dxi,vfrac,deposition_scale_factor,reg_cell_vol,
           WeightFunc,flagsarr,arr_G_gk_fp,idx_G,nspecies_s,nreactions,
           ptile_data,nrp,p_species_id_s,p_species_id_g,nspecies_g,p_nproducts,
           p_nreactants,p_reactants_id,p_reactants_coeffs,p_reactants_phases,
           p_products_id,p_products_coeffs,p_products_phases,p_types,p_phases,
           p_nphases,Solid,Fluid,p_MW_sn,p_MW_gk,Heterogeneous,InvalidIdx]
        AMREX_GPU_DEVICE (int p_id) noexcept
        {
          const ParticleType& p = pstruct[p_id];

          int i(0); int j(0); int k(0);

          GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

          WeightFunc(plo, dx, dxi, flagsarr, p.pos(), p_realarray[SoArealData::radius][p_id], i, j, k, weights,
                     deposition_scale_factor);

          // Pointer to this particle's species rate of formation
          GpuArray<GpuArray<Real,REACTIONS::NMAX>,SPECIES::NMAX> G_sn_pg_q;
          
          for (int n_s(0); n_s < nspecies_s; n_s++) {
            for (int q(0); q < nreactions; q++)
              G_sn_pg_q[n_s][q] = 0;
          } 
          
          for (int n_s(0); n_s < nspecies_s; n_s++) {
            for (int q(0); q < nreactions; q++) {
              const int idx = idx_G + n_s*nreactions + q;
              G_sn_pg_q[n_s][q] = ptile_data.m_runtime_rdata[idx][p_id];
            }
          }

          // Create R_q for storing chemical reactions Rates
          GpuArray<Real,REACTIONS::NMAX> R_q;

          for (int q(0); q < nreactions; q++)
            R_q[q] = 0;

          // Create G_gk_fp_q for computing fluid density rate of change
          GpuArray<Real,SPECIES::NMAX> G_gk_fp;

          for (int n_g(0); n_g < nspecies_g; n_g++)
            G_gk_fp[n_g] = 0;

          for (int q(0); q < nreactions; q++)
          {
            // Do something only if reaction is heterogeneous and contains
            // a solid compound
            if (p_types[q] == Heterogeneous and
                MFIXfind(p_phases[q], p_nphases[q], Solid) != InvalidIdx)
            {
              Real stoc_coeff(0);
              int current_species_id(InvalidIdx);

              // Get the position of the first reactant which is Solid
              const int n_s_react = MFIXfind(p_reactants_phases[q], p_nreactants[q], Solid);

              // If a Solid reactant it was found (it may happen to not find it)
              if (n_s_react != InvalidIdx)
              {
                // Get the unique identifier ID for the first solid species
                // among the reactants of current reaction q
                current_species_id = p_reactants_id[q][n_s_react];

                // Add the contribution to the total stoichiometric coefficient
                stoc_coeff += p_reactants_coeffs[q][n_s_react];

                // Then look for the same compound among the products as well
                const int n_s_prod = MFIXfind(p_products_id[q], p_nproducts[q],
                    current_species_id);

                // If the current species is found among the products
                if (n_s_prod != InvalidIdx) {
                  // If the found species is also Solid
                  if (p_products_phases[q][n_s_prod] == Solid) {
                    // Add the contribution to the total stoichiometric coefficient
                    stoc_coeff += p_products_coeffs[q][n_s_prod];
                  }
                }
              }
              else { // No solid species among reactants, look among products
                // Get the position of the first product which is Solid
                const int n_s_prod = MFIXfind(p_products_phases[q], p_nproducts[q], Solid);

                BL_ASSERT(n_s_prod != InvalidIdx);

                current_species_id = p_products_id[q][n_s_prod];

                // Add the contribution to the total stoichiometric coefficient
                stoc_coeff += p_products_coeffs[q][n_s_prod];
              }

              // Get the index of the current solid species
              const int n_s = MFIXfind(p_species_id_s, nspecies_s, current_species_id);

              // Extract the density reaction rate for the current solid species
              // n_s and the current reaction q
              R_q[q] = G_sn_pg_q[n_s][q] / (stoc_coeff*p_MW_sn[n_s]);
            }
            else {
              R_q[q] = 0;
            }
          }

          // Compute fluid species densities reaction rates
          for (int n_g(0); n_g < nspecies_g; n_g++)
          {
            const int current_species_id = p_species_id_g[n_g];

            for (int q(0); q < nreactions; q++)
            {
              // Do something only if reaction is heterogeneous and contains
              // a fluid compound
              if (p_types[q] == Heterogeneous and
                  MFIXfind(p_phases[q], p_nphases[q], Fluid) != InvalidIdx)
              {
                Real stoc_coeff(0);

                // 
                const int n_g_react = 
                  MFIXfind(p_reactants_id[q], p_nreactants[q], current_species_id);

                // If current species is among the reactants
                if (n_g_react != InvalidIdx)
                {
                  // IF the current species phase is Fluid
                  if (p_reactants_phases[q][n_g_react] == Fluid)
                    // 
                    stoc_coeff += p_reactants_coeffs[q][n_g_react];
                }

                // 
                const int n_g_prod = 
                  MFIXfind(p_products_id[q], p_nproducts[q], current_species_id);

                // If current species is among the products
                if (n_g_prod != InvalidIdx)
                {
                  // IF the current species phase is Fluid
                  if (p_products_phases[q][n_g_prod] == Fluid)
                    // 
                    stoc_coeff += p_products_coeffs[q][n_g_prod];
                }

                G_gk_fp[n_g] += stoc_coeff * p_MW_gk[n_g] * R_q[q];
              }
            }
          }

          // Deposition
          for (int ii = -1; ii <= 0; ++ii) {
            for (int jj = -1; jj <= 0; ++jj) {
              for (int kk = -1; kk <= 0; ++kk) {
                if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                  continue;

                Real weight_vol = weights[ii+1][jj+1][kk+1] /
                                         vfrac(i+ii,j+jj,k+kk);

                for (int n_g(0); n_g < nspecies_g; n_g++)
                {
                  // Deposition of Rrates
                  amrex::Gpu::Atomic::Add(&arr_G_gk_fp(i+ii,j+jj,k+kk,n_g),
                                          weight_vol*G_gk_fp[n_g]);
                }
              }
            }
          }
        });

#ifdef _OPENMP
        if (Gpu::notInLaunchRegion())
        {
          fab_G_gk_fp.atomicAdd<RunOn::Host>(local_fab_txfr, tile_box, tile_box,
              0, 0, ncomp);
        }
#endif

      }
    }
  }
}
