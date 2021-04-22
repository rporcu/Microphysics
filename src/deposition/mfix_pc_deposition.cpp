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
