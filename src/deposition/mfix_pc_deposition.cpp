#include <AMReX.H>
#include <AMReX_Particles.H>
#include <mfix_pc.H>

#include <mfix_deposition_K.H>
#include <mfix_dem.H>
#include <mfix_species.H>
#include <mfix_fluid.H>
#include <mfix_reactions.H>
#include <mfix_algorithm.H>

using namespace amrex;

void MFIXParticleContainer::
SolidsVolumeDeposition (int lev,
                        MultiFab & mf_to_be_filled,
                        const MultiFab * volfrac,
                        const amrex::FabArray<EBCellFlagFab>* flags)
{
  if (mfix::m_deposition_scheme == DepositionScheme::trilinear) {

    SolidsVolumeDeposition(TrilinearDeposition(), lev, mf_to_be_filled, volfrac, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::square_dpvm) {

    SolidsVolumeDeposition(TrilinearDPVMSquareDeposition(), lev, mf_to_be_filled, volfrac, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::true_dpvm) {

    SolidsVolumeDeposition(TrueDPVMDeposition(), lev, mf_to_be_filled, volfrac, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::centroid) {

    SolidsVolumeDeposition(CentroidDeposition(), lev, mf_to_be_filled, volfrac, flags);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");

  }
}


template <typename F>
void MFIXParticleContainer::
SolidsVolumeDeposition (F WeightFunc,
                        int lev,
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
           reg_cell_vol,WeightFunc,flagsarr,local_cg_dem=m_dem.cg_dem()]
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

                  HostDevice::Atomic::Add(&volarr(i+ii,j+jj,k+kk), weight_vol*pvol);
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
                          std::map<PairIndex, Gpu::DeviceVector<Real>>& aux)
{
  if (mfix::m_deposition_scheme == DepositionScheme::trilinear) {

    InterphaseTxfrDeposition(TrilinearDeposition(), lev, mf_tmp_eps, txfr_mf,
                             volfrac, flags, aux);

  } else if (mfix::m_deposition_scheme == DepositionScheme::square_dpvm) {

    InterphaseTxfrDeposition(TrilinearDPVMSquareDeposition(), lev, mf_tmp_eps,
                             txfr_mf, volfrac, flags, aux);

  } else if (mfix::m_deposition_scheme == DepositionScheme::true_dpvm) {

    InterphaseTxfrDeposition(TrueDPVMDeposition(), lev, mf_tmp_eps, txfr_mf,
                             volfrac, flags, aux);

  } else if (mfix::m_deposition_scheme == DepositionScheme::centroid) {

    InterphaseTxfrDeposition(CentroidDeposition(), lev, mf_tmp_eps, txfr_mf,
                             volfrac, flags, aux);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");

  }
}


template <typename F>
void MFIXParticleContainer::
InterphaseTxfrDeposition (F WeightFunc,
                          int lev,
                          MultiFab & eps_mf,
                          MultiFab & txfr_mf,
                          const MultiFab * volfrac,
                          const amrex::FabArray<EBCellFlagFab>* flags,
                          std::map<PairIndex, Gpu::DeviceVector<Real>>& aux)
{
  BL_PROFILE("MFIXParticleContainer::InterphaseTxfrDeposition()");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      plo = gm.ProbLoArray();
  const auto      dx  = gm.CellSizeArray();
  const auto      dxi = gm.InvCellSizeArray();

  const auto      reg_cell_vol = dx[0]*dx[1]*dx[2];

  const int nspecies_g = fluid.nspecies();
  const int solve_reactions = reactions.solve();

  const int idx_mass_txfr = m_runtimeRealData.mass_txfr;
  const int idx_vel_txfr = m_runtimeRealData.vel_txfr;
  const int idx_h_txfr = m_runtimeRealData.h_txfr;

  InterphaseTxfrIndexes txfr_idxs(fluid.nspecies(), reactions.nreactions());

  const int idx_velx_txfr = txfr_idxs.vel+0;
  const int idx_vely_txfr = txfr_idxs.vel+1;
  const int idx_velz_txfr = txfr_idxs.vel+2;
  const int idx_drag_txfr = txfr_idxs.drag_coeff;
  const int idx_gammaTp_txfr = txfr_idxs.gammaTp;
  const int idx_convection_coeff_txfr = txfr_idxs.convection_coeff;

  const int idx_Xg_txfr   = txfr_idxs.chem_ro_gk;
  const int idx_velg_txfr = txfr_idxs.chem_vel;
  const int idx_hg_txfr   = txfr_idxs.chem_h;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox local_eps;
    FArrayBox local_txfr;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

      PairIndex index(pti.index(), pti.LocalTileIndex());
      auto& ptile = GetParticles(lev)[index];

      //Access to added variables
      auto ptile_data = ptile.getParticleTileData();

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();

      FArrayBox dummy_fab;

      FArrayBox& eps_fab  = eps_mf[pti];
      FArrayBox& txfr_fab = txfr_mf[pti];

      Real* aux_ptr = aux[index].dataPtr();

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered) {

        auto        txfr_arr = txfr_fab.array();
        auto         eps_arr = eps_fab.array();
        const auto& flagsarr = (*flags)[pti].const_array();
        const auto&    vfrac = (*volfrac)[pti].const_array();

        const Real deposition_scale_factor = mfix::m_deposition_scale_factor;

#ifdef _OPENMP
        Box eps_tile_box = box;
        {
          const int ncomp = eps_mf.nComp();

          if (Gpu::notInLaunchRegion()) {
            eps_tile_box.grow(eps_mf.nGrow());
            local_eps.resize(eps_tile_box, ncomp);
            local_eps.setVal<RunOn::Host>(0.0);
            eps_arr = local_eps.array();
          }
        }

        Box txfr_tile_box = box;
        {
          const int ncomp = txfr_mf.nComp();

          if (Gpu::notInLaunchRegion()) {
            txfr_tile_box.grow(txfr_mf.nGrow());
            local_txfr.resize(txfr_tile_box, ncomp);
            local_txfr.setVal<RunOn::Host>(0.0);
            txfr_arr = local_txfr.array();
          }
        }
#endif

        const int solve_enthalpy = fluid.solve_enthalpy();

        amrex::ParallelFor(nrp,
            [pstruct,p_realarray,plo,dx,dxi,vfrac,eps_arr,deposition_scale_factor,
             nrp,reg_cell_vol,WeightFunc,flagsarr,txfr_arr,solve_enthalpy,
             ptile_data,nspecies_g,solve_reactions,idx_mass_txfr,idx_vel_txfr,aux_ptr,
             idx_h_txfr,idx_Xg_txfr,idx_velg_txfr,idx_hg_txfr, idx_velx_txfr,
             idx_vely_txfr, idx_velz_txfr, idx_drag_txfr, idx_gammaTp_txfr,
             idx_convection_coeff_txfr, local_cg_dem=m_dem.cg_dem()]
          AMREX_GPU_DEVICE (int ip) noexcept
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

          Real pgamma(0.);

          if (solve_enthalpy)
            pgamma = statwt * p_realarray[SoArealData::convection][ip] / reg_cell_vol;

          if (local_cg_dem) {
            pvol = pvol / statwt;
            pbeta = pbeta / statwt;
          }

          Real pvx = p_realarray[SoArealData::velx][ip] * pbeta;
          Real pvy = p_realarray[SoArealData::vely][ip] * pbeta;
          Real pvz = p_realarray[SoArealData::velz][ip] * pbeta;

          Real pTp(0.);

          if (solve_enthalpy)
            pTp = p_realarray[SoArealData::temperature][ip] * pgamma;

          // Chemical reactions deposition terms
          GpuArray<Real, MFIXSpecies::NMAX> ro_chem_txfr;
          ro_chem_txfr.fill(0.);

          if (solve_reactions) {
            for (int n_g(0); n_g < nspecies_g; ++n_g) {
              ro_chem_txfr[n_g] = statwt * aux_ptr[n_g*nrp + ip] / reg_cell_vol;
            }
          }

          Real velx_chem_txfr(0.);
          Real vely_chem_txfr(0.);
          Real velz_chem_txfr(0.);

          Real h_chem_txfr(0);

          if (solve_reactions) {
            Real psigma = statwt * ptile_data.m_runtime_rdata[idx_vel_txfr][ip] / reg_cell_vol;

            velx_chem_txfr = p_realarray[SoArealData::velx][ip] * psigma;
            vely_chem_txfr = p_realarray[SoArealData::vely][ip] * psigma;
            velz_chem_txfr = p_realarray[SoArealData::velz][ip] * psigma;

            h_chem_txfr = statwt * ptile_data.m_runtime_rdata[idx_h_txfr][ip] / reg_cell_vol;
          }

          // Deposition
          for (int ii = -1; ii <= 0; ++ii) {
            for (int jj = -1; jj <= 0; ++jj) {
              for (int kk = -1; kk <= 0; ++kk) {
                if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                  continue;

                Real weight_vol = weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk);

                HostDevice::Atomic::Add(&eps_arr(i+ii,j+jj,k+kk), weight_vol*pvol);

                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_velx_txfr), weight_vol*pvx);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vely_txfr), weight_vol*pvy);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_velz_txfr), weight_vol*pvz);

                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_drag_txfr), weight_vol*pbeta);

                if (solve_enthalpy) {
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_gammaTp_txfr         ), weight_vol*pTp);
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_convection_coeff_txfr), weight_vol*pgamma);
                }

                if (solve_reactions) {
                  for (int n_g(0); n_g < nspecies_g; ++n_g) {
                    HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_Xg_txfr+n_g), weight_vol*ro_chem_txfr[n_g]);
                  }

                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_velg_txfr+0), weight_vol*velx_chem_txfr);
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_velg_txfr+1), weight_vol*vely_chem_txfr);
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_velg_txfr+2), weight_vol*velz_chem_txfr);

                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_hg_txfr), weight_vol*h_chem_txfr);
                }
              }
            }
          }
        });

#ifdef _OPENMP
        if (Gpu::notInLaunchRegion()) {
            const int ncomp = eps_mf.nComp();
            eps_fab.atomicAdd<RunOn::Host>(local_eps, eps_tile_box, eps_tile_box, 0, 0, ncomp);
        }

        if (Gpu::notInLaunchRegion()) {
            const int ncomp = txfr_mf.nComp();
            txfr_fab.atomicAdd<RunOn::Host>(local_txfr, txfr_tile_box, txfr_tile_box, 0, 0, ncomp);
        }
#endif
      }
    }
  }
}
