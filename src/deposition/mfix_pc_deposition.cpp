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
                          MultiFab & chem_txfr_mf,
                          const MultiFab * volfrac,
                          const amrex::FabArray<EBCellFlagFab>* flags,
                          const FluidPhase& fluid,
                          const int advect_enthalpy)
{
  if (mfix::m_deposition_scheme == DepositionScheme::trilinear) {

    InterphaseTxfrDeposition(TrilinearDeposition(), lev, mf_tmp_eps, txfr_mf,
                             chem_txfr_mf, volfrac, flags, fluid, advect_enthalpy);

  } else if (mfix::m_deposition_scheme == DepositionScheme::square_dpvm) {

    InterphaseTxfrDeposition(TrilinearDPVMSquareDeposition(), lev, mf_tmp_eps,
                             txfr_mf, chem_txfr_mf, volfrac, flags, fluid, advect_enthalpy);

  } else if (mfix::m_deposition_scheme == DepositionScheme::true_dpvm) {

    InterphaseTxfrDeposition(TrueDPVMDeposition(), lev, mf_tmp_eps, txfr_mf,
                             chem_txfr_mf, volfrac, flags, fluid, advect_enthalpy);

  } else if (mfix::m_deposition_scheme == DepositionScheme::centroid) {

    InterphaseTxfrDeposition(CentroidDeposition(), lev, mf_tmp_eps, txfr_mf,
                             chem_txfr_mf, volfrac, flags, fluid, advect_enthalpy);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");

  }
}


template <typename F>
void MFIXParticleContainer::
InterphaseTxfrDeposition (F WeightFunc,
                          int lev,
                          MultiFab & mf_tmp_eps,
                          MultiFab & txfr_mf,
                          MultiFab & chem_txfr_mf,
                          const MultiFab * volfrac,
                          const amrex::FabArray<EBCellFlagFab>* flags,
                          const FluidPhase& fluid,
                          const int advect_enthalpy)
{
  BL_PROFILE("MFIXParticleContainer::InterphaseTxfrDeposition()");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      plo = gm.ProbLoArray();
  const auto      dx  = gm.CellSizeArray();
  const auto      dxi = gm.InvCellSizeArray();

  const auto      reg_cell_vol = dx[0]*dx[1]*dx[2];

  const int nspecies_g = fluid.nspecies;
  const int solve_reactions = reactions.solve;

  const int idx_X_txfr = m_runtimeRealData.species_txfr;
  const int idx_vel_txfr = m_runtimeRealData.vel_txfr;
  const int idx_h_txfr = m_runtimeRealData.h_txfr;

  ChemTransfer chem_txfr_idxs(fluid.nspecies, reactions.nreactions);
  const int idx_Xg_txfr = chem_txfr_idxs.ro_gk_txfr;
  const int idx_velg_txfr = chem_txfr_idxs.vel_g_txfr;
  const int idx_hg_txfr = chem_txfr_idxs.h_g_txfr;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox local_txfr;
    FArrayBox local_chem_txfr;

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

      FArrayBox& eps_fab  = mf_tmp_eps[pti];
      FArrayBox& txfr_fab = txfr_mf[pti];
      FArrayBox& chem_txfr_fab = solve_reactions ? chem_txfr_mf[pti] : dummy_fab;

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered) {

        Array4<Real> dummy_arr;

        auto        txfr_arr = txfr_fab.array();
        auto   chem_txfr_arr = solve_reactions ? chem_txfr_fab.array() : dummy_arr;
        const auto&   volarr = eps_fab.array();
        const auto& flagsarr = (*flags)[pti].array();
        const auto&    vfrac = (*volfrac)[pti].array();

        const Real deposition_scale_factor = mfix::m_deposition_scale_factor;

#ifdef _OPENMP
        Box tile_box = box;

        {
          const int ncomp = txfr_mf.nComp();

          if (Gpu::notInLaunchRegion())
          {
            tile_box.grow(txfr_mf.nGrow());
            local_txfr.resize(tile_box, ncomp);
            local_txfr.setVal<RunOn::Host>(0.0);
            txfr_arr = local_txfr.array();
          }
        }

        if (solve_reactions) {
          const int ncomp = chem_txfr_mf.nComp();

          if (Gpu::notInLaunchRegion())
          {
            tile_box.grow(chem_txfr_mf.nGrow());
            local_chem_txfr.resize(tile_box, ncomp);
            local_chem_txfr.setVal<RunOn::Host>(0.0);
            chem_txfr_arr = local_chem_txfr.array();
          }
        }
#endif

        amrex::ParallelFor(nrp,
            [pstruct,p_realarray,plo,dx,dxi,vfrac,volarr,deposition_scale_factor,
             reg_cell_vol,WeightFunc,flagsarr,txfr_arr,chem_txfr_arr,advect_enthalpy,
             ptile_data,nspecies_g,solve_reactions,idx_X_txfr,idx_vel_txfr,
             idx_h_txfr,idx_Xg_txfr,idx_velg_txfr,idx_hg_txfr,local_cg_dem=DEM::cg_dem]
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

          Real pgamma = advect_enthalpy ?
            statwt * p_realarray[SoArealData::convection][ip] / reg_cell_vol : 0;

          if (local_cg_dem){
             pvol = pvol / statwt;
             pbeta = pbeta / statwt;
          }

          Real pvx = p_realarray[SoArealData::velx][ip] * pbeta;
          Real pvy = p_realarray[SoArealData::vely][ip] * pbeta;
          Real pvz = p_realarray[SoArealData::velz][ip] * pbeta;

          Real pTp = advect_enthalpy ?
            p_realarray[SoArealData::temperature][ip] * pgamma : 0;

          // Chemical reactions deposition terms
          GpuArray<Real,SPECIES::NMAX> X_gk_txfr;
          X_gk_txfr.fill(0.);

          if (solve_reactions) {
            for (int n_g(0); n_g < nspecies_g; ++n_g) {
              X_gk_txfr[n_g] = statwt * ptile_data.m_runtime_rdata[idx_X_txfr+n_g][ip] / reg_cell_vol;
            }
          }

          RealVect vel_g_txfr(0.);
          vel_g_txfr[0]= solve_reactions ?
            statwt * ptile_data.m_runtime_rdata[idx_vel_txfr+0][ip] / reg_cell_vol : 0.;
          vel_g_txfr[1]= solve_reactions ?
            statwt * ptile_data.m_runtime_rdata[idx_vel_txfr+1][ip] / reg_cell_vol : 0.;
          vel_g_txfr[2]= solve_reactions ?
            statwt * ptile_data.m_runtime_rdata[idx_vel_txfr+2][ip] / reg_cell_vol : 0.;

          Real h_g_txfr = solve_reactions ?
            statwt * ptile_data.m_runtime_rdata[idx_h_txfr][ip] / reg_cell_vol : 0.;

          Real pvx_chem = p_realarray[SoArealData::velx][ip] * vel_g_txfr[0];
          Real pvy_chem = p_realarray[SoArealData::vely][ip] * vel_g_txfr[1];
          Real pvz_chem = p_realarray[SoArealData::velz][ip] * vel_g_txfr[2];

          // Deposition
          for (int ii = -1; ii <= 0; ++ii) {
            for (int jj = -1; jj <= 0; ++jj) {
              for (int kk = -1; kk <= 0; ++kk) {
                if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                  continue;

                Real weight_vol = weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk);

                HostDevice::Atomic::Add(&volarr(i+ii,j+jj,k+kk), weight_vol*pvol);

                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,Transfer::velx), weight_vol*pvx);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,Transfer::vely), weight_vol*pvy);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,Transfer::velz), weight_vol*pvz);

                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,Transfer::beta), weight_vol*pbeta);

                if (advect_enthalpy) {
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,Transfer::gammaTp), weight_vol*pTp);
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,Transfer::gamma), weight_vol*pgamma);
                }

                if (solve_reactions) {
                  for (int n_g(0); n_g < nspecies_g; ++n_g) {
                    HostDevice::Atomic::Add(&chem_txfr_arr(i+ii,j+jj,k+kk,idx_Xg_txfr+n_g), weight_vol*X_gk_txfr[n_g]);
                }

                  HostDevice::Atomic::Add(&chem_txfr_arr(i+ii,j+jj,k+kk,idx_velg_txfr+0), weight_vol*pvx_chem);
                  HostDevice::Atomic::Add(&chem_txfr_arr(i+ii,j+jj,k+kk,idx_velg_txfr+1), weight_vol*pvy_chem);
                  HostDevice::Atomic::Add(&chem_txfr_arr(i+ii,j+jj,k+kk,idx_velg_txfr+2), weight_vol*pvz_chem);

                  HostDevice::Atomic::Add(&chem_txfr_arr(i+ii,j+jj,k+kk,idx_hg_txfr), weight_vol*h_g_txfr);
                }
              }
            }
          }
        });

#ifdef _OPENMP
        if (Gpu::notInLaunchRegion())
        {
          {
            const int ncomp = txfr_mf.nComp();
            txfr_fab.atomicAdd<RunOn::Host>(local_txfr, tile_box, tile_box, 0, 0, ncomp);
          }

          if (solve_reactions) {
            const int ncomp = chem_txfr_mf.nComp();
            chem_txfr_fab.atomicAdd<RunOn::Host>(local_chem_txfr, tile_box, tile_box, 0, 0, ncomp);
          }
        }
#endif

      }
    }
  }
}
