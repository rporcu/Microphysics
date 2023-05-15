#include <mfix_deposition_op.H>
#include <mfix_pc.H>
#include <mfix_deposition_K.H>
#include <mfix_dem.H>
#include <mfix_species.H>
#include <mfix_fluid.H>
#include <mfix_reactions.H>
#include <mfix_algorithm.H>

#include <AMReX.H>
#include <AMReX_Particles.H>


using namespace amrex;


void
MFIXSolidsVolume::deposit (int lev,
                           const Geometry& geom,
                           MFIXParticleContainer* pc,
                           const MultiFab* volfrac,
                           const amrex::FabArray<EBCellFlagFab>* flags,
                           MultiFab* mf_to_be_filled,
                           MultiFab*)
{
  if (mfix::m_deposition_scheme == DepositionScheme::trilinear) {

    deposit(TrilinearDeposition(), lev, geom, pc, volfrac, flags, mf_to_be_filled);

  } else if (mfix::m_deposition_scheme == DepositionScheme::square_dpvm) {

    deposit(TrilinearDPVMSquareDeposition(), lev, geom, pc, volfrac, flags, mf_to_be_filled);

  } else if (mfix::m_deposition_scheme == DepositionScheme::true_dpvm) {

    deposit(TrueDPVMDeposition(), lev, geom, pc, volfrac, flags, mf_to_be_filled);

  } else if (mfix::m_deposition_scheme == DepositionScheme::centroid) {

    deposit(CentroidDeposition(), lev, geom, pc, volfrac, flags, mf_to_be_filled);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");

  }
}


template <typename F>
void
MFIXSolidsVolume::deposit (F WeightFunc,
                           int lev,
                           const Geometry& geom,
                           MFIXParticleContainer* pc,
                           const MultiFab* volfrac,
                           const amrex::FabArray<EBCellFlagFab>* flags,
                           MultiFab* mf_to_be_filled,
                           MultiFab*)
{
  BL_PROFILE("MFIXParticleContainer::SolidsVolumeDeposition()");

  // We always use the coarse dx
  const auto      plo = geom.ProbLoArray();
  const auto      dx  = geom.CellSizeArray();
  const auto      dxi = geom.InvCellSizeArray();

  const auto      reg_cell_vol = dx[0]*dx[1]*dx[2];

  const int solve_pic = pc->get_pic().solve();
  const int cg_dem = pc->get_dem().cg_dem();
  const int idx_statwt = pc->m_runtimeRealData.statwt;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox local_fab_to_be_filled;

    for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {

      MFIXParticleContainer::PairIndex index(pti.index(), pti.LocalTileIndex());
      auto& plev  = pc->GetParticles(lev);
      auto& ptile = plev[index];
      auto ptile_data = ptile.getParticleTileData();

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();
      FArrayBox& fab_to_be_filled = (*mf_to_be_filled)[pti];

      const Box& bx  = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(bx) != FabType::covered ) {

        auto volarr = fab_to_be_filled.array();
        const auto& flagsarr = (*flags)[pti].array();
        const auto& vfrac = (*volfrac)[pti].array();

#ifdef _OPENMP
        const int ncomp = mf_to_be_filled->nComp();
        Box tile_box = bx;

        if(Gpu::notInLaunchRegion())
        {
          tile_box.grow(mf_to_be_filled->nGrow());
          local_fab_to_be_filled.resize(tile_box, ncomp);
          local_fab_to_be_filled.setVal<RunOn::Host>(0.0);
          volarr = local_fab_to_be_filled.array();
        }
#endif

        const Real deposition_scale_factor = mfix::m_deposition_scale_factor;

        const auto local_cg_dem = pc->get_dem().cg_dem();

        amrex::ParallelFor(nrp,
            [pstruct,p_realarray,plo,dx,dxi,vfrac,deposition_scale_factor,volarr,
             reg_cell_vol,WeightFunc,flagsarr,local_cg_dem,ptile_data,idx_statwt,
             solve_pic,cg_dem]
          AMREX_GPU_DEVICE (int ip) noexcept
        {
          const ParticleType& p = pstruct[ip];

          int i;
          int j;
          int k;

          GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

          const Real pradius = p_realarray[SoArealData::radius][ip];

          WeightFunc(plo, dx, dxi, flagsarr, p.pos(), pradius, i, j, k, weights,
              deposition_scale_factor);

          const Real pvolume = SoArealData::volume(pradius);

          Real pvol = pvolume / reg_cell_vol;

          if (solve_pic) {
            pvol *= ptile_data.m_runtime_rdata[idx_statwt][ip];
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


void
MFIXInterphaseTxfr::deposit (int lev,
                             const Geometry& geom,
                             MFIXParticleContainer* pc,
                             const MultiFab* volfrac,
                             const amrex::FabArray<EBCellFlagFab>* flags,
                             MultiFab* txfr_mf,
                             MultiFab* eps_mf)
{
  if (mfix::m_deposition_scheme == DepositionScheme::trilinear) {

    deposit(TrilinearDeposition(), lev, geom, pc, volfrac, flags, txfr_mf, eps_mf);

  } else if (mfix::m_deposition_scheme == DepositionScheme::square_dpvm) {

    deposit(TrilinearDPVMSquareDeposition(), lev, geom, pc, volfrac, flags, txfr_mf, eps_mf);

  } else if (mfix::m_deposition_scheme == DepositionScheme::true_dpvm) {

    deposit(TrueDPVMDeposition(), lev, geom, pc, volfrac, flags, txfr_mf, eps_mf);

  } else if (mfix::m_deposition_scheme == DepositionScheme::centroid) {

    deposit(CentroidDeposition(), lev, geom, pc, volfrac, flags, txfr_mf, eps_mf);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");

  }
}


template <typename F>
void
MFIXInterphaseTxfr::deposit (F WeightFunc,
                             int lev,
                             const Geometry& geom,
                             MFIXParticleContainer* pc,
                             const MultiFab* volfrac,
                             const amrex::FabArray<EBCellFlagFab>* flags,
                             MultiFab* txfr_mf,
                             MultiFab* eps_mf)
{
  BL_PROFILE("MFIXParticleContainer::InterphaseTxfrDeposition()");

  // We always use the coarse dx
  const auto plo = geom.ProbLoArray();
  const auto dx  = geom.CellSizeArray();
  const auto dxi = geom.InvCellSizeArray();

  const auto reg_cell_vol = dx[0]*dx[1]*dx[2];

  const int nspecies_g = pc->get_fluid().nspecies();
  const int solve_reactions = pc->get_reactions().solve();

  const int idx_mass_txfr = pc->m_runtimeRealData.mass_txfr;
  const int idx_vel_txfr = pc->m_runtimeRealData.vel_txfr;
  const int idx_h_txfr = pc->m_runtimeRealData.h_txfr;
  const int idx_statwt = pc->m_runtimeRealData.statwt;

  InterphaseTxfrIndexes txfr_idxs(pc->get_fluid().nspecies(), pc->get_reactions().nreactions());

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

    for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {

      PairIndex index(pti.index(), pti.LocalTileIndex());
      auto& ptile = pc->GetParticles(lev)[index];

      //Access to added variables
      auto ptile_data = ptile.getParticleTileData();

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();

      FArrayBox dummy_fab;

      FArrayBox& eps_fab  = (*eps_mf)[pti];
      FArrayBox& txfr_fab = (*txfr_mf)[pti];

      auto aux_iterator = m_aux[lev].find(index);

      AMREX_ASSERT(aux_iterator != m_aux[lev].end());

      auto& aux_vector = aux_iterator->second;
      Real* aux_ptr = aux_vector.dataPtr();

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
          const int ncomp = eps_mf->nComp();

          if (Gpu::notInLaunchRegion()) {
            eps_tile_box.grow(eps_mf->nGrow());
            local_eps.resize(eps_tile_box, ncomp);
            local_eps.setVal<RunOn::Host>(0.0);
            eps_arr = local_eps.array();
          }
        }

        Box txfr_tile_box = box;
        {
          const int ncomp = txfr_mf->nComp();

          if (Gpu::notInLaunchRegion()) {
            txfr_tile_box.grow(txfr_mf->nGrow());
            local_txfr.resize(txfr_tile_box, ncomp);
            local_txfr.setVal<RunOn::Host>(0.0);
            txfr_arr = local_txfr.array();
          }
        }
#endif

        const int solve_enthalpy = pc->get_fluid().solve_enthalpy();
        const int solve_pic = pc->get_pic().solve();
        const int cg_dem = pc->get_dem().cg_dem();

        const auto local_cg_dem = pc->get_dem().cg_dem();

        amrex::ParallelFor(nrp,
            [pstruct,p_realarray,plo,dx,dxi,vfrac,eps_arr,deposition_scale_factor,
             nrp,reg_cell_vol,WeightFunc,flagsarr,txfr_arr,solve_enthalpy,
             ptile_data,nspecies_g,solve_reactions,idx_mass_txfr,idx_vel_txfr,aux_ptr,
             idx_h_txfr,idx_Xg_txfr,idx_velg_txfr,idx_hg_txfr,idx_velx_txfr,
             idx_vely_txfr,idx_velz_txfr,idx_drag_txfr,idx_gammaTp_txfr,cg_dem,
             idx_convection_coeff_txfr,local_cg_dem,idx_statwt,solve_pic]
          AMREX_GPU_DEVICE (int ip) noexcept
        {
          const ParticleType& p = pstruct[ip];

          int i;
          int j;
          int k;

          GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

          Real pradius = p_realarray[SoArealData::radius][ip];

          WeightFunc(plo, dx, dxi, flagsarr, p.pos(), pradius, i, j, k, weights,
                     deposition_scale_factor);

          Real pvol = SoArealData::volume(pradius) / reg_cell_vol;
          Real pbeta = p_realarray[SoArealData::dragcoeff][ip] / reg_cell_vol;

          if (solve_pic) {
            pvol *= ptile_data.m_runtime_rdata[idx_statwt][ip];
            pbeta *= ptile_data.m_runtime_rdata[idx_statwt][ip];
          }

          Real pgamma(0.);

          if (solve_enthalpy) {
            pgamma = p_realarray[SoArealData::convection][ip] / reg_cell_vol;

            if (solve_pic)
              pgamma *= ptile_data.m_runtime_rdata[idx_statwt][ip];
          }

          Real pvx = p_realarray[SoArealData::velx][ip] * pbeta;
          Real pvy = p_realarray[SoArealData::vely][ip] * pbeta;
          Real pvz = p_realarray[SoArealData::velz][ip] * pbeta;

          Real pTp(0.);

          if (solve_enthalpy)
            pTp = p_realarray[SoArealData::temperature][ip] * pgamma;

          // Chemical reactions deposition terms
          if (solve_reactions) {
            for (int n_g(0); n_g < nspecies_g; ++n_g) {
              aux_ptr[n_g*nrp + ip] = aux_ptr[n_g*nrp + ip] / reg_cell_vol;
            }

            if (solve_pic)
              for (int n_g(0); n_g < nspecies_g; ++n_g)
                aux_ptr[n_g*nrp + ip] *= ptile_data.m_runtime_rdata[idx_statwt][ip];
          }

          Real velx_chem_txfr(0.);
          Real vely_chem_txfr(0.);
          Real velz_chem_txfr(0.);

          Real h_chem_txfr(0);

          if (solve_reactions) {
            // Note: ptile_data.m_runtime_rdata[idx_vel_txfr][ip] currently
            // contains G_m_p_heterogeneous
            const Real G_m_g_heterogeneous = -1*ptile_data.m_runtime_rdata[idx_vel_txfr][ip];
            const Real coeff = amrex::max(0., G_m_g_heterogeneous);

            Real psigma = coeff / reg_cell_vol;

            if (solve_pic)
              psigma *= ptile_data.m_runtime_rdata[idx_statwt][ip];

            velx_chem_txfr = p_realarray[SoArealData::velx][ip] * psigma;
            vely_chem_txfr = p_realarray[SoArealData::vely][ip] * psigma;
            velz_chem_txfr = p_realarray[SoArealData::velz][ip] * psigma;

            // Note: ptile_data.m_runtime_rdata[idx_h_txfr][ip] currently
            // contains the opposite of G_H_g_heterogeneous
            const Real G_H_g_heterogeneous = -1*ptile_data.m_runtime_rdata[idx_h_txfr][ip];

            h_chem_txfr = G_H_g_heterogeneous / reg_cell_vol;

            if (solve_pic)
              h_chem_txfr *= ptile_data.m_runtime_rdata[idx_statwt][ip];
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
                    HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_Xg_txfr+n_g), weight_vol*aux_ptr[n_g*nrp + ip]);
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
            const int ncomp = eps_mf->nComp();
            eps_fab.atomicAdd<RunOn::Host>(local_eps, eps_tile_box, eps_tile_box, 0, 0, ncomp);
        }

        if (Gpu::notInLaunchRegion()) {
            const int ncomp = txfr_mf->nComp();
            txfr_fab.atomicAdd<RunOn::Host>(local_txfr, txfr_tile_box, txfr_tile_box, 0, 0, ncomp);
        }
#endif
      }
    }
  }
}
