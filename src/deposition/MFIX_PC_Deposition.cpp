#include <AMReX.H>
#include "AMReX_Particles.H"
#include <MFIXParticleContainer.H>

#include <mfix_deposition_K.H>
#include <mfix.H>

using namespace amrex;
using namespace std;

void MFIXParticleContainer::
ScalarDeposition(int lev,
                 amrex::MultiFab & mf_to_be_filled,
                 const amrex::MultiFab * volfrac,
                 const amrex::FabArray<EBCellFlagFab>* flags)
{

  if (mfix::m_deposition_scheme == DepositionScheme::trilinear) {

    ScalarDeposition(TrilinearDeposition(),
                     lev, mf_to_be_filled, volfrac, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::square_dpvm) {

    ScalarDeposition(TrilinearDPVMSquareDeposition(),
                     lev, mf_to_be_filled, volfrac, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::true_dpvm) {

    ScalarDeposition(TrueDPVMDeposition(),
                     lev, mf_to_be_filled, volfrac, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::centroid) {

    ScalarDeposition(CentroidDeposition(),
                     lev, mf_to_be_filled, volfrac, flags);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");

  }


}


template <typename F>
void MFIXParticleContainer::
ScalarDeposition(F WeightFunc, int lev,
                 amrex::MultiFab & mf_to_be_filled,
                 const amrex::MultiFab * volfrac,
                 const amrex::FabArray<EBCellFlagFab>* flags)
{
  BL_PROFILE("MFIXParticleContainer::ScalarDeposition()");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      plo = gm.ProbLoArray();
  const auto      dx  = gm.CellSizeArray();
  const auto      dxi = gm.InvCellSizeArray();

  const auto      reg_cell_vol = dx[0]*dx[1]*dx[2];


  using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      const long nrp = pti.numParticles();
      FArrayBox& fab = mf_to_be_filled[pti];

      const Box& bx  = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(bx) != FabType::covered ) {

        const auto& volarr = fab.array();
        const auto& flagsarr = (*flags)[pti].array();
        const auto& vfrac = (*volfrac)[pti].array();

        const amrex::Real deposition_scale_factor =
          mfix::m_deposition_scale_factor;

        amrex::ParallelFor(nrp,
          [pstruct,plo,dx,dxi,vfrac,deposition_scale_factor,volarr,
           reg_cell_vol,WeightFunc,flagsarr]
          AMREX_GPU_DEVICE (int ip) noexcept
          {
            const ParticleType& p = pstruct[ip];

            int i;
            int j;
            int k;

            amrex::Real weights[2][2][2];

            WeightFunc(plo, dx, dxi, flagsarr, p, i, j, k, weights,
                deposition_scale_factor);

            amrex::Real pvol = p.rdata(realData::volume) / reg_cell_vol;

            for (int ii = -1; ii <= 0; ++ii) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int kk = -1; kk <= 0; ++kk) {
                  if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                    continue;

                  amrex::Real weight_vol = weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk);

                  amrex::Gpu::Atomic::Add(&volarr(i+ii,j+jj,k+kk), weight_vol*pvol);
                }
              }
            }
          });
      }
      Gpu::synchronize();
    }
  }
}



void MFIXParticleContainer::
FluidDragForceDeposition(int lev,
                         amrex::MultiFab & mf_tmp_eps,
                         amrex::MultiFab & drag_mf,
                         const amrex::MultiFab * volfrac,
                         const amrex::FabArray<EBCellFlagFab>* flags)
{

  if (mfix::m_deposition_scheme == DepositionScheme::trilinear) {

    FluidDragForceDeposition(TrilinearDeposition(),
                             lev, mf_tmp_eps, drag_mf, volfrac, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::square_dpvm) {

    FluidDragForceDeposition(TrilinearDPVMSquareDeposition(),
                             lev, mf_tmp_eps, drag_mf, volfrac, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::true_dpvm) {

    FluidDragForceDeposition(TrueDPVMDeposition(),
                             lev, mf_tmp_eps, drag_mf, volfrac, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::centroid) {

    FluidDragForceDeposition(CentroidDeposition(),
                             lev, mf_tmp_eps, drag_mf, volfrac, flags);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");

  }

}



template <typename F>
void MFIXParticleContainer::
FluidDragForceDeposition(F WeightFunc, int lev,
                         amrex::MultiFab & mf_tmp_eps,
                         amrex::MultiFab & drag_mf,
                         const amrex::MultiFab * volfrac,
                         const amrex::FabArray<EBCellFlagFab>* flags)
{
  BL_PROFILE("MFIXParticleContainer::FluidDragForceDeposition()");

  // We always use the coarse dx
  const Geometry& gm          = Geom(0);
  const auto      plo         = gm.ProbLoArray();
  const auto      dx          = gm.CellSizeArray();
  const auto      dxi         = gm.InvCellSizeArray();

  const auto      reg_cell_vol = dx[0]*dx[1]*dx[2];

  using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {

    for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();
      const long nrp = pti.numParticles();

      FArrayBox& eps_fab  = mf_tmp_eps[pti];
      FArrayBox& drag_fab = drag_mf[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered ) {

        const auto& drag_arr = drag_fab.array();
        const auto&   volarr =  eps_fab.array();
        const auto& flagsarr = (*flags)[pti].array();
        const auto&    vfrac = (*volfrac)[pti].array();

        const amrex::Real deposition_scale_factor =
          mfix::m_deposition_scale_factor;

        amrex::ParallelFor(nrp,
          [pstruct,plo,dx,dxi,vfrac,volarr,deposition_scale_factor,
           reg_cell_vol,WeightFunc,flagsarr,drag_arr]
          AMREX_GPU_DEVICE (int ip) noexcept
          {
            const ParticleType& p = pstruct[ip];

            int i;
            int j;
            int k;

            amrex::Real weights[2][2][2];

            WeightFunc(plo, dx, dxi, flagsarr, p, i, j, k, weights,
                       deposition_scale_factor);

            amrex::Real pvol = p.rdata(realData::volume) / reg_cell_vol;

            amrex::Real pbeta = p.rdata(realData::dragx) / reg_cell_vol;
            amrex::Real pvx   = p.rdata(realData::velx) * pbeta;
            amrex::Real pvy   = p.rdata(realData::vely) * pbeta;
            amrex::Real pvz   = p.rdata(realData::velz) * pbeta;

            for (int ii = -1; ii <= 0; ++ii) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int kk = -1; kk <= 0; ++kk) {
                  if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                    continue;

                  amrex::Real weight_vol = weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk);

                  amrex::Gpu::Atomic::Add(&volarr(i+ii,j+jj,k+kk), weight_vol*pvol);

                  amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,0),weight_vol*pvx);
                  amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,1),weight_vol*pvy);
                  amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,2),weight_vol*pvz);
                  amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,3),weight_vol*pbeta);
                }
              }
            }
          });
      }
      Gpu::synchronize();
    }
  }
}
