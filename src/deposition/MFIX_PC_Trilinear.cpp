#include <AMReX.H>
#include "AMReX_Particles.H"
#include <MFIXParticleContainer.H>

using namespace amrex;
using namespace std;

void MFIXParticleContainer::
TrilinearDepositionScalar(int lev,
                          amrex::MultiFab & mf_to_be_filled,
                          const amrex::MultiFab * volfrac,
                          const amrex::FabArray<EBCellFlagFab>* flags,
                          int fortran_particle_comp )
{
  BL_PROFILE("MFIXParticleContainer::TrilinearDepositionScalar()");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      plo = gm.ProbLoArray();
  const auto      dx  = gm.CellSizeArray();
  const auto      dxi = gm.InvCellSizeArray();

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

        auto volarr = fab.array();
        auto flagsarr = (*flags)[pti].array();
        auto vfrac = (*volfrac)[pti].array();

        AMREX_FOR_1D ( nrp, ip, {

            const ParticleType& p = pstruct[ip];

            amrex::Real reg_cell_vol = dx[0]*dx[1]*dx[2];

            amrex::Real x = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
            amrex::Real y = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
            amrex::Real z = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

            int i = std::floor(x);
            int j = std::floor(y);
            int k = std::floor(z);

            amrex::Real wx_hi = x - i;
            amrex::Real wy_hi = y - j;
            amrex::Real wz_hi = z - k;

            amrex::Real wx_lo = 1.0 - wx_hi;
            amrex::Real wy_lo = 1.0 - wy_hi;
            amrex::Real wz_lo = 1.0 - wz_hi;

            amrex::Real weights[2][2][2];

            weights[0][0][0] = wx_lo * wy_lo * wz_lo;
            weights[0][0][1] = wx_lo * wy_lo * wz_hi;
            weights[0][1][0] = wx_lo * wy_hi * wz_lo;
            weights[0][1][1] = wx_lo * wy_hi * wz_hi;
            weights[1][0][0] = wx_hi * wy_lo * wz_lo;
            weights[1][0][1] = wx_hi * wy_lo * wz_hi;
            weights[1][1][0] = wx_hi * wy_hi * wz_lo;
            weights[1][1][1] = wx_hi * wy_hi * wz_hi;

            amrex::Real total_weight = 0.0;
            for (int ii = 0; ii <= 1; ++ii)
                for (int jj = 0; jj <= 1; ++jj)
                    for (int kk = 0; kk <= 1; ++kk)
                        total_weight += weights[ii][jj][kk] * vfrac(i-1+ii,j-1+jj,k-1+kk);

            for (int ii = 0; ii <= 1; ++ii)
                for (int jj = 0; jj <= 1; ++jj)
                    for (int kk = 0; kk <= 1; ++kk)
                        weights[ii][jj][kk] /= total_weight;

            amrex::Real pvol = p.rdata(realData::volume) / reg_cell_vol;

            for (int ii = -1; ii <= 0; ++ii) {
                for (int jj = -1; jj <= 0; ++jj) {
                    for (int kk = -1; kk <= 0; ++kk) {
                        if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                            continue;
                        amrex::Gpu::Atomic::Add(&volarr(i+ii,j+jj,k+kk),
                                                weights[ii+1][jj+1][kk+1]*pvol);
                    }
                }
            }
          });

        Gpu::synchronize();
      }
    }
  }
}


void MFIXParticleContainer::
TrilinearDepositionFluidDragForce(int lev,
                                  amrex::MultiFab & drag_mf,
                                  const amrex::MultiFab * volfrac,
                                  const amrex::FabArray<EBCellFlagFab>* flags,
                                  int fortran_beta_comp, int fortran_vel_comp)
{
  BL_PROFILE("MFIXParticleContainer::TrilinearDepositionFluidDragForce()");

  // We always use the coarse dx
  const Geometry& gm          = Geom(0);
  const auto      plo         = gm.ProbLoArray();
  const auto      dx          = gm.CellSizeArray();
  const auto      dxi         = gm.InvCellSizeArray();

  using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {

    for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();
      const long nrp = pti.numParticles();

      FArrayBox& drag_fab = drag_mf[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered ) {

        auto drag_arr = drag_fab.array();
        auto flagsarr = (*flags)[pti].array();
        auto vfrac = (*volfrac)[pti].array();

        AMREX_FOR_1D ( nrp, ip, {

            const ParticleType& p = pstruct[ip];

            amrex::Real reg_cell_vol = dx[0]*dx[1]*dx[2];

            amrex::Real x = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
            amrex::Real y = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
            amrex::Real z = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

            int i = std::floor(x);
            int j = std::floor(y);
            int k = std::floor(z);

            amrex::Real wx_hi = x - i;
            amrex::Real wy_hi = y - j;
            amrex::Real wz_hi = z - k;

            amrex::Real wx_lo = 1.0 - wx_hi;
            amrex::Real wy_lo = 1.0 - wy_hi;
            amrex::Real wz_lo = 1.0 - wz_hi;

            amrex::Real weights[2][2][2];

            weights[0][0][0] = wx_lo * wy_lo * wz_lo;
            weights[0][0][1] = wx_lo * wy_lo * wz_hi;
            weights[0][1][0] = wx_lo * wy_hi * wz_lo;
            weights[0][1][1] = wx_lo * wy_hi * wz_hi;
            weights[1][0][0] = wx_hi * wy_lo * wz_lo;
            weights[1][0][1] = wx_hi * wy_lo * wz_hi;
            weights[1][1][0] = wx_hi * wy_hi * wz_lo;
            weights[1][1][1] = wx_hi * wy_hi * wz_hi;

            amrex::Real total_weight = 0.0;
            for (int ii = 0; ii <= 1; ++ii)
              for (int jj = 0; jj <= 1; ++jj)
                for (int kk = 0; kk <= 1; ++kk)
                  total_weight += weights[ii][jj][kk] * vfrac(i-1+ii,j-1+jj,k-1+kk);

            for (int ii = 0; ii <= 1; ++ii)
              for (int jj = 0; jj <= 1; ++jj)
                for (int kk = 0; kk <= 1; ++kk)
                  weights[ii][jj][kk] /= total_weight;

            amrex::Real pbeta = p.rdata(realData::dragx) / reg_cell_vol;
            amrex::Real pvx   = p.rdata(realData::velx) * pbeta;
            amrex::Real pvy   = p.rdata(realData::vely) * pbeta;
            amrex::Real pvz   = p.rdata(realData::velz) * pbeta;

            for (int ii = -1; ii <= 0; ++ii) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int kk = -1; kk <= 0; ++kk) {
                  if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                    continue;

                  amrex::Real weight_vol = weights[ii+1][jj+1][kk+1];

                  amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,0),weight_vol*pvx);
                  amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,1),weight_vol*pvy);
                  amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,2),weight_vol*pvz);
                  amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,3),weight_vol*pbeta);
                }
              }
            }
          });

        Gpu::synchronize();

      }
    }
  }
}
