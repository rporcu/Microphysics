#include <AMReX.H>
#include "AMReX_Particles.H"
#include <MFIXParticleContainer.H>

#include <mfix_deposition_K.H>
#include <mfix.H>

void MFIXParticleContainer::
MFIX_PC_SolidsVelocityDeposition (int lev,
                                  amrex::MultiFab & vel_s_mf,
                                  const amrex::FabArray<EBCellFlagFab>* flags)
{

  if (mfix::m_deposition_scheme == DepositionScheme::trilinear) {

    MFIX_PC_SolidsVelocityDeposition(TrilinearDeposition(),
                             lev, vel_s_mf, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::square_dpvm) {

    MFIX_PC_SolidsVelocityDeposition(TrilinearDPVMSquareDeposition(),
                             lev, vel_s_mf, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::true_dpvm) {

    MFIX_PC_SolidsVelocityDeposition(TrueDPVMDeposition(),
                             lev, vel_s_mf, flags);

  } else if (mfix::m_deposition_scheme == DepositionScheme::centroid) {

    MFIX_PC_SolidsVelocityDeposition(CentroidDeposition(),
                             lev, vel_s_mf, flags);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");

  }

}


template <typename F>
void MFIXParticleContainer::
MFIX_PC_SolidsVelocityDeposition (F WeightFunc, int lev,
                                  amrex::MultiFab & vel_s_mf,
                                  const amrex::FabArray<EBCellFlagFab>* flags)
{
  BL_PROFILE("(MFIXParticleContainer::MFIX_PC_SolidsVelocityDeposition)");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      plo = gm.ProbLoArray();
  const auto      dx  = gm.CellSizeArray();
  const auto      dxi = gm.InvCellSizeArray();

  const auto      inv_reg_cell_vol = 1.0 / (dx[0]*dx[1]*dx[2]);

  using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {

    for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();
      const long nrp = pti.numParticles();

      FArrayBox& vel_s_fab = vel_s_mf[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered ) {

        const auto& vel_s_arr = vel_s_fab.array();
        const auto& flagsarr = (*flags)[pti].array();

        const amrex::Real deposition_scale_factor =
          mfix::m_deposition_scale_factor;

        amrex::ParallelFor(nrp,
          [pstruct,plo,dx,dxi,deposition_scale_factor,
           inv_reg_cell_vol,WeightFunc,flagsarr,vel_s_arr]
          AMREX_GPU_DEVICE (int ip) noexcept
          {
            const ParticleType& p = pstruct[ip];

            int i;
            int j;
            int k;

            amrex::Real weights[2][2][2];

            WeightFunc(plo, dx, dxi, flagsarr, p, i, j, k, weights,
                       deposition_scale_factor);

            amrex::Real pmass = p.rdata(realData::statwt)*p.rdata(realData::mass);
            amrex::Real pvelx = pmass*p.rdata(realData::velx);
            amrex::Real pvely = pmass*p.rdata(realData::vely);
            amrex::Real pvelz = pmass*p.rdata(realData::velz);

            for (int ii = -1; ii <= 0; ++ii) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int kk = -1; kk <= 0; ++kk) {
                  if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                    continue;

                  amrex::Real weight = weights[ii+1][jj+1][kk+1];

                  amrex::Gpu::Atomic::Add(&vel_s_arr(i+ii,j+jj,k+kk,0),weight*pvelx);
                  amrex::Gpu::Atomic::Add(&vel_s_arr(i+ii,j+jj,k+kk,1),weight*pvely);
                  amrex::Gpu::Atomic::Add(&vel_s_arr(i+ii,j+jj,k+kk,2),weight*pvelz);
                  amrex::Gpu::Atomic::Add(&vel_s_arr(i+ii,j+jj,k+kk,3),weight*pmass);
                }
              }
            }
          });
      }
    }
  }
}
