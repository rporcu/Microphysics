#include <AMReX.H>
#include "AMReX_Particles.H"
#include <mfix_pc.H>

#include <mfix_deposition_K.H>
#include <mfix.H>

void MFIXParticleContainer::
MFIX_PC_SolidsVelocityDeposition (int lev,
                                  amrex::MultiFab & vel_s_mf,
                                  const amrex::FabArray<EBCellFlagFab>* flags)
{
  BL_PROFILE("(MFIXParticleContainer::MFIX_PC_SolidsVelocityDeposition)");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      plo = gm.ProbLoArray();
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

      FArrayBox& vel_s_fab = vel_s_mf[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered ) {

        const auto& vel_s_arr = vel_s_fab.array();

        amrex::ParallelFor(nrp,
          [pstruct,plo,dxi,vel_s_arr]
          AMREX_GPU_DEVICE (int ip) noexcept
          {
            const ParticleType& p = pstruct[ip];

            const amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
            const amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
            const amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

            const int i = static_cast<int>(amrex::Math::floor(lx));
            const int j = static_cast<int>(amrex::Math::floor(ly));
            const int k = static_cast<int>(amrex::Math::floor(lz));

            const amrex::Real wx_hi(lx - static_cast<amrex::Real>(i));
            const amrex::Real wy_hi(ly - static_cast<amrex::Real>(j));
            const amrex::Real wz_hi(lz - static_cast<amrex::Real>(k));

            const amrex::Real wx_lo(1.0 - wx_hi);
            const amrex::Real wy_lo(1.0 - wy_hi);
            const amrex::Real wz_lo(1.0 - wz_hi);

            const amrex::Real pmass = p.rdata(realData::statwt)*p.rdata(realData::mass);


            {// Deposition of x velocity -- x-face deposition

              const amrex::Real pvelx = pmass*p.rdata(realData::velx);

              const amrex::Real lxc = (p.pos(0) - plo[0]) * dxi[0];

              const int ii = static_cast<int>(amrex::Math::floor(lxc));

              const amrex::Real wu_hi(lxc - static_cast<amrex::Real>(ii));
              const amrex::Real wu_lo(1.0 - wu_hi);

              amrex::Gpu::Atomic::Add(&vel_s_arr(ii,   j-1, k-1, 0),wu_lo*wy_lo*wz_lo*pvelx);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii,   j-1, k  , 0),wu_lo*wy_lo*wz_hi*pvelx);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii,   j,   k-1, 0),wu_lo*wy_hi*wz_lo*pvelx);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii,   j,   k  , 0),wu_lo*wy_hi*wz_hi*pvelx);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii+1, j-1, k-1, 0),wu_hi*wy_lo*wz_lo*pvelx);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii+1, j-1, k  , 0),wu_hi*wy_lo*wz_hi*pvelx);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii+1, j,   k-1, 0),wu_hi*wy_hi*wz_lo*pvelx);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii+1, j,   k  , 0),wu_hi*wy_hi*wz_hi*pvelx);

              amrex::Gpu::Atomic::Add(&vel_s_arr(ii,   j-1, k-1, 3),wu_lo*wy_lo*wz_lo*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii,   j-1, k  , 3),wu_lo*wy_lo*wz_hi*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii,   j,   k-1, 3),wu_lo*wy_hi*wz_lo*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii,   j,   k  , 3),wu_lo*wy_hi*wz_hi*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii+1, j-1, k-1, 3),wu_hi*wy_lo*wz_lo*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii+1, j-1, k  , 3),wu_hi*wy_lo*wz_hi*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii+1, j,   k-1, 3),wu_hi*wy_hi*wz_lo*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(ii+1, j,   k  , 3),wu_hi*wy_hi*wz_hi*pmass);
            }


            {// Deposition of y velocity -- y-face deposition

              const amrex::Real pvely = pmass*p.rdata(realData::vely);

              const amrex::Real lyc = (p.pos(1) - plo[1]) * dxi[1];

              const int jj = static_cast<int>(amrex::Math::floor(lyc));

              const amrex::Real wv_hi(lyc - static_cast<amrex::Real>(jj));
              const amrex::Real wv_lo(1.0 - wv_hi);

              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, jj,   k-1, 1),wx_lo*wv_lo*wz_lo*pvely);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, jj,   k  , 1),wx_lo*wv_lo*wz_hi*pvely);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, jj+1, k-1, 1),wx_lo*wv_hi*wz_lo*pvely);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, jj+1, k  , 1),wx_lo*wv_hi*wz_hi*pvely);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   jj,   k-1, 1),wx_hi*wv_lo*wz_lo*pvely);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   jj,   k  , 1),wx_hi*wv_lo*wz_hi*pvely);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   jj+1, k-1, 1),wx_hi*wv_hi*wz_lo*pvely);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   jj+1, k  , 1),wx_hi*wv_hi*wz_hi*pvely);

              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, jj,   k-1, 4),wx_lo*wv_lo*wz_lo*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, jj,   k  , 4),wx_lo*wv_lo*wz_hi*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, jj+1, k-1, 4),wx_lo*wv_hi*wz_lo*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, jj+1, k  , 4),wx_lo*wv_hi*wz_hi*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   jj,   k-1, 4),wx_hi*wv_lo*wz_lo*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   jj,   k  , 4),wx_hi*wv_lo*wz_hi*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   jj+1, k-1, 4),wx_hi*wv_hi*wz_lo*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   jj+1, k  , 4),wx_hi*wv_hi*wz_hi*pmass);
            }


            {// Deposition of z velocity -- z-face deposition

              const amrex::Real pvelz = pmass*p.rdata(realData::velz);
              const amrex::Real lzc = (p.pos(2) - plo[2]) * dxi[2];

              const int kk = static_cast<int>(amrex::Math::floor(lzc));

              const amrex::Real ww_hi(lzc - static_cast<amrex::Real>(kk));
              const amrex::Real ww_lo(1.0 - ww_hi);

              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j-1, kk  , 2),wx_lo*wy_lo*ww_lo*pvelz);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j-1, kk+1, 2),wx_lo*wy_lo*ww_hi*pvelz);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j,   kk  , 2),wx_lo*wy_hi*ww_lo*pvelz);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j,   kk+1, 2),wx_lo*wy_hi*ww_hi*pvelz);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j-1, kk  , 2),wx_hi*wy_lo*ww_lo*pvelz);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j-1, kk+1, 2),wx_hi*wy_lo*ww_hi*pvelz);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j,   kk  , 2),wx_hi*wy_hi*ww_lo*pvelz);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j,   kk+1, 2),wx_hi*wy_hi*ww_hi*pvelz);

              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j-1, kk  , 5),wx_lo*wy_lo*ww_lo*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j-1, kk+1, 5),wx_lo*wy_lo*ww_hi*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j,   kk  , 5),wx_lo*wy_hi*ww_lo*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j,   kk+1, 5),wx_lo*wy_hi*ww_hi*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j-1, kk  , 5),wx_hi*wy_lo*ww_lo*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j-1, kk+1, 5),wx_hi*wy_lo*ww_hi*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j,   kk  , 5),wx_hi*wy_hi*ww_lo*pmass);
              amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j,   kk+1, 5),wx_hi*wy_hi*ww_hi*pmass);
            }
          });
      }
    }
  }
}
