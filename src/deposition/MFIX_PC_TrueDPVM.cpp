#include <AMReX.H>
#include "AMReX_Particles.H"
#include <MFIXParticleContainer.H>

using namespace amrex;
using namespace std;

void MFIXParticleContainer::
TrueDPVMDepositionScalar(int lev,
                          amrex::MultiFab & mf_to_be_filled,
                          const amrex::MultiFab * volfrac,
                          const amrex::FabArray<EBCellFlagFab>* flags,
                          int fortran_particle_comp )
{
  BL_PROFILE("MFIXParticleContainer::TrueDPVMDepositionScalar()");

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

            amrex::Real wx[2];
            amrex::Real wy[2];
            amrex::Real wz[2];

            amrex::Real weights[2][2][2];

            // We now know which side of the cell (low vs high) the particle sits.
            //

            amrex::Real hh;
            Real rp = p.rdata(realData::radius);

            const int ip = std::floor(x-0.5);
            const int jp = std::floor(y-0.5);
            const int kp = std::floor(z-0.5);

            amrex::Real ratio;
            amrex::Real overlap[3];

            if ( i == ip ) {

              // particle sits on low side of i index
              overlap[0] = amrex::max(0.0, (i*dx[0]) - p.pos(0) + rp);

              wx[0] = (overlap[0] == 0.0) ? 0.0 : fface(overlap[0], rp);
              wx[1] = 1.0 - wx[0];

            } else {

              // particle sits on high side of i-1 index
              overlap[0] = amrex::max(0.0, p.pos(0) + rp - (i*dx[0]));

              wx[1] = (overlap[0] == 0.0) ? 0.0 : fface(overlap[0], rp);
              wx[0] = 1.0 - wx[1];

            }

            if ( j == jp ) {

              // particle sits on low side of j index
              overlap[1] = amrex::max(0.0, (j*dx[1]) - p.pos(1) + rp);

              wy[0] = (overlap[1] == 0.0) ? 0.0 : fface(overlap[1], rp);
              wy[1] = 1.0 - wy[0];

            } else {

              // particle sits on high side of j-1 index
              overlap[1] = amrex::max(0.0, p.pos(1) + rp - (j*dx[1]));

              wy[1] = (overlap[1] == 0.0) ? 0.0 : fface(overlap[1], rp);
              wy[0] = 1.0 - wy[1];

            }


            if ( k == kp ) {

              // particle sits on low side of k index
              overlap[2] = amrex::max(0.0, (k*dx[2]) - p.pos(2) + rp);

              wz[0] = (overlap[2] == 0.0) ? 0.0 : fface(overlap[2], rp);
              wz[1] = 1.0 - wz[0];

            } else {

              // particle sits on high side of k-1 index
              overlap[2] = amrex::max(0.0, p.pos(2) + rp - (k*dx[2]));

              wz[1] = (overlap[2] == 0.0) ? 0.0 : fface(overlap[2], rp);
              wz[0] = 1.0 - wz[1];

            }


            int cuts = 0;
            for ( int ii=0; ii<=2; ++ii) {
              if( overlap[ii] != 0. ) cuts += 1;
            }


            const int i0 = i-ip;
            const int i1 = std::abs(1-i0);
            const int j0 = j-jp;
            const int j1 = std::abs(1-j0);
            const int k0 = k-kp;
            const int k1 = std::abs(1-k0);


            // if ( cuts == 0 || cuts == 1 ) we don't need to do anything as the
            // calculation of the initial 'weights' makes that all work out.

            if ( cuts == 2 ) {

              // This has two sub-cases that have to be addressed. The first is when
              // we only intersect the two faces and need to deposit into three cells;
              // the second is when we intersect the corner and need to deposit
              // into four cells.

              if ( overlap[0] == 0.0 ) {

                weights[i0][0][0] = 0.0;
                weights[i0][0][1] = 0.0;
                weights[i0][1][0] = 0.0;
                weights[i0][1][1] = 0.0;

                amrex::Real ahat = 1.0 - overlap[1] / rp;
                amrex::Real bhat = 1.0 - overlap[2] / rp;

                weights[i1][j0][k0] =  fedge(ahat, bhat);
                weights[i1][j1][k0] =  wy[j0] - weights[i1][j0][k0];
                weights[i1][j0][k1] =  wz[k0] - weights[i1][j0][k0];
                weights[i1][j1][k1] =  1.0 - weights[i1][j0][k0] - weights[i1][j0][k1] - weights[i1][j1][k0];

              } else if ( overlap[1] == 0.0 ) {

                weights[0][j0][0] = 0.0;
                weights[0][j0][1] = 0.0;
                weights[1][j0][0] = 0.0;
                weights[1][j0][1] = 0.0;

                amrex::Real ahat = 1.0 - overlap[0] / rp;
                amrex::Real bhat = 1.0 - overlap[2] / rp;

                weights[i0][j1][k0] =  fedge(ahat, bhat);
                weights[i1][j1][k0] =  wx[i0] - weights[i0][j1][k0];
                weights[i0][j1][k1] =  wz[k0] - weights[i0][j1][k0];
                weights[i1][j1][k1] =  1.0   - weights[i0][j1][k0] - weights[i1][j1][k0] - weights[i0][j1][k1];

              } else {

                weights[0][0][k0] = 0.0;
                weights[0][1][k0] = 0.0;
                weights[1][0][k0] = 0.0;
                weights[1][1][k0] = 0.0;

                amrex::Real ahat = 1.0 - overlap[0] / rp;
                amrex::Real bhat = 1.0 - overlap[1] / rp;

                weights[i0][j0][k1] =  fedge(ahat, bhat);
                weights[i1][j0][k1] =  wx[i0] - weights[i0][j0][k1];
                weights[i0][j1][k1] =  wy[j0] - weights[i0][j0][k1];
                weights[i1][j1][k1] =  1.0   - weights[i0][j0][k1] - weights[i1][j0][k1] - weights[i0][j1][k1];

              }

            } else if ( cuts == 3 ) {

              amrex::Real ahat = 1.0 - overlap[0] / rp;
              amrex::Real bhat = 1.0 - overlap[1] / rp;
              amrex::Real chat = 1.0 - overlap[2] / rp;

              amrex::Real edge12 = fedge(ahat, bhat);
              amrex::Real edge13 = fedge(ahat, chat);
              amrex::Real edge23 = fedge(bhat, chat);

              amrex::Real corner = fcorner(ahat, bhat, chat, edge12);

              // corner
              weights[i0][j0][k0] =  corner;

              // face overlaps
              weights[i0][j1][k1] =  wx[i0] - (edge12 + edge13) + corner;
              weights[i1][j0][k1] =  wy[j0] - (edge12 + edge23) + corner;
              weights[i1][j1][k0] =  wz[k0] - (edge13 + edge23) + corner;

              // edge overlaps
              weights[i0][j0][k1] =  edge12 - corner;
              weights[i0][j1][k0] =  edge13 - corner;
              weights[i1][j0][k0] =  edge23 - corner;

              // centroid
              weights[i1][j1][k1] =  1.0 - (wx[i0] + wy[j0] + wz[k0] -
                                            edge12 - edge13 - edge23 + corner);

            } else {

              weights[0][0][0] = wx[0] * wy[0] * wz[0];
              weights[0][0][1] = wx[0] * wy[0] * wz[1];
              weights[0][1][0] = wx[0] * wy[1] * wz[0];
              weights[0][1][1] = wx[0] * wy[1] * wz[1];
              weights[1][0][0] = wx[1] * wy[0] * wz[0];
              weights[1][0][1] = wx[1] * wy[0] * wz[1];
              weights[1][1][0] = wx[1] * wy[1] * wz[0];
              weights[1][1][1] = wx[1] * wy[1] * wz[1];

            }


            // amrex::Real pvol = p.rdata(realData::volume) / reg_cell_vol;
            amrex::Real pvol = p.rdata(realData::volume);

            for (int ii = -1; ii <= 0; ++ii) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int kk = -1; kk <= 0; ++kk) {
                  if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                    continue;
                  amrex::Real eps = pvol / (reg_cell_vol * vfrac(i+ii,j+jj,k+kk));

                  amrex::Gpu::Atomic::Add(&volarr(i+ii,j+jj,k+kk),
                                          weights[ii+1][jj+1][kk+1]*eps);
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
TrueDPVMDepositionFluidDragForce(int lev,
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
