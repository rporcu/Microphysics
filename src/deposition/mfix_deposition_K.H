#include "MFIXParticleContainer.H"


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fface (const amrex::Real& overlap, const amrex::Real& rad)
{
  amrex::Real del = overlap / rad;
  return 0.25*( 3.0 - del)*del*del;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fedge (const amrex::Real& ahat,
                   const amrex::Real& bhat)
{

  amrex::Real xhat_sq = 1.0 - ahat*ahat - bhat*bhat;
  if ( xhat_sq <= 0.0 ) {
    return 0.0;
  } else {

    amrex::Real xhat = sqrt(xhat_sq);

    return (1./(4.*M_PI)) * ( 2.0*ahat*bhat*xhat -
                            ( 3.*ahat - ahat*ahat*ahat) * std::atan(xhat/bhat) -
                            ( 3.*bhat - bhat*bhat*bhat) * std::atan(xhat/ahat) +
                            ( 2. * std::atan(xhat*ahat/bhat)) +
                            ( 2. * std::atan(xhat*bhat/ahat)));
  }

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fcorner (const amrex::Real& ahat,
                     const amrex::Real& bhat,
                     const amrex::Real& chat,
                     const amrex::Real& fe12)
{

  if ( ahat*ahat + bhat*bhat + chat*chat >= 1.0 ) {
    return 0.0;

  } else {

    amrex::Real Ahat = sqrt(1.0 - ahat*ahat - chat*chat);
    amrex::Real Bhat = sqrt(1.0 - bhat*bhat - chat*chat);

    return 0.5*fe12 - (1./(8.*M_PI)) *
      (6.0*ahat*bhat*chat - 2.*ahat*Ahat*chat - 2.*bhat*Bhat*chat -
       (3.*bhat - bhat*bhat*bhat)*std::atan(chat/Bhat) -
       (3.*ahat - ahat*ahat*ahat)*std::atan(chat/Ahat) +
       (3.*chat - chat*chat*chat)*(std::atan(Ahat/ahat) - std::atan(bhat/Bhat)) +
       2.*(std::atan(chat*(ahat/Ahat)) + std::atan(chat*(bhat/Bhat))));

  }

}

struct TrilinearDeposition
{
  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  void operator() (amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                   amrex::Array4<amrex::Real const> const& vfrac,
                   const MFIXParticleContainer::ParticleType& p,
                   int& i, int& j, int& k, amrex::Real weights[2][2][2])
  {

    amrex::Real x = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real y = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
    amrex::Real z = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

    i = std::floor(x);
    j = std::floor(y);
    k = std::floor(z);

    amrex::Real wx[2];
    amrex::Real wy[2];
    amrex::Real wz[2];

    wx[1] = x - i;
    wy[1] = y - j;
    wz[1] = z - k;

    wx[0] = 1.0 - wx[1];
    wy[0] = 1.0 - wy[1];
    wz[0] = 1.0 - wz[1];

    amrex::Real total_weight = 0.0;

    for (int ii = 0; ii <= 1; ++ii)
      for (int jj = 0; jj <= 1; ++jj)
        for (int kk = 0; kk <= 1; ++kk){
          weights[ii][jj][kk] = wx[ii]*wy[jj]*wz[kk];
          total_weight += weights[ii][jj][kk] * vfrac(i-1+ii,j-1+jj,k-1+kk);
        }

    for (int ii = 0; ii <= 1; ++ii)
      for (int jj = 0; jj <= 1; ++jj)
        for (int kk = 0; kk <= 1; ++kk)
          weights[ii][jj][kk] /= total_weight;

  }

};
