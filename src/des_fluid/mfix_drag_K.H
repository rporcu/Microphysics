#include "MFIXParticleContainer.H"

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void trilinear_interp (MFIXParticleContainer::ParticleType const& p, 
                       amrex::Real* velfp,
                       amrex::Array4<amrex::Real const> const& vel,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::Real lx = (p.pos(0) - plo[0])*dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1])*dxi[1] + 0.5;
    amrex::Real lz = (p.pos(2) - plo[2])*dxi[2] + 0.5;

    int i = std::floor(lx);
    int j = std::floor(ly);
    int k = std::floor(lz);

    // Weights
    amrex::Real sx_hi = lx - i;  Real sx_lo = 1.0 - sx_hi;
    amrex::Real sy_hi = ly - j;  Real sy_lo = 1.0 - sy_hi;
    amrex::Real sz_hi = lz - k;  Real sz_lo = 1.0 - sz_hi;

    for (int n = 0; n < 3; n++)
       velfp[n] = sx_lo*sy_lo*sz_lo*vel(i-1, j-1, k-1,n) +
                  sx_lo*sy_lo*sz_hi*vel(i-1, j-1, k  ,n) +
                  sx_lo*sy_hi*sz_lo*vel(i-1, j  , k-1,n) +
                  sx_lo*sy_hi*sz_hi*vel(i-1, j  , k  ,n) +
                  sx_hi*sy_lo*sz_lo*vel(i  , j-1, k-1,n) +
                  sx_hi*sy_lo*sz_hi*vel(i  , j-1, k  ,n) +
                  sx_hi*sy_hi*sz_lo*vel(i  , j  , k-1,n) +
                  sx_hi*sy_hi*sz_hi*vel(i  , j  , k  ,n);
}
