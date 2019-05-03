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

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void normal_from_ls (amrex::Real* normal, amrex::Real* pos,
                     amrex::Array4<amrex::Real const> const& phi,
                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::Real lx = (pos[0] - plo[0])*dxi[0];
    amrex::Real ly = (pos[1] - plo[1])*dxi[1];
    amrex::Real lz = (pos[2] - plo[2])*dxi[2];

    int i = floor(lx);
    int j = floor(ly);
    int k = floor(lz);

    amrex::Real wx_hi = lx - i;
    amrex::Real wy_hi = ly - j;
    amrex::Real wz_hi = lz - k;

    amrex::Real wx_lo = 1.0 - wx_hi;
    amrex::Real wy_lo = 1.0 - wy_hi;
    amrex::Real wz_lo = 1.0 - wz_hi;

    normal[0] = - phi(i,   j,   k  )*dxi[0] * wy_lo * wz_lo
                + phi(i+1, j,   k  )*dxi[0] * wy_lo * wz_lo
                - phi(i,   j+1, k  )*dxi[0] * wy_hi * wz_lo
                + phi(i+1, j+1, k  )*dxi[0] * wy_hi * wz_lo
                - phi(i,   j,   k+1)*dxi[0] * wy_lo * wz_hi
                + phi(i+1, j,   k+1)*dxi[0] * wy_lo * wz_hi
                - phi(i,   j+1, k+1)*dxi[0] * wy_hi * wz_hi
                + phi(i+1, j+1, k+1)*dxi[0] * wy_hi * wz_hi;

    normal[1] = - phi(i,   j,   k  )*dxi[1] * wx_lo * wz_lo
                + phi(i,   j+1, k  )*dxi[1] * wx_lo * wz_lo
                - phi(i+1, j,   k  )*dxi[1] * wx_hi * wz_lo
                + phi(i+1, j+1, k  )*dxi[1] * wx_hi * wz_lo
                - phi(i,   j,   k+1)*dxi[1] * wx_lo * wz_hi
                + phi(i,   j+1, k+1)*dxi[1] * wx_lo * wz_hi
                - phi(i+1, j,   k+1)*dxi[1] * wx_hi * wz_hi
                + phi(i+1, j+1, k+1)*dxi[1] * wx_hi * wz_hi;

    normal[2] = - phi(i,   j,   k  )*dxi[2] * wx_lo * wy_lo
                + phi(i,   j,   k+1)*dxi[2] * wx_lo * wy_lo
                - phi(i+1, j,   k  )*dxi[2] * wx_hi * wy_lo
                + phi(i+1, j,   k+1)*dxi[2] * wx_hi * wy_lo
                - phi(i,   j+1, k  )*dxi[2] * wx_lo * wy_hi
                + phi(i,   j+1, k+1)*dxi[2] * wx_lo * wy_hi
                - phi(i+1, j+1, k  )*dxi[2] * wx_hi * wy_hi
                + phi(i+1, j+1, k+1)*dxi[2] * wx_hi * wy_hi;
}
