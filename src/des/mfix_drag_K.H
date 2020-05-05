#ifndef _MFIX_DRAG_K_H_
#define _MFIX_DRAG_K_H_

#include "MFIXParticleContainer.H"

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void trilinear_interp (const amrex::RealVect& pos, 
                       amrex::Real* velfp,
                       amrex::Array4<amrex::Real const> const& vel,
                       const amrex::RealVect& plo,
                       const amrex::RealVect& dxi)
{
    const amrex::RealVect lx = (pos - plo)*dxi + 0.5;
    const amrex::IntVect ijk = lx.floor();

    int i = ijk[0]; int j = ijk[1]; int k = ijk[2];

    // Weights
    const amrex::RealVect sx_hi = lx - ijk;
    const amrex::RealVect sx_lo = 1 - sx_hi;

    for (int n = 0; n < 3; n++)
       velfp[n] = sx_lo[0]*sx_lo[1]*sx_lo[2]*vel(i-1, j-1, k-1,n) +
                  sx_lo[0]*sx_lo[1]*sx_hi[2]*vel(i-1, j-1, k  ,n) +
                  sx_lo[0]*sx_hi[1]*sx_lo[2]*vel(i-1, j  , k-1,n) +
                  sx_lo[0]*sx_hi[1]*sx_hi[2]*vel(i-1, j  , k  ,n) +
                  sx_hi[0]*sx_lo[1]*sx_lo[2]*vel(i  , j-1, k-1,n) +
                  sx_hi[0]*sx_lo[1]*sx_hi[2]*vel(i  , j-1, k  ,n) +
                  sx_hi[0]*sx_hi[1]*sx_lo[2]*vel(i  , j  , k-1,n) +
                  sx_hi[0]*sx_hi[1]*sx_hi[2]*vel(i  , j  , k  ,n);
}

#endif
