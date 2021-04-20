#ifndef _MFIX_INTERP_K_H_
#define _MFIX_INTERP_K_H_

#include "mfix_pc.H"

/****************************************************************
 *                                                              *
 *                                                              *
 *                                                              *
 *                                                              *
 ***************************************************************/
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void trilinear_interp (const amrex::RealVect& pos,
                       amrex::Real* interp_vals,
                       amrex::Array4<amrex::Real const> const& interp_array,
                       const amrex::RealVect& plo,
                       const amrex::RealVect& dxi,
                       const int interp_comp)
{

    const amrex::RealVect lx = (pos - plo)*dxi + 0.5;
    const amrex::IntVect ijk = lx.floor();

    int i = ijk[0]; int j = ijk[1]; int k = ijk[2];

    // Weights
    const amrex::RealVect sx_hi = lx - ijk;
    const amrex::RealVect sx_lo = 1 - sx_hi;

    for (int n = 0; n < interp_comp; n++)
       interp_vals[n] = sx_lo[0]*sx_lo[1]*sx_lo[2]*interp_array(i-1, j-1, k-1,n) +
                        sx_lo[0]*sx_lo[1]*sx_hi[2]*interp_array(i-1, j-1, k  ,n) +
                        sx_lo[0]*sx_hi[1]*sx_lo[2]*interp_array(i-1, j  , k-1,n) +
                        sx_lo[0]*sx_hi[1]*sx_hi[2]*interp_array(i-1, j  , k  ,n) +
                        sx_hi[0]*sx_lo[1]*sx_lo[2]*interp_array(i  , j-1, k-1,n) +
                        sx_hi[0]*sx_lo[1]*sx_hi[2]*interp_array(i  , j-1, k  ,n) +
                        sx_hi[0]*sx_hi[1]*sx_lo[2]*interp_array(i  , j  , k-1,n) +
                        sx_hi[0]*sx_hi[1]*sx_hi[2]*interp_array(i  , j  , k  ,n);
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void trilinear_interp (const amrex::RealVect& pos,
                       int& ip,
                       int& jp,
                       int& kp,
                       amrex::GpuArray<amrex::GpuArray<amrex::GpuArray<Real,2>,2>,2>& weights,
                       amrex::Real* interp_vals,
                       amrex::Array4<amrex::Real const> const& interp_array,
                       const amrex::RealVect& plo,
                       const amrex::RealVect& dxi,
                       const int interp_comp)
{
    const amrex::RealVect lx = (pos - plo)*dxi + 0.5;
    const amrex::IntVect ijk = lx.floor();

    ip = ijk[0]; jp = ijk[1]; kp = ijk[2];

    // Weights
    const amrex::RealVect sx_hi = lx - ijk;
    const amrex::RealVect sx_lo = 1 - sx_hi;

    weights[0][0][0] = sx_lo[0]*sx_lo[1]*sx_lo[2];
    weights[0][0][1] = sx_lo[0]*sx_lo[1]*sx_hi[2];
    weights[0][1][0] = sx_lo[0]*sx_hi[1]*sx_lo[2];
    weights[0][1][1] = sx_lo[0]*sx_hi[1]*sx_hi[2];
    weights[1][0][0] = sx_hi[0]*sx_lo[1]*sx_lo[2];
    weights[1][0][1] = sx_hi[0]*sx_lo[1]*sx_hi[2];
    weights[1][1][0] = sx_hi[0]*sx_hi[1]*sx_lo[2];
    weights[1][1][1] = sx_hi[0]*sx_hi[1]*sx_hi[2];

    for (int nn(0); nn < interp_comp; ++nn)
    for (int kk(0); kk < 2; ++kk)
    for (int jj(0); jj < 2; ++jj)
    for (int ii(0); ii < 2; ++ii)
      interp_vals[nn] += weights[ii][jj][kk]*interp_array(ip+ii-1,jp+jj-1,kp+kk-1,nn);
}

#endif
