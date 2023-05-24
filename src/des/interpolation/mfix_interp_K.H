#ifndef _MFIX_INTERP_K_H_
#define _MFIX_INTERP_K_H_

#include <mfix_pc.H>

/****************************************************************
 *                                                              *
 *                                                              *
 *                                                              *
 *                                                              *
 ***************************************************************/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void trilinear_interp (const amrex::RealVect& pos,
                       amrex::Real* interp_vals,
                       amrex::Array4<amrex::Real const> const& interp_array,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
                       const int interp_comp)
{
    const amrex::RealVect lx((pos[0] - plo[0])*dxi[0] + 0.5,
                             (pos[1] - plo[1])*dxi[1] + 0.5,
                             (pos[2] - plo[2])*dxi[2] + 0.5);

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


/****************************************************************
 *                                                              *
 *                                                              *
 *                                                              *
 *                                                              *
 ***************************************************************/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void trilinear_interp (const amrex::RealVect& pos,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
                       amrex::Array4<amrex::Real const> const& interp_array1,
                       amrex::Real* interp_vals1,
                       const int interp_comp1,
                       amrex::Array4<amrex::Real const> const& interp_array2,
                       amrex::Real* interp_vals2,
                       const int interp_comp2,
                       const int np,
                       const int p_id2)
{
    const amrex::RealVect lx((pos[0] - plo[0])*dxi[0] + 0.5,
                             (pos[1] - plo[1])*dxi[1] + 0.5,
                             (pos[2] - plo[2])*dxi[2] + 0.5);

    const amrex::IntVect ijk = lx.floor();

    int i = ijk[0]; int j = ijk[1]; int k = ijk[2];

    // Weights
    const amrex::RealVect sx_hi = lx - ijk;
    const amrex::RealVect sx_lo = 1 - sx_hi;

    amrex::GpuArray<amrex::Real,8> weights;
    weights[0] = sx_lo[0]*sx_lo[1]*sx_lo[2];
    weights[1] = sx_lo[0]*sx_lo[1]*sx_hi[2];
    weights[2] = sx_lo[0]*sx_hi[1]*sx_lo[2];
    weights[3] = sx_lo[0]*sx_hi[1]*sx_hi[2];
    weights[4] = sx_hi[0]*sx_lo[1]*sx_lo[2];
    weights[5] = sx_hi[0]*sx_lo[1]*sx_hi[2];
    weights[6] = sx_hi[0]*sx_hi[1]*sx_lo[2];
    weights[7] = sx_hi[0]*sx_hi[1]*sx_hi[2];

    for (int n = 0; n < interp_comp1; n++)
      interp_vals1[n] = weights[0]*interp_array1(i-1, j-1, k-1,n) +
                        weights[1]*interp_array1(i-1, j-1, k  ,n) +
                        weights[2]*interp_array1(i-1, j  , k-1,n) +
                        weights[3]*interp_array1(i-1, j  , k  ,n) +
                        weights[4]*interp_array1(i  , j-1, k-1,n) +
                        weights[5]*interp_array1(i  , j-1, k  ,n) +
                        weights[6]*interp_array1(i  , j  , k-1,n) +
                        weights[7]*interp_array1(i  , j  , k  ,n);

    for (int n = 0; n < interp_comp2; n++)
      interp_vals2[n*np+p_id2] = weights[0]*interp_array2(i-1, j-1, k-1,n) +
                                 weights[1]*interp_array2(i-1, j-1, k  ,n) +
                                 weights[2]*interp_array2(i-1, j  , k-1,n) +
                                 weights[3]*interp_array2(i-1, j  , k  ,n) +
                                 weights[4]*interp_array2(i  , j-1, k-1,n) +
                                 weights[5]*interp_array2(i  , j-1, k  ,n) +
                                 weights[6]*interp_array2(i  , j  , k-1,n) +
                                 weights[7]*interp_array2(i  , j  , k  ,n);
}

#endif
