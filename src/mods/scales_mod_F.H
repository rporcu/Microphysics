#ifndef __SCALES_MOD_F_H_
#define __SCALES_MOD_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

amrex::Real get_p_ref ();

amrex::Real get_p_scale ();

#ifdef __cplusplus
}
#endif

AMREX_GPU_HOST_DEVICE
amrex::Real scale_pressure_cpp (const amrex::Real XXX,
                                const amrex::Real P_ref,
                                const amrex::Real P_scale);

AMREX_GPU_HOST_DEVICE
amrex::Real unscale_pressure_cpp (const amrex::Real XXX,
                                  const amrex::Real P_ref,
                                  const amrex::Real P_scale);

#endif
