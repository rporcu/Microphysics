#ifndef MFIX_DEM_PARMS_H_
#define MFIX_DEM_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace DEM
{
    // tangential component of spring constant
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w;

    // normal component of spring constant
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kn;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kn_w;
    
    // assumed to be a constant factor relating these for linear spring-dashpot
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt_fac;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w_fac;

    void Initialize ();
}

#endif
