#ifndef MFIX_DEM_PARMS_H_
#define MFIX_DEM_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace DEM
{

    enum COLLISIONMODEL {
        HERTZIAN = 0,
        LSD
    };

    extern AMREX_GPU_DEVICE_MANAGED COLLISIONMODEL CollisionModel;

    // tangential component of spring constant
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w;

    // normal component of spring constant
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kn;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kn_w;
    
    // assumed to be a constant factor relating these for linear spring-dashpot
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt_fac;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w_fac;

    // normal and tangential components of the damping coefficients
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real etan;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real etan_w;

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real etat;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real etat_w;

    // coefficients of restitution, normal and tangential
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real en_input[1];
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real et_input[1];

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real en_w_input[1];
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real et_w_input[1];
       
    void Initialize ();

    void InitializeLSD ();

    void InitializeHertzian ();    
}

#endif
