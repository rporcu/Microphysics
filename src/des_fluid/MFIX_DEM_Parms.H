#ifndef MFIX_DEM_PARMS_H_
#define MFIX_DEM_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace DEMParams
{
    constexpr int NMAX = 10;
    extern AMREX_GPU_DEVICE_MANAGED int NPHASE;
    
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

    // friction coefficients for particle-particle and particle-wall collisions
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real mew;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real mew_w;
    
    // assumed to be a constant factor relating these for linear spring-dashpot
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt_fac;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w_fac;

    // normal and tangential components of the damping coefficients
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real etan[NMAX][NMAX];
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real etan_w[NMAX];

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real etat[NMAX][NMAX];
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real etat_w[NMAX];
    
    // coefficients of restitution, normal and tangential
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real en_input[NMAX+NMAX*(NMAX-1)/2];
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real et_input[NMAX+NMAX*(NMAX-1)/2];
    
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real en_w_input[NMAX];
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real et_w_input[NMAX];
    
    void Initialize ();
    
    void InitializeLSD ();
    
    void InitializeHertzian ();
}

#endif
