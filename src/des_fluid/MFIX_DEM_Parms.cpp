#include <AMReX.H>
#include <AMReX_Arena.H>
#include <MFIX_DEM_Parms.H>

#include <mfix_des_F.H>

namespace DEMParams
{
    AMREX_GPU_DEVICE_MANAGED COLLISIONMODEL CollisionModel = LSD;
    AMREX_GPU_DEVICE_MANAGED int NPHASE;
    
    AMREX_GPU_DEVICE_MANAGED amrex::Real kt;
    AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w;
    
    AMREX_GPU_DEVICE_MANAGED amrex::Real kn;
    AMREX_GPU_DEVICE_MANAGED amrex::Real kn_w;
    
    AMREX_GPU_DEVICE_MANAGED amrex::Real kt_fac = 2.0/7.0;
    AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w_fac = 2.0/7.0;
    
    // normal and tangential components of the damping coefficients
    AMREX_GPU_DEVICE_MANAGED amrex::Real etan[NMAX][NMAX];
    AMREX_GPU_DEVICE_MANAGED amrex::Real etan_w[NMAX];
    
    AMREX_GPU_DEVICE_MANAGED amrex::Real etat[NMAX][NMAX];
    AMREX_GPU_DEVICE_MANAGED amrex::Real etat_w[NMAX];
    
    // coefficients of restitution, normal and tangential
    AMREX_GPU_DEVICE_MANAGED amrex::Real en_input[NMAX+NMAX*(NMAX-1)/2];
    AMREX_GPU_DEVICE_MANAGED amrex::Real et_input[NMAX+NMAX*(NMAX-1)/2];
    
    AMREX_GPU_DEVICE_MANAGED amrex::Real en_w_input[NMAX];
    AMREX_GPU_DEVICE_MANAGED amrex::Real et_w_input[NMAX];
    
    void Initialize ()
    {
        if      (CollisionModel == LSD     ) InitializeLSD();
        else if (CollisionModel == HERTZIAN) InitializeHertzian();
        else amrex::Abort("DEM collision model not recognized");
    }
    
    void InitializeLSD () 
    {
        get_collision_coefficients(&NPHASE, &kt, &kt_w, &kn, &kn_w, 
                                   &etan[0][0], &etan_w[0], &etat[0][0], &etat_w[0]);
    }
    
    void InitializeHertzian ()
    {
        amrex::Abort("The Hertzian collision model is currently not recognized");
    }
}
