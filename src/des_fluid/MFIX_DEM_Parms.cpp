#include <AMReX.H>
#include <MFIX_DEM_Parms.H>

namespace DEM
{

    AMREX_GPU_DEVICE_MANAGED COLLISIONMODEL CollisionModel = LSD;

    AMREX_GPU_DEVICE_MANAGED amrex::Real kt;
    AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w;

    AMREX_GPU_DEVICE_MANAGED amrex::Real kn;
    AMREX_GPU_DEVICE_MANAGED amrex::Real kn_w;

    AMREX_GPU_DEVICE_MANAGED amrex::Real kt_fac = 2.0/7.0;
    AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w_fac = 2.0/7.0;

    void Initialize ()
    {
        if      (CollisionModel == LSD     ) InitializeLSD();
        else if (CollisionModel == HERTZIAN) InitializeHertzian();
        else amrex::Abort("DEM collision model not recognized");
    }    

    void InitializeLSD () {}

    void InitializeHertzian () {}
}
