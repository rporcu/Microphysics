#include <AMReX.H>
#include <AMReX_Arena.H>
#include <MFIX_DEM_Parms.H>
#include <AMReX_Print.H>

#include <mfix_des_F.H>

namespace DEM
{
    AMREX_GPU_DEVICE_MANAGED COLLISIONMODEL CollisionModel = LSD;
    AMREX_GPU_DEVICE_MANAGED int NPHASE;

    AMREX_GPU_DEVICE_MANAGED amrex::Real kt;
    AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w;

    AMREX_GPU_DEVICE_MANAGED amrex::Real kn;
    AMREX_GPU_DEVICE_MANAGED amrex::Real kn_w;

    AMREX_GPU_DEVICE_MANAGED amrex::Real kt_fac = 2.0/7.0;
    AMREX_GPU_DEVICE_MANAGED amrex::Real kt_w_fac = 2.0/7.0;

    AMREX_GPU_DEVICE_MANAGED amrex::Real mew;
    AMREX_GPU_DEVICE_MANAGED amrex::Real mew_w;

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

    AMREX_GPU_DEVICE_MANAGED amrex::Real small_number = 1.0e-15;
    AMREX_GPU_DEVICE_MANAGED amrex::Real large_number = 1.0e32;
    AMREX_GPU_DEVICE_MANAGED amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();

    AMREX_GPU_DEVICE_MANAGED amrex::Real neighborhood = 1.0;

    void Initialize ()
    {

        get_lsd_collision_coefficients(&NPHASE,
                                       &kt, &kt_w, &kn, &kn_w,
                                       &mew, &mew_w,
                                       &etan[0][0], &etan_w[0],
                                       &etat[0][0], &etat_w[0],
                                       &neighborhood);
    }


}
