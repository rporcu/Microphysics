#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <mfix_des_F.H>

#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <MFIX_DEM_Parms.H>

namespace DEM
{

    int solve = 0;

    AMREX_GPU_DEVICE_MANAGED COLLISIONMODEL CollisionModel = LSD;
    AMREX_GPU_DEVICE_MANAGED int NPHASE;

    AMREX_GPU_DEVICE_MANAGED amrex::Real mew;
    AMREX_GPU_DEVICE_MANAGED amrex::Real mew_w;

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
    AMREX_GPU_DEVICE_MANAGED amrex::Real en_w_input[NMAX];

    AMREX_GPU_DEVICE_MANAGED amrex::Real eta_fac;
    AMREX_GPU_DEVICE_MANAGED amrex::Real eta_w_fac;

    AMREX_GPU_DEVICE_MANAGED amrex::Real small_number = 1.0e-15;
    AMREX_GPU_DEVICE_MANAGED amrex::Real large_number = 1.0e32;
    AMREX_GPU_DEVICE_MANAGED amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();

    AMREX_GPU_DEVICE_MANAGED amrex::Real neighborhood = 1.0;

    // Tangential damping factor ratio
    amrex::Real eta_fac_pp = 0.5;
    amrex::Real eta_fac_pw = 0.5;

    void Initialize ()
    {

      amrex::ParmParse pp("dem");

      pp.get("solve", solve);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(solve >= 0, "dem.solve must be >= 0");
      NPHASE = solve;

      if( solve )
        {

         // Read MEW
         pp.get("friction_coeff_pp", mew);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= mew && mew <= 1.0,
              "Invalid value: dem.friction_coeff_pp must be in [0.0, 1.0]");

         // Read MEW_W
         pp.get("friction_coeff_pw", mew_w);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= mew_w && mew_w <= 1.0,
              "Invalid value: dem.friction_coeff_pw must be in [0.0, 1.0]");

         // Read KN
         pp.get("spring_const_pp", kn);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kn > 0.0,
              "Invalid value: dem.spring_const_pp must be > 0.0");

         // Read KN_w
         pp.get("spring_const_pw", kn_w);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kn_w > 0.0,
              "Invalid value: dem.spring_const_pw must be > 0.0");


         // Read DES_EN_INPUT
         amrex::Vector<amrex::Real> rest_coeff_pp_in;
         pp.queryarr("restitution_coeff_pp", rest_coeff_pp_in);

         // We know that we should have an upper-triangular matrix worth
         // of entries. (1-1, 1-2, 2-2, ...) for NPHASEs
         int req_coeffs = NPHASE+NPHASE*(NPHASE-1)/2;
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rest_coeff_pp_in.size() == req_coeffs,
              "Invalid number of entries: dem.restitution_coeff_pp");

         // Copy into the gloabl array that has a fixed size
         for (int idx=0; idx < req_coeffs; idx++){
           en_input[idx] = rest_coeff_pp_in[idx];
           amrex::Print() << "EN STUFF " << idx << "  " << en_input[idx] << "  " << rest_coeff_pp_in[idx] << "\n";
         }

         // Read DES_EN_WALL_INPUT
         amrex::Vector<amrex::Real> rest_coeff_pw_in;
         pp.queryarr("restitution_coeff_pw", rest_coeff_pw_in);

         // We know that we should have an entry for each NPHASE
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rest_coeff_pw_in.size() == NPHASE,
              "Invalid number of entries: dem.restitution_coeff_pw");

         // Copy into the gloabl array that has a fixed size
         for (int idx=0; idx < NPHASE; idx++)
           en_w_input[idx] = rest_coeff_pw_in[idx];


         // Read KT_FAC
        pp.get("spring_tang_fac_pp", kt_fac);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= kt_fac && kt_fac <= 1.0,
             "Invalid value: dem.spring_tang_fac_pp must be in [0.0, 1.0]");

        //Read KT_W_FAC
        pp.get("spring_tang_fac_pw", kt_w_fac);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= kt_w_fac && kt_w_fac <= 1.0,
             "Invalid value: dem.spring_tang_fac_pw must be in [0.0, 1.0]");

        // Read DES_ETA_FAC
        pp.get("damping_tang_fac_pp", eta_fac);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= eta_fac && eta_fac <= 1.0,
             "Invalid value: dem.damping_tang_fac_pp must be in [0.0, 1.0]");

        // Read DES_ETA_W_FAC
        pp.get("damping_tang_fac_pw", eta_w_fac);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= eta_w_fac && eta_w_fac <= 1.0,
             "Invalid value: dem.damping_tang_fac_pw must be in [0.0, 1.0]");


        // Calculate the tangential spring stiffness
        kt   = kt_fac   * kn;
        kt_w = kt_w_fac * kn_w;


        set_lsd_collision_coefficients(&NPHASE,
                                       &mew,         &mew_w,
                                       &kn,          &kn_w,
                                       &kt,          &kt_w,
                                       &en_input[0], &en_w_input[0],
                                       &kt_fac,      &kt_w_fac,
                                       &eta_fac,     &eta_w_fac);

      }

    }

}
