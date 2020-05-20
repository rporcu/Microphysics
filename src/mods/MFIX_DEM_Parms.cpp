#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <MFIX_DEM_Parms.H>

namespace DEM
{
    int solve = 0;

    AMREX_GPU_DEVICE_MANAGED COLLISIONMODEL CollisionModel = LSD;
    AMREX_GPU_DEVICE_MANAGED int NPHASE;

    AMREX_GPU_DEVICE_MANAGED amrex::Real dtsolid;

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

    // Names of the solids used to build input regions.
    amrex::Vector<std::string> names;

    // Particle species
    amrex::Vector<std::string> species_dem;

    // Particle species fractions
    amrex::Vector<amrex::Real> spec_frac_dem;
    // Number of species at each particle
    AMREX_GPU_DEVICE_MANAGED int nspecies_dem = 0;

    // Coarse-grain DEM
    AMREX_GPU_DEVICE_MANAGED int cg_dem = 0;

    void Initialize ()
    {

      amrex::ParmParse pp("dem");

      pp.queryarr("solve", names);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(names.size() >= 1,
           "DEM solver not specified: Input dem.solve is undefined!");

      solve = 1;
      for(int lc=0; lc < names.size(); ++lc){
        if (names[0] == "None" ||
            names[0] == "none" ||
            names[0] == "NONE" ||
            names[0] == "0" ) solve = 0;
      }

      // You can't name a solids "None" or "0" -- you just can't
      if( solve == 0 && names.size() > 1 ){
        amrex::Abort("Invalid input: One or more DEM solids defined"
                     "but, the solver is disabled!");
      }




      //TODO: Add check to prevent using the same name twice.


      if( solve ) {

        // Store the total number of solids
        NPHASE = names.size();

        // Read MEW
        pp.get("friction_coeff.pp", mew);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= mew && mew <= 1.0,
             "Invalid value: dem.friction_coeff.pp must be in [0.0, 1.0]");

        // Read MEW_W
        pp.get("friction_coeff.pw", mew_w);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= mew_w && mew_w <= 1.0,
             "Invalid value: dem.friction_coeff.pw must be in [0.0, 1.0]");

        // Read KN
        pp.get("spring_const.pp", kn);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kn > 0.0,
             "Invalid value: dem.spring_const.pp must be > 0.0");

        // Read KN_w
        pp.get("spring_const.pw", kn_w);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kn_w > 0.0,
             "Invalid value: dem.spring_const.pw must be > 0.0");

         // Read KT_FAC
        pp.get("spring_tang_fac.pp", kt_fac);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= kt_fac && kt_fac <= 1.0,
             "Invalid value: dem.spring_tang_fac.pp must be in [0.0, 1.0]");

        //Read KT_W_FAC
        pp.get("spring_tang_fac.pw", kt_w_fac);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= kt_w_fac && kt_w_fac <= 1.0,
             "Invalid value: dem.spring_tang_fac.pw must be in [0.0, 1.0]");

        // Read DES_ETA_FAC
        pp.get("damping_tang_fac.pp", eta_fac);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= eta_fac && eta_fac <= 1.0,
             "Invalid value: dem.damping_tang_fac.pp must be in [0.0, 1.0]");

        // Read DES_ETA_W_FAC
        pp.get("damping_tang_fac.pw", eta_w_fac);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= eta_w_fac && eta_w_fac <= 1.0,
             "Invalid value: dem.damping_tang_fac.pw must be in [0.0, 1.0]");


        // Calculate the tangential spring stiffness
        kt   = kt_fac   * kn;
        kt_w = kt_w_fac * kn_w;

        // Read species inputs -----------------------------------------//
        int species = 0;
        pp.query("species", species);
        if (species > 0)
        {
             pp.getarr("species.names", species_dem);
             nspecies_dem = species_dem.size();
             pp.getarr("species.fractions", spec_frac_dem);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE( species_dem.size() == spec_frac_dem.size(),
             "Species fraction number does not match species number");
        }

        // Read coarse-grain DEM
        pp.query("coarse_grain", cg_dem);

        //// We know that we should have an upper-triangular matrix worth
        //// of entries. (1-1, 1-2, 2-2, ...) for NPHASEs
        //int req_coeffs = NPHASE+NPHASE*(NPHASE-1)/2;


        // Read restitution coefficients. These can be given two different ways:
        //    dem.restitution_coeff.solid1.solid2 = coeff
        //    dem.restitution_coeff.solid2.solid1 = coeff
        //
        // We want to make sure that at least one is given. If both are given
        // then they must be equal.  The values get stroed in en_input


        {
          int lc_pp = 0;
          int lc_pw = 0;

          amrex::ParmParse ppRC("dem.restitution_coeff");
          for (int idx0=0; idx0 < NPHASE; idx0++){
            for (int idx1=idx0; idx1 < NPHASE; idx1++){

              std::string pp01 = DEM::names[idx0]+"."+DEM::names[idx1];
              amrex::Real coeff01 = -1.0;
              ppRC.query(pp01.c_str(), coeff01);

              std::string pp10 = DEM::names[idx1]+"."+DEM::names[idx0];
              amrex::Real coeff10 = -1.0;
              ppRC.query(pp10.c_str(), coeff10);

              // Set the temp variable to something we can check against.
              amrex::Real rest_coeff(-1.0);

              // Take either one if they are the same. Otherwise, take
              // the "other one" if we see that one isn't set.
              if( coeff01 == coeff10 ){
                rest_coeff = coeff01;

              } else if ( coeff01 == -1.0 ) {
                rest_coeff = coeff10;

              } else if ( coeff10 == -1.0 ) {
                rest_coeff = coeff01;

              }
              // There is no need for an 'else' here. This is covered
              // by initializating the variable to an invalid value.

              AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rest_coeff >= 0.0 && rest_coeff <= 1.0,
                   "Invalid restitution coefficient.");

              en_input[lc_pp] = rest_coeff;
              lc_pp += 1;

            }

            std::string pp01 = DEM::names[idx0]+".wall";
            amrex::Real coeff01 = -1.0;
            ppRC.query(pp01.c_str(), coeff01);

            std::string pp10 = "wall."+DEM::names[idx0];
            amrex::Real coeff10 = -1.0;
            ppRC.query(pp10.c_str(), coeff10);

            // Set the temp variable to something we can check against.
            amrex::Real rest_coeff(-1.0);

            // Take either one if they are the same. Otherwise, take
            // the "other one" if we see that one isn't set.
            if( coeff01 == coeff10 ){
              rest_coeff = coeff01;

            } else if ( coeff01 == -1.0 ) {
              rest_coeff = coeff10;

            } else if ( coeff10 == -1.0 ) {
              rest_coeff = coeff01;

            }
            // There is no need for an 'else' here. This is covered
            // by initializating the variable to an invalid value.

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rest_coeff >= 0.0 && rest_coeff <= 1.0,
                 "Invalid restitution coefficient.");

            en_w_input[lc_pw] = rest_coeff;
            lc_pw += 1;

          }
        }

      }

    }

}
