#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <mfix_dem_parms.H>

namespace DEM
{
    int solve = 0;

    COLLISIONMODEL CollisionModel = LSD;
    int NPHASE;

    amrex::Real dtsolid;

    amrex::Real mew;
    amrex::Real mew_w;

    amrex::Real kt;
    amrex::Real kt_w;

    amrex::Real kn;
    amrex::Real kn_w;

    amrex::Real kt_fac = 2.0/7.0;
    amrex::Real kt_w_fac = 2.0/7.0;

    // normal and tangential components of the damping coefficients
    A2D etan;
    A1D etan_w;

    A2D etat;
    A1D etat_w;

    // coefficients of restitution, normal and tangential
    A2D en;
    A1D en_w;

    amrex::Real eta_fac;
    amrex::Real eta_w_fac;

    amrex::Real small_number = 1.0e-15;
    amrex::Real large_number = 1.0e32;
    amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();

    amrex::Real neighborhood = 1.0;

    // Tangential damping factor ratio
    amrex::Real eta_fac_pp = 0.5;
    amrex::Real eta_fac_pw = 0.5;

    // Coarse-grain DEM
    int cg_dem = 0;

    int restart_from_PIC = 0;

    void Initialize ()
    {

      amrex::ExecOnFinalize(Finalize);

      etan.alloc();
      etan_w.alloc();
      etat.alloc();
      etat_w.alloc();
      en.alloc();
      en_w.alloc();

      amrex::ParmParse ppDEM("dem");

      // Names of the solids used to build input regions.
      amrex::Vector<std::string> names;

      ppDEM.queryarr("solve", names);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(names.size() >= 1,
           "DEM solver not specified: Input dem.solve is undefined!");

      solve = 1;
      for(int lc=0; lc < names.size(); ++lc){
        if (amrex::toLower(names[0]).compare("none") == 0 ||
          (names[0]).compare("0") == 0) solve = 0;
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
        ppDEM.get("friction_coeff.pp", mew);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= mew && mew <= 1.0,
             "Invalid value: dem.friction_coeff.pp must be in [0.0, 1.0]");

        // Read MEW_W
        ppDEM.get("friction_coeff.pw", mew_w);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= mew_w && mew_w <= 1.0,
             "Invalid value: dem.friction_coeff.pw must be in [0.0, 1.0]");

        // Read KN
        ppDEM.get("spring_const.pp", kn);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kn > 0.0,
             "Invalid value: dem.spring_const.pp must be > 0.0");

        // Read KN_w
        ppDEM.get("spring_const.pw", kn_w);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kn_w > 0.0,
             "Invalid value: dem.spring_const.pw must be > 0.0");

         // Read KT_FAC
        ppDEM.get("spring_tang_fac.pp", kt_fac);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= kt_fac && kt_fac <= 1.0,
             "Invalid value: dem.spring_tang_fac.pp must be in [0.0, 1.0]");

        //Read KT_W_FAC
        ppDEM.get("spring_tang_fac.pw", kt_w_fac);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= kt_w_fac && kt_w_fac <= 1.0,
             "Invalid value: dem.spring_tang_fac.pw must be in [0.0, 1.0]");

        // Read DES_ETA_FAC
        ppDEM.get("damping_tang_fac.pp", eta_fac);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= eta_fac && eta_fac <= 1.0,
             "Invalid value: dem.damping_tang_fac.pp must be in [0.0, 1.0]");

        // Read DES_ETA_W_FAC
        ppDEM.get("damping_tang_fac.pw", eta_w_fac);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= eta_w_fac && eta_w_fac <= 1.0,
             "Invalid value: dem.damping_tang_fac.pw must be in [0.0, 1.0]");


        // Calculate the tangential spring stiffness
        kt   = kt_fac   * kn;
        kt_w = kt_w_fac * kn_w;

        //// We know that we should have an upper-triangular matrix worth
        //// of entries. (1-1, 1-2, 2-2, ...) for NPHASEs
        //int req_coeffs = NPHASE+NPHASE*(NPHASE-1)/2;


        // Read restitution coefficients. These can be given two different ways:
        //    dem.restitution_coeff.solid1.solid2 = coeff
        //    dem.restitution_coeff.solid2.solid1 = coeff
        //
        // We want to make sure that at least one is given. If both are given
        // then they must be equal.  The values get stored in en.
        {
          A2D::array_type host_en;
          A1D::array_type host_en_w;

          amrex::ParmParse ppRC("dem.restitution_coeff");
          for (int idx0=0; idx0 < NPHASE; idx0++){
            for (int idx1=idx0; idx1 < NPHASE; idx1++){

              std::string pp01 = names[idx0]+"."+names[idx1];
              amrex::Real coeff01 = -1.0;
              ppRC.query(pp01.c_str(), coeff01);

              std::string pp10 = names[idx1]+"."+names[idx0];
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

              host_en(idx0,idx1) = rest_coeff;
              host_en(idx1,idx0) = rest_coeff;

            }

            std::string pp01 = names[idx0]+".wall";
            amrex::Real coeff01 = -1.0;
            ppRC.query(pp01.c_str(), coeff01);

            std::string pp10 = "wall."+names[idx0];
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

            host_en_w(idx0) = rest_coeff;

          }
#ifdef AMREX_USE_GPU
          amrex::Gpu::htod_memcpy_async(DEM::en.arrayPtr(), &host_en, sizeof(DEM::A2D::array_type));
          amrex::Gpu::htod_memcpy_async(DEM::en_w.arrayPtr(), &host_en_w, sizeof(DEM::A1D::array_type));
          amrex::Gpu::synchronize();
#else
          *(DEM::en.arrayPtr()) = host_en;
          *(DEM::en_w.arrayPtr()) = host_en_w;
#endif
        }


        // Read coarse-grain DEM
        ppDEM.query("coarse_grain", cg_dem);

        // Read PIC to DEM parameters
        ppDEM.query("restart_from_PIC", restart_from_PIC);

        if (restart_from_PIC) {

          amrex::ParmParse ppAMR("amr");
          std::string restart_file {""};

          ppAMR.get("restart", restart_file);

          const int is_restarting = !(restart_file.empty());

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(is_restarting,
              "Invalid attempt to restart from PIC from an empty chk file");
        }
      }

    }

    void Finalize ()
    {
      etan.free();
      etan_w.free();
      etat.free();
      etat_w.free();
      en.free();
      en_w.free();
    }
}
