#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <mfix_dem.H>

using namespace amrex;


MFIXDEM::MFIXDEM ()
  : m_collision_model(CollisionModel::LSD)
  , m_solve(0)
  , m_kt_fac(2.0/7.0)
  , m_kt_w_fac(2.0/7.0)
  , m_eta_fac_pp(0.5)
  , m_eta_fac_pw(0.5)
  , m_small_number(1.0e-15)
  , m_large_number(1.0e32)
  , m_eps(std::numeric_limits<amrex::Real>::epsilon())
  , m_neighborhood(1.0)
  , m_cg_dem(0)
  , m_restart_from_PIC(0)
{}


MFIXDEM::~MFIXDEM ()
{
  m_etan.free();
  m_etan_w.free();
  m_etat.free();
  m_etat_w.free();
  m_en.free();
  m_en_w.free();
}


void
MFIXDEM::Initialize ()
{
  m_etan.alloc();
  m_etan_w.alloc();
  m_etat.alloc();
  m_etat_w.alloc();
  m_en.alloc();
  m_en_w.alloc();

  // Fluid conductivity hack
  amrex::ParmParse ppKG("fluid");
  ppKG.query("thermal_conductivity.constant", m_k_g_dem);

  amrex::ParmParse ppDEM("dem");

  // Polydisperse neighbor search
  ppDEM.query("PolyNeighSearch" , m_pneigh_flag);
  if (m_pneigh_flag) {
#if !(MFIX_POLYDISPERSE)
      amrex::Abort("MFIX was not built with POLYDISPERSE support");
#endif
      ppDEM.query("PolyNumTypes"    , m_nptypes);
      m_prefrats.resize(m_nptypes);
      ppDEM.queryarr("PolyRefRatios", m_prefrats);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_prefrats.back() == 1,
                                       "The largest polydisperse refinement ratio must be 1.");
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_nptypes >= 2,
                                       "Cannot use polydisperse neighbor algorithm with less than 2 types.");
  }


  // Names of the solids used to build input regions.
  amrex::Vector<std::string> names;

  ppDEM.queryarr("solve", names);

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(names.size() >= 1,
       "DEM solver not specified: Input dem.solve() is undefined!");

  m_solve = 1;
  for(int lc=0; lc < names.size(); ++lc){
    if (amrex::toLower(names[0]).compare("none") == 0 || (names[0]).compare("0") == 0)
      m_solve = 0;
  }

  // You can't name a solids "None" or "0" -- you just can't
  if (m_solve == 0 && names.size() > 1) {
    amrex::Abort("Invalid input: One or more DEM solids defined"
                 "but, the solver is disabled!");
  }

  //TODO: Add check to prevent using the same name twice.

  if (m_solve) {

    // Store the total number of solids
    m_NPHASE = names.size();

    // Read MEW
    ppDEM.get("friction_coeff.pp", m_mew);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= m_mew && m_mew <= 1.0,
         "Invalid value: dem.friction_coeff.pp must be in [0.0, 1.0]");

    // Read MEW_W
    ppDEM.get("friction_coeff.pw", m_mew_w);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= m_mew_w && m_mew_w <= 1.0,
         "Invalid value: dem.friction_coeff.pw must be in [0.0, 1.0]");

    // Read KN
    ppDEM.get("spring_const.pp", m_kn);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_kn > 0.0,
         "Invalid value: dem.spring_const.pp must be > 0.0");

    // Read KN_w
    ppDEM.get("spring_const.pw", m_kn_w);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_kn_w > 0.0,
         "Invalid value: dem.spring_const.pw must be > 0.0");

     // Read KT_FAC
    ppDEM.get("spring_tang_fac.pp", m_kt_fac);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= m_kt_fac && m_kt_fac <= 1.0,
         "Invalid value: dem.spring_tang_fac.pp must be in [0.0, 1.0]");

    //Read KT_W_FAC
    ppDEM.get("spring_tang_fac.pw", m_kt_w_fac);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= m_kt_w_fac && m_kt_w_fac <= 1.0,
         "Invalid value: dem.spring_tang_fac.pw must be in [0.0, 1.0]");

    // Read DES_ETA_FAC
    ppDEM.get("damping_tang_fac.pp", m_eta_fac);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= m_eta_fac && m_eta_fac <= 1.0,
         "Invalid value: dem.damping_tang_fac.pp must be in [0.0, 1.0]");

    // Read DES_ETA_W_FAC
    ppDEM.get("damping_tang_fac.pw", m_eta_w_fac);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0.0 <= m_eta_w_fac && m_eta_w_fac <= 1.0,
         "Invalid value: dem.damping_tang_fac.pw must be in [0.0, 1.0]");


    // Calculate the tangential spring stiffness
    m_kt   = m_kt_fac   * m_kn;
    m_kt_w = m_kt_w_fac * m_kn_w;

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
      for (int idx0=0; idx0 < m_NPHASE; idx0++){
        for (int idx1=idx0; idx1 < m_NPHASE; idx1++){

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
          if (coeff01 == coeff10) {
            rest_coeff = coeff01;

          } else if (coeff01 == -1.0) {
            rest_coeff = coeff10;

          } else if (coeff10 == -1.0) {
            rest_coeff = coeff01;

          }
          // There is no need for an 'else' here. This is covered
          // by initializing the variable to an invalid value.

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
        if (coeff01 == coeff10) {
          rest_coeff = coeff01;

        } else if (coeff01 == -1.0) {
          rest_coeff = coeff10;

        } else if (coeff10 == -1.0) {
          rest_coeff = coeff01;

        }
        // There is no need for an 'else' here. This is covered
        // by initializing the variable to an invalid value.

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rest_coeff >= 0.0 && rest_coeff <= 1.0,
             "Invalid restitution coefficient.");

        host_en_w(idx0) = rest_coeff;

      }
#ifdef AMREX_USE_GPU
      amrex::Gpu::htod_memcpy_async(m_en.arrayPtr(), &host_en, sizeof(A2D::array_type));
      amrex::Gpu::htod_memcpy_async(m_en_w.arrayPtr(), &host_en_w, sizeof(A1D::array_type));
      amrex::Gpu::synchronize();
#else
      *(m_en.arrayPtr()) = host_en;
      *(m_en_w.arrayPtr()) = host_en_w;
#endif
    }

    // Read coarse-grain DEM
    ppDEM.query("coarse_grain", m_cg_dem);

    // Read PIC to DEM parameters
    ppDEM.query("restart_from_PIC", m_restart_from_PIC);

    if (m_restart_from_PIC) {

      std::string restart_file {""};
      ParmParse("mfix").query("restart", restart_file);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!restart_file.empty(),
          "Invalid attempt to restart from PIC from an empty chk file");
    }
  }
}
