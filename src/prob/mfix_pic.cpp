#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <mfix_pic.H>

using namespace amrex;


MFIXPIC::MFIXPIC ()
  : m_solve(0)
  , m_verbose(0)
  , m_vel_ref_frame(0.5)
  , m_small_number(1.e-7)
  , m_damping_factor(0.4)
  , m_damping_factor_wall_normal(0.3)
  , m_damping_factor_wall_tangent(0.99)
  , m_advance_vel_p(0.5)
  , m_max_iter(3)
  , m_initial_step(InitialStepType::Invalid)
  , m_restart_refinement(1)
{}


void
MFIXPIC::Initialize ()
{
  amrex::ParmParse pp("pic");

  // Names of the solids used to build input regions.
  amrex::Vector<std::string> names;
  pp.queryarr("solve", names);

  m_solve = names.size();

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_solve >= 0, "pic.m_solve must be >= 0");

  for(int lc=0; lc < names.size(); ++lc){
    if (amrex::toLower(names[0]).compare("none") == 0 ||
      (names[0]).compare("0") == 0) m_solve = 0;
  }

  // You can't name a solids "None" or "0" -- you just can't
  if (m_solve == 0 && names.size() > 1) {
    amrex::Abort("Invalid input: One or more PIC solids defined"
                  "but, the m_solver is disabled!");
  }

  if (m_solve) {
    // Store the total number of solids
    m_NPHASE = names.size();
  }

  if(m_solve)
  {
    // verbosity level
    pp.query("verbose", m_verbose);

    // Read Ps
    pp.get("pressure_coefficient", m_Ps);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_Ps > 0,
        "Invalid value: pic.Ps() must be > 0");

    // Read beta
    pp.get("beta", m_beta);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_beta >= 2. && m_beta <= 5.,
        "Invalid value: pic.beta() must be in [2.0, 5.0]");

    // Read close_pack coefficient
    pp.get("close_pack", m_ep_cp);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_ep_cp > 0. && m_ep_cp < 1.,
        "Invalid value: pic.close_pack must be in [0.0, 1.0]");

    pp.get("damping_factor", m_damping_factor);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_damping_factor >= 0. && m_damping_factor <= 1.,
         "Invalid value: pic.damping_factor() must be in [0.0, 1.0]");

    pp.get("damping_factor_wall_normal", m_damping_factor_wall_normal);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_damping_factor_wall_normal >= 0. && m_damping_factor_wall_normal <= 1.,
         "Invalid value: pic.damping_factor()_wall_normal must be in [0.0, 1.0]");

    pp.get("damping_factor_wall_tangent", m_damping_factor_wall_tangent);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_damping_factor_wall_tangent >= 0. && m_damping_factor_wall_tangent <= 1.,
         "Invalid value: pic.damping_factor()_wall_tangent must be in [0.0, 1.0]");

    // Read small number
    pp.query("small_number", m_small_number);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_small_number > 0,
        "Invalid value: pic.small_number() must be > 0");

    // Solids slip velocity factor
    pp.query("velocity_reference_frame", m_vel_ref_frame);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_vel_ref_frame >= 0. && m_vel_ref_frame <= 1.,
        "Invalid value: pic.velocity_reference_frame must be in [0.0, 1.0]");

    pp.query("advance_vel_p", m_advance_vel_p);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_advance_vel_p >= 0.0 && m_advance_vel_p <= 1.0,
        "Invalid value: pic.advance_vel_p() must be in [0.0, 1.0]");

    pp.query("max_iter", m_max_iter);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_max_iter >= 1,
             "Invalid value: pic.max_iter() must be greater than or equal to 1");

    std::string l_initial_step = "nth_eps";
    pp.query("initial_step_type", l_initial_step);
    if (l_initial_step == "nth_eps") {
      m_initial_step = InitialStepType::nth_eps;
    } else if (l_initial_step == "zero_eps") {
      m_initial_step = InitialStepType::zero_eps;
    } else if (l_initial_step == "taylor_approx") {
      m_initial_step = InitialStepType::taylor_approx;
    } else {
      amrex::Abort("pic.initial_step()_type must be nth_eps, zero_eps or taylor_approx");
    }

    pp.query("restart_refinement", m_restart_refinement);

    if (m_restart_refinement > 1) {

      amrex::ParmParse ppAMR("amr");
      std::string restart_file {""};

      ppAMR.get("restart", restart_file);

      const int is_restarting = !(restart_file.empty());

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(is_restarting,
          "Invalid attempt to set a PIC restart refinement without actually restarting from a chk file");
    }
  }
}
