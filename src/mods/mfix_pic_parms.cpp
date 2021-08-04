#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <mfix_pic_parms.H>

namespace PIC
{
  int solve = 0;

//  COLLISIONMODEL CollisionModel = LSD;
  int NPHASE;

  int verbose(0);

  amrex::Real Ps;
  amrex::Real beta;
  amrex::Real ep_cp;

  amrex::Real vel_ref_frame(0.5);

  amrex::Real small_number(1.e-7);

  amrex::Real damping_factor(0.4);
  amrex::Real damping_factor_wall_normal(0.3);
  amrex::Real damping_factor_wall_tangent(0.99);

  amrex::Real advance_vel_p = 0.5;

  int max_iter(3);
  InitialStepType initial_step(InitialStepType::Invalid);

  // Names of the solids used to build input regions.
  amrex::Vector<std::string> names;


  void Initialize ()
  {
    amrex::ParmParse pp("pic");

    pp.queryarr("solve", names);

    solve = names.size();

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(solve >= 0, "pic.solve must be >= 0");

    for(int lc=0; lc < names.size(); ++lc){
      if (amrex::toLower(names[0]).compare("none") == 0 ||
        (names[0]).compare("0") == 0) solve = 0;
    }

    // You can't name a solids "None" or "0" -- you just can't
    if( solve == 0 && names.size() > 1 ) {
      amrex::Abort("Invalid input: One or more PIC solids defined"
                    "but, the solver is disabled!");
    }

    NPHASE = solve;

    if(solve)
    {
      // verbosity level
      pp.query("verbose", verbose);

      // Read Ps
      pp.get("pressure_coefficient", Ps);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Ps > 0,
          "Invalid value: pic.Ps must be > 0");

      // Read beta
      pp.get("beta", beta);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(beta >= 2. && beta <= 5.,
          "Invalid value: pic.beta must be in [2.0, 5.0]");

      // Read close_pack coefficient
      pp.get("close_pack", ep_cp);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_cp > 0. && ep_cp < 1.,
          "Invalid value: pic.close_pack must be in [0.0, 1.0]");

      pp.get("damping_factor", damping_factor);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(damping_factor >= 0. && damping_factor <= 1.,
           "Invalid value: pic.damping_factor must be in [0.0, 1.0]");

      pp.get("damping_factor_wall_normal", damping_factor_wall_normal);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(damping_factor_wall_normal >= 0. && damping_factor_wall_normal <= 1.,
           "Invalid value: pic.damping_factor_wall_normal must be in [0.0, 1.0]");

      pp.get("damping_factor_wall_tangent", damping_factor_wall_tangent);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(damping_factor_wall_tangent >= 0. && damping_factor_wall_tangent <= 1.,
           "Invalid value: pic.damping_factor_wall_tangent must be in [0.0, 1.0]");

      // Read small number
      pp.query("small_number", small_number);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(small_number > 0,
          "Invalid value: pic.small_number must be > 0");

      // Solids slip velocity factor
      pp.query("velocity_reference_frame", vel_ref_frame);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vel_ref_frame >= 0. && vel_ref_frame <= 1.,
          "Invalid value: pic.velocity_reference_frame must be in [0.0, 1.0]");

      pp.query("advance_vel_p", advance_vel_p);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(advance_vel_p >= 0.0 && advance_vel_p <= 1.0,
          "Invalid value: pic.advance_vel_p must be in [0.0, 1.0]");

      pp.query("max_iter", max_iter);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_iter >= 1,
               "Invalid value: pic.max_iter must be greater than or equal to 1");

      std::string l_initial_step = "nth_eps";
      pp.query("initial_step_type", l_initial_step);
      if (l_initial_step == "nth_eps") {
        initial_step = InitialStepType::nth_eps;
      } else if (l_initial_step == "zero_eps") {
        initial_step = InitialStepType::zero_eps;
      } else if (l_initial_step == "taylor_approx") {
        initial_step = InitialStepType::taylor_approx;
      } else {
        amrex::Abort("pic.initial_step_type must be nth_eps, zero_eps or taylor_approx");
      }
    }
}
}
