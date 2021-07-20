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

  amrex::Real Ps;
  amrex::Real beta;
  amrex::Real ep_cp;

  amrex::Real velfac = 0.8;
  amrex::Real velfac_alpha = 0.0625;
  amrex::Real velfac_intercept = 0.85;

  amrex::Real small_number = 1.e-7;

  amrex::Real damping_factor = 0.4;
  amrex::Real damping_factor_wall_normal  = 0.30;
  amrex::Real damping_factor_wall_tangent = 0.99;

  amrex::Real sigmoidal_offset;

  amrex::Real advance_eps_step = 0.5;

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
      pp.query("velfac", velfac);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(velfac >= 0. && velfac <= 1.,
          "Invalid value: pic.velfac must be in [0.0, 1.0]");

      pp.query("velfac_alpha", velfac_alpha);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(velfac_alpha > 0. and velfac_alpha <= 0.5,
          "Invalid value: pic.velfac_alpha must be in (0.0, 0.5]");

      pp.query("velfac_intercept", velfac_intercept);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(velfac_intercept >= velfac && velfac_intercept < 1.0,
          "Invalid value: pic.velfac_intercept must be in (velfac, 1.0)");

      {
        const amrex::Real b = -9.0*(2.0*(velfac_intercept - velfac)/(1.0 - velfac) - 1.0);
        const amrex::Real d =  3.0*b;

        const amrex::Real p = -b/3.0;
        const amrex::Real q = p*p*p + (27*b-3.0*d)/6.0;
        const amrex::Real r = 9.0;

        const amrex::Real qrp = std::sqrt(q*q + (r-p*p)*(r-p*p)*(r-p*p));

        sigmoidal_offset = (std::cbrt(q + qrp) + std::cbrt(q - qrp) + p)/3.0;
      }

      pp.query("advance_eps_step", advance_eps_step);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(advance_eps_step >= 0.0 && advance_eps_step <= 1.0,
          "Invalid value: pic.advance_eps_step must be in [0.0, 1.0]");

    }
  }
}
