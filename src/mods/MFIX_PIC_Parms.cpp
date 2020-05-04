#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>

#include <MFIX_PIC_Parms.H>

namespace PIC
{
  int solve = 0;

//  AMREX_GPU_DEVICE_MANAGED COLLISIONMODEL CollisionModel = LSD;
  AMREX_GPU_DEVICE_MANAGED int NPHASE;

  AMREX_GPU_DEVICE_MANAGED amrex::Real Ps;
  AMREX_GPU_DEVICE_MANAGED amrex::Real beta;
  AMREX_GPU_DEVICE_MANAGED amrex::Real ep_cp;

  AMREX_GPU_DEVICE_MANAGED amrex::Real velfac;

  AMREX_GPU_DEVICE_MANAGED amrex::Real small_number = 1.e-7;

  AMREX_GPU_DEVICE_MANAGED amrex::Real damping_factor;

  // Names of the solids used to build input regions.
  amrex::Vector<std::string> names;

  void Initialize ()
  {
    amrex::ParmParse pp("pic");

    pp.queryarr("solve", names);

    solve = names.size();

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(solve >= 0, "pic.solve must be >= 0");
    NPHASE = solve;

    if(solve)
    {

      // Read Ps
      pp.get("pressure_coefficient", Ps);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Ps > 0,
          "Invalid value: pic.Ps must be > 0");

      // Read beta
      pp.get("beta", beta);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(beta >= 2. and beta <= 5.,
          "Invalid value: pic.beta must be in [2.0, 5.0]");

      // Read close_pack coefficient
      pp.get("close_pack", ep_cp);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_cp > 0. and ep_cp < 1.,
          "Invalid value: pic.close_pack must be in [0.0, 1.0]");

      // Read small number
      pp.query("small_number", small_number, 1.e-7);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(small_number > 0,
          "Invalid value: pic.small_number must be > 0");

      // Solids slip velocity factor
      velfac = 1.0;
      pp.query("velfac", velfac);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(velfac >= 0. and velfac <= 1.,
          "Invalid value: pic.velfac must be in [0.0, 1.0]");

      damping_factor = 0.85;
      pp.get("damping_factor", damping_factor);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(damping_factor >= 0. and damping_factor <= 1.,
           "Invalid value: pic.damping_factor must be in [0.0, 1.0]");


    }
  }
}
