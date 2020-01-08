#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>
#include <MFIX_FLUID_Parms.H>

namespace FLUID
{

  int solve;

  // Specified constant gas density
  amrex::Real ro_g0;

  // Specified constant tracer value
  amrex::Real trac_0 = 0.0;

  // Specified constant gas viscosity
  amrex::Real mu_g0;

  // Average molecular weight of gas
  amrex::Real mw_avg;

  // Name to later reference when building inputs for IC/BC regions
  std::string name;


  void Initialize ()
  {

    amrex::ParmParse pp("fluid");

    std::vector<std::string> fluid_name;
    pp.queryarr("solve", fluid_name);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid_name.size() == 1,
       "Fluid solver not specified. fluid.sove = ? ");

    // Disable the fluid solver if the fluid is defined as "None"
    if (fluid_name[0] == "None" ||
        fluid_name[0] == "none" ||
        fluid_name[0] == "NONE" ||
        fluid_name[0] == "0" ) {
      solve = 0;
    } else {
      solve = 1;
      name = fluid_name[0];
    }


    if( solve ) {

      amrex::ParmParse ppFluid(name.c_str());

      // Get density inputs ------------------------------------//
      std::string density_model;
      ppFluid.get("density", density_model );

      if (density_model == "constant") {

        DENSITYMODEL DensityModel = ConstantDensity;
        ppFluid.get("density.constant",   ro_g0 );

      } else {

        amrex::Abort("Unknown fluid density model!");

      }


      // Get viscosity inputs ----------------------------------//
      std::string viscosity_model;
      ppFluid.get("viscosity", viscosity_model );

      if (density_model == "constant") {

        VISCOSITYMODEL ViscosityModel = ConstantViscosity;
        ppFluid.get("viscosity.constant", mu_g0 );


      } else {

        amrex::Abort("Unknown fluid viscosity model!");

      }

      pp.query("trac0",  trac_0 );
      pp.query("mw_avg", mw_avg );

    }

  }

}
