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

  void Initialize ()
  {

    amrex::ParmParse pp("fluid");

    pp.get("solve", solve);

    if( solve ) {

      // Get density inputs ------------------------------------//
      std::string density_model;
      pp.get("density_model", density_model );

      if (density_model == "constant") {

        DENSITYMODEL DensityModel = ConstantDensity;
        pp.get("constant_density",   ro_g0 );

      } else {

        amrex::Abort("Unknown fluid density model!");

      }


      // Get viscosity inputs ----------------------------------//
      std::string viscosity_model;
      pp.get("viscosity_model", viscosity_model );

      if (density_model == "constant") {

        VISCOSITYMODEL ViscosityModel = ConstantViscosity;
        pp.get("constant_viscosity", mu_g0 );


      } else {

        amrex::Abort("Unknown fluid viscosity model!");

      }

      pp.query("trac0",  trac_0 );
      pp.query("mw_avg", mw_avg );

    }

  }

}
