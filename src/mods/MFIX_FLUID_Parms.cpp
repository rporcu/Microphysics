#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>
#include <MFIX_FLUID_Parms.H>

namespace FLUID
{

  bool solve;

  // Specified constant gas density
  amrex::Real ro_g0;

  // Specified constant tracer value
  amrex::Real trac_0;

  // Specified constant gas viscosity
  amrex::Real mu_g0;

  // Average molecular weight of gas
  amrex::Real mw_avg;

  void Initialize ()
  {

    amrex::ParmParse pp("fluid");

    pp.get("solve",   solve );

    if( solve ) {

      pp.get("density",   ro_g0 );
      pp.get("viscosity", mu_g0 );
      pp.query("trac0",  trac_0 );
      pp.query("mw_avg", mw_avg );

    }

  }

}
