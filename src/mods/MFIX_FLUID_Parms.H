#ifndef MFIX_FLUID_PARMS_H_
#define MFIX_FLUID_PARMS_H_

#include <AMReX_REAL.H>

namespace FLUID
{

  enum DENSITYMODEL {ConstantDensity, IdealGas};

  enum VISCOSITYMODEL {ConstantViscosity, Sutherland};

  // Flag to solve fluid equations
  extern int solve;

  // Specified constant gas density
  extern amrex::Real ro_g0;

  // Specified constant tracer value
  extern amrex::Real trac_0;

  // Specified constant gas viscosity
  extern amrex::Real mu_g0;

  // Average molecular weight of gas
  extern amrex::Real mw_avg;

  // Name to later reference when building input regions
  extern std::string name;

  struct FLUID_t {

    amrex::Real volfrac;
    amrex::Real temperature;
    amrex::Real pressure;
    amrex::Vector<amrex::Real> velocity;

    amrex::Real pressure_defined; // Flag to check if pressure is defined

    FLUID_t() :
      volfrac(-1.0),
      pressure(-1.0),
      temperature(0.0) { }

  };

  void Initialize ();

}

#endif
