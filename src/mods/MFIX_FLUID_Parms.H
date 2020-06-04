#ifndef MFIX_FLUID_PARMS_H_
#define MFIX_FLUID_PARMS_H_

#include <AMReX_REAL.H>

#include <MFIX_SPECIES_Parms.H>

namespace FLUID
{

  enum DENSITYMODEL {ConstantDensity, IdealGas};

  enum VISCOSITYMODEL {ConstantViscosity, Sutherland};

  enum SPECIFICHEATMODEL {ConstantSpecificHeat, NASA9Polynomials};

  enum THERMALCONDUCTIVITYMODEL {ConstantThermalConductivity};

  // Flag to solve fluid equations
  extern int solve;

  // Specified constant gas temperature
  extern amrex::Real T_g0;

  // Specified constant gas specific heat
  extern amrex::Real cp_g0;

  // Gas phase thermal conductivity coefficient
  extern amrex::Real k_g0;

  // Specified constant gas density
  extern amrex::Real ro_g0;

  // Specified constant tracer value
  extern amrex::Real trac_0;

  // Specified constant gas viscosity
  extern amrex::Real mu_g0;

  // Average molecular weight of gas
  extern amrex::Real mw_avg;

  // Flag to solve species fluid equations
  extern int solve_species;

  // Fluid phase species names
  extern std::vector<std::string> species_g;

  // Total number of fluid species
  extern int nspecies_g;

  // Fluid phase species diffusion coefficients
  extern std::vector<amrex::Real> D_g0;

  // Fluid phase species mass fraction
  extern std::vector<amrex::Real> X_g0;

  // Name to later reference when building input regions
  extern std::string name;

  struct FLUID_t {

    amrex::Real volfrac;
    amrex::Real temperature;
    amrex::Real pressure;
    amrex::Vector<amrex::Real> velocity;

    amrex::Real pressure_defined; // Flag to check if pressure is defined

    SPECIES::SPECIES_t species;

    FLUID_t() :
      volfrac(-1.0),
      temperature(0.0),
      pressure(-1.0),
      species(FLUID::nspecies_g)
    {}

  };

  void Initialize ();

}

#endif
