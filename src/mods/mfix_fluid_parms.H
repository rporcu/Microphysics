#ifndef MFIX_FLUID_PARMS_H_
#define MFIX_FLUID_PARMS_H_

#include <AMReX_REAL.H>

#include <mfix_species_parms.H>

namespace FLUID
{
  // List of recognized fluid density models
  struct DENSITYMODEL { enum {Invalid=-1, Constant, IdealGas}; };
  
  extern int DensityModel;

  // List of recognized fluid molecular weight models
  struct MOLECULARWEIGHTMODEL { enum {Invalid=-1, Constant, Mixture}; };
  
  extern int MolecularWeightModel;

  // List of recognized fluid viscosity models
  struct VISCOSITYMODEL { enum {Invalid=-1, Constant, Sutherland}; };
  
  extern int ViscosityModel;

  // List of recognized fluid specific heat models
  struct SPECIFICHEATMODEL { enum {Invalid=-1, Constant, NASA9Polynomials}; };
  
  extern int SpecificHeatModel;

  // List of recognized fluid thermal conductivity models
  struct THERMALCONDUCTIVITYMODEL{ enum {Invalid=-1, Constant}; };
  
  extern int ThermalConductivityModel;

  static constexpr amrex::Real R = 8.31446261815324;

  // Flag to solve fluid equations
  extern int solve;

  // Specified constant gas density
  extern amrex::Real ro_g0;

  // Specified constant molecular weight of gas
  extern amrex::Real MW_g0;

  // Specified constant tracer value
  extern amrex::Real trac_0;

  // Specified constant gas viscosity
  extern amrex::Real mu_g0;

  // Average molecular weight of gas
  extern amrex::Real mw_avg;

  // Flag to solve enthalpy fluid equation
  extern int solve_enthalpy;

  // Specified constant gas temperature
  extern amrex::Real T_g0;

  // Specified constant gas specific heat
  extern amrex::Real cp_g0;

  // Specified constant gas phase thermal conductivity coefficient
  extern amrex::Real k_g0;

  // Flag to solve species fluid equations
  extern int solve_species;

  // Fluid phase species names
  extern amrex::Vector<std::string> species;

  // Species unique identifying code
  extern std::vector<int> species_id;

  // Total number of fluid species
  extern int nspecies;

  // Specified constant gas phase species molecular weight
  extern amrex::Vector<amrex::Real> MW_gk0;

  // Specified constant gas phase species diffusion coefficients
  extern std::vector<amrex::Real> D_gk0;

  // Flag to understand if fluid is a mixture of fluid species
  extern int is_a_mixture;

  // Specified constant gas phase species specific heat
  extern std::vector<amrex::Real> cp_gk0;

  // Name to later reference when building input regions
  extern std::string name;

  struct FLUID_t {

    amrex::Real volfrac;
    amrex::Real temperature;
    amrex::Real pressure;
    amrex::Vector<amrex::Real> velocity;

    amrex::Real pressure_defined; // Flag to check if pressure is defined

    amrex::Vector<SPECIES::SPECIES_t> species;

    FLUID_t() :
      volfrac(-1.0),
      temperature(0.0),
      pressure(-1.0),
      species(0)
    {}

  };

  void Initialize ();

}

#endif
