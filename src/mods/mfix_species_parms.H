#ifndef MFIX_SPECIES_PARMS_H_
#define MFIX_SPECIES_PARMS_H_

#include <AMReX_REAL.H>

namespace SPECIES
{
  // List of recognized species diffusivity models
  struct DIFFUSIVITYMODEL { enum {Invalid=-1, Constant}; };

  extern int DiffusivityModel;

  // List of recognized species specific heat models
  struct SPECIFICHEATMODEL { enum {Invalid=-1, Constant}; };

  extern int SpecificHeatModel;

  // Flag to solve species equations
  extern int solve;

  // Total number of species
  extern int nspecies;

  // Specified species names
  extern std::vector<std::string> species;

  // Specified species constant molecular weights
  extern std::vector<amrex::Real> MW_k0;

  // Specified species constant diffusion coefficients
  extern std::vector<amrex::Real> D_k0;

  // Flag to solve enthalpy species equation
  extern int solve_enthalpy;

  // Specified constant species specific heat
  extern std::vector<amrex::Real> cp_k0;

  struct SPECIES_t {

    std::vector<amrex::Real> mass_fractions;

    SPECIES_t() {}

    SPECIES_t(const int nspec) :
      mass_fractions(nspec, 0.0)
    {}

  };

  void Initialize ();

}

#endif
