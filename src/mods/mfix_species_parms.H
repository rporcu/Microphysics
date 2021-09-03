#ifndef MFIX_SPECIES_PARMS_H_
#define MFIX_SPECIES_PARMS_H_

#include <AMReX_REAL.H>

namespace SPECIES
{
  constexpr int NMAX = 15;

  // List of recognized species diffusivity models
  struct DIFFUSIVITYMODEL { enum {Invalid=-1, Constant}; };

  extern int DiffusivityModel;

  // List of recognized species specific heat models
  struct SPECIFICHEATMODEL { enum {Invalid=-1, Constant, NASA7Polynomials}; };

  extern int SpecificHeatModel;

  // Flag to solve species equations
  extern int solve;

  // Total number of species
  extern int nspecies;

  // Specified species names
  extern amrex::Vector<std::string> species;

  // Species unique identifying code (at the moment = their index in the input
  // entries)
  extern amrex::Vector<int> species_id;

  // Specified species constant molecular weights
  extern amrex::Vector<amrex::Real> MW_k0;

  // Specified species constant diffusion coefficients
  extern amrex::Vector<amrex::Real> D_k0;

  // Flag to solve enthalpy species equation
  extern int solve_enthalpy;

  // Specified parameters species specific heat
  extern amrex::Vector<amrex::Real> cp_k0;

  // Enthalpy of formation
  extern amrex::Vector<amrex::Real> H_fk0;

  struct SPECIES_t {
    // Constructor
    SPECIES_t() :
      mass_fraction(0.)
    {}

    // Species mass fraction, used for ICs and BCs
    amrex::Real mass_fraction;
  };

  void Initialize ();

}

#endif
