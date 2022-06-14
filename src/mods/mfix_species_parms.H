#ifndef MFIX_SPECIES_PARMS_H_
#define MFIX_SPECIES_PARMS_H_

#include <AMReX_REAL.H>

struct DIFFUSIVITYMODEL { enum {Invalid=-1, Constant}; };
struct SPECIFICHEATMODEL { enum {Invalid=-1, Constant, NASA7Polynomials}; };

class Species
{
  public:
    // Constructor
    Species();

    static constexpr int NMAX = 15;

    // List of recognized species diffusivity models
    int DiffusivityModel;

    // List of recognized species specific heat models
    int SpecificHeatModel;

    // Flag to solve species equations
    int solve;

    // Total number of species
    int nspecies;

    // Specified species names
    amrex::Vector<std::string> names;

    // Species unique identifying code (at the moment = their index in the input
    // entries)
    amrex::Vector<int> IDs;

    // Specified species constant molecular weights
    amrex::Vector<amrex::Real> MW_k0;

    // Specified species constant diffusion coefficients
    amrex::Real D_0;

    // Specified parameters species specific heat
    amrex::Vector<amrex::Real> cp_k0;

    // Enthalpy of formation
    amrex::Vector<amrex::Real> H_fk0;

    int is_initialized;

    void Initialize ();
};


struct SPECIES_t {
  // Constructor
  SPECIES_t() :
    mass_fraction(0.)
  {}

  // Species mass fraction, used for ICs and BCs
  amrex::Real mass_fraction;
};

#endif
