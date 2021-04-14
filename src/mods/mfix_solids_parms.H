#ifndef MFIX_SOLIDS_PARMS_H_
#define MFIX_SOLIDS_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <mfix_species_parms.H>

namespace SOLIDS
{
  constexpr int NMAX = 10;

  extern int NPHASE;

  // List of recognized fluid molecular weight models
  struct MOLECULARWEIGHTMODEL { enum {Invalid=-1, Constant, Mixture}; };

  extern int MolecularWeightModel;

  // List of recognized fluid specific heat models
  struct SPECIFICHEATMODEL { enum {Invalid=-1, Constant, NASA9Polynomials}; };

  extern int SpecificHeatModel;

  // List of recognized fluid thermal conductivity models
  struct THERMALCONDUCTIVITYMODEL{ enum {Invalid=-1, Constant}; };

  extern int ThermalConductivityModel;

  // Specified constant molecular weight of solids
  extern amrex::Real MW_s0;

  // Names of solids used in IC/BC setups
  extern amrex::Vector<std::string> names;

  // Flag to solve species solids equations
  extern int solve_species;

  // Species unique identifying code
  extern amrex::Vector<int> species_id;

  // Total number of dem species
  extern int nspecies;

  // Particle species
  extern amrex::Vector<std::string> species;

  // Specified constant solid phase species molecular weight
  extern amrex::Vector<amrex::Real> MW_sn0;

  // Specified constant specific heat
  extern amrex::Vector<amrex::Real> cp_p0;

  extern amrex::Real enthalpy_source;

  struct INPUT_DIST_t {
    std::string distribution;
    amrex::Real mean; // Mean value of distriubtion (or constant value)
    amrex::Real std;  // Standard deviation
    amrex::Real min;  // Minimum value of distribution
    amrex::Real max;  // Max value of distribution

    INPUT_DIST_t() :
      mean(-1.0),
      std(-1.0),
      min(-1.0),
      max(-1.0) { }

  };

  struct SOLIDS_t {

    std::string name;
    amrex::Real volfrac;
    amrex::Real statwt;
    amrex::Real temperature;
    amrex::Vector<amrex::Real> velocity;

    amrex::Vector<SPECIES::SPECIES_t> species;

    INPUT_DIST_t diameter;
    INPUT_DIST_t density;

    SOLIDS_t() :
      volfrac(0.0),
      statwt(0.0),
      temperature(0.0),
      velocity(0),
      species(0),
      diameter(),
      density()
    {}
  };

  void Initialize ();

}

#endif
