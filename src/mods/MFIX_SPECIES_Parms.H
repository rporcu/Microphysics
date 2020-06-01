#ifndef MFIX_SPECIES_PARMS_H_
#define MFIX_SPECIES_PARMS_H_

#include <AMReX_REAL.H>

namespace SPECIES
{

  enum SPECIESDIFFUSIVITYMODEL {ConstantSpeciesDiffusivity};

  // Total number of species
  extern int nspecies;

  // Species names
  extern std::vector<std::string> species;

  // Species diffusion coefficients
  extern std::vector<amrex::Real> D_0;

  struct SPECIES_t {

    std::vector<amrex::Real> mass_fractions;

    SPECIES_t() {}

    SPECIES_t(const int nspecies) :
      mass_fractions(nspecies, 0.0)
    {}

  };

  void Initialize ();

}

#endif
