#ifndef _MFIX_SPECIES_H_
#define _MFIX_SPECIES_H_

#include <AMReX_REAL.H>

class MFIXSpecies
{
  public:
    // 
    struct DiffusivityModel { enum {Invalid=-1, Constant}; };

    //
    struct SpecificHeatModel { enum {Invalid=-1, Constant, NASA7Polynomials}; };

    // Constructor
    MFIXSpecies();

    void Initialize ();

    // List of recognized species diffusivity models
    int diffusivity_model () const { return m_diffusivity_model; };

    // List of recognized species specific heat models
    int specific_heat_model () const { return m_specific_heat_model; };

    // Flag to solve species equations
    int solve () const { return m_solve; };

    // Total number of species
    int nspecies () const { return m_nspecies; };

    // Specified species names
    const amrex::Vector<std::string>& names () const { return m_names; };
    const std::string& names (const int k) const { return m_names[k]; };

    // Species unique identifying code (at the moment = their index in the input
    // entries)
    const amrex::Vector<int>& IDs () const { return m_IDs; };
    int IDs (const int k) const { return m_IDs[k]; };

    // Specified species constant molecular weights
    const amrex::Vector<amrex::Real>& MW_k () const { return m_MW_k; };
    amrex::Real MW_k (const int k) const { return m_MW_k[k]; };

    // Specified species constant diffusion coefficients
    amrex::Real diffusivity () const { return m_D; };

    // Specified parameters species specific heat
    const amrex::Vector<amrex::Real>& cp_k () const { return m_cp_k; };
    amrex::Real cp_k (const int k) const { return m_cp_k[k]; };

    // Enthalpy of formation
    const amrex::Vector<amrex::Real>& H_fk () const { return m_H_fk; };
    amrex::Real H_fk (const int k) const { return m_H_fk[k]; };

    int isInitialized () const { return m_is_initialized; };

  private:
    // List of recognized species diffusivity models
    int m_diffusivity_model;

    // List of recognized species specific heat models
    int m_specific_heat_model;

    // Flag to solve species equations
    int m_solve;

    // Total number of species
    int m_nspecies;

    // Specified species names
    amrex::Vector<std::string> m_names;

    // Species unique identifying code (at the moment = their index in the input
    // entries)
    amrex::Vector<int> m_IDs;

    // Specified species constant molecular weights
    amrex::Vector<amrex::Real> m_MW_k;

    // Specified species constant diffusion coefficients
    amrex::Real m_D;

    // Specified parameters species specific heat
    amrex::Vector<amrex::Real> m_cp_k;

    // Enthalpy of formation
    amrex::Vector<amrex::Real> m_H_fk;

    int m_is_initialized;
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
