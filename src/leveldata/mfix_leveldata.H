#ifndef _MFIX_LEVELDATA_H_
#define _MFIX_LEVELDATA_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabFactory.H>
#include <AMReX_MultiFab.H>

#include <mfix_fluid_parms.H>

struct Transfer {
  enum {
    velx = 0,
    vely,
    velz,
    beta,
    gammaTp,
    gamma,
    count
  };
};

struct ChemTransfer
{
  ChemTransfer (const int nspecies_g, const int nreactions)
    : ro_gk_txfr(0)
    , vel_g_txfr(ro_gk_txfr + nspecies_g*(nreactions>0))
    , h_g_txfr(vel_g_txfr + 3*(nreactions>0))
    , count (h_g_txfr + 1*(nreactions>0))
  {}

  const int ro_gk_txfr;
  const int vel_g_txfr;
  const int h_g_txfr;
  const int count;
};

class LevelData
{
  public:
    // Default Constructor
    LevelData ()
      : level_allocated(0)
      , solve_enthalpy(0)
      , solve_species(0)
      , nspecies_g(0)
      , solve_reactions(0)
      , nreactions(0)
    {}

    // Constructor
    LevelData (amrex::BoxArray const& ba,
               amrex::DistributionMapping const& dm,
               const int nghost,
               amrex::FabFactory<amrex::FArrayBox> const& fact,
               const int solve_enthalpy,
               const int solve_species,
               const int nspecies_g,
               const int solve_reactions,
               const int nreactions);

    // Reset initial values
    void resetValues (const amrex::Real covered_val);

    // Destructor
    ~LevelData ();

    // Cell-centered MultiFabs
    // Void fraction
    amrex::MultiFab* ep_g;

    // Gas perturbational pressure
    amrex::MultiFab* p_g;
    amrex::MultiFab* p_go;

    // Gas density
    amrex::MultiFab* ro_g;
    amrex::MultiFab* ro_go;

    // Tracer in gas
    amrex::MultiFab* trac;
    amrex::MultiFab* trac_o;

    // Gas velocity
    amrex::MultiFab* vel_g;
    amrex::MultiFab* vel_go;

    // Base state pressure
    amrex::MultiFab* p0_g;

    // Pressure gradients
    amrex::MultiFab* gp;

    // Gas thermodynamic pressure
    amrex::Real pressure_g;
    amrex::Real pressure_go;

    // Gas temperature
    amrex::MultiFab* T_g;
    amrex::MultiFab* T_go;
    amrex::MultiFab* T_g_on_eb;

    // Gas enthalpy
    amrex::MultiFab* h_g;
    amrex::MultiFab* h_go;

    // Gas species mass fractions
    amrex::MultiFab* X_gk;
    amrex::MultiFab* X_gko;

    // Cell-based
    amrex::MultiFab* vort;

    // interphase transfer (drag, heat, etc)
    amrex::MultiFab* txfr;

    // Interphase transfer due to heterogeneous reactions
    amrex::MultiFab* chem_txfr;

    // div (ep_g * u)
    amrex::MultiFab* diveu;

    // Solution for MAC projection
    amrex::MultiFab* mac_phi;

    amrex::MultiFab* divtau_o;

    const int level_allocated;

    const int solve_enthalpy;
    const int solve_species;
    const int nspecies_g;
    const int solve_reactions;
    const int nreactions;
};

#endif
