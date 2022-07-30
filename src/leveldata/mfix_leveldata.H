#ifndef _MFIX_LEVELDATA_H_
#define _MFIX_LEVELDATA_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabFactory.H>
#include <AMReX_MultiFab.H>

#include <mfix_eb.H>
#include <mfix_fluid.H>
#include <mfix_reactions.H>


struct Transfer {

  Transfer (const int nspecies_g, const int nreactions)
    : velx(0)
    , vely(1)
    , velz(2)
    , drag_coeff(3)
    , gammaTp(4)
    , convection_coeff(5)
    , countNoChem(convection_coeff+1) // Num txfr terms (no rxns)
    , ro_gk_txfr(countNoChem)
    , vel_g_txfr(ro_gk_txfr + nspecies_g*(nreactions>0))
    , h_g_txfr(vel_g_txfr + 3*(nreactions>0))
    , countChem ((nspecies_g + 3 + 1)*(nreactions>0)) // Num of chem rxns txfr terms
    , count (countNoChem + countChem)
  {}

  const int velx;
  const int vely;
  const int velz;
  const int drag_coeff;
  const int gammaTp;
  const int convection_coeff;
  const int countNoChem;
  const int ro_gk_txfr;
  const int vel_g_txfr;
  const int h_g_txfr;
  const int countChem;
  const int count;
};



class LevelData {

  public:
    // Default Constructor
    LevelData ()
      : level_allocated(0)
    {}

    // Constructor
    LevelData (amrex::BoxArray const& ba,
               amrex::DistributionMapping const& dm,
               const int nghost,
               amrex::FabFactory<amrex::FArrayBox> const& fact,
               MFIXEmbeddedBoundaries* embedded_boundaries,
               MFIXFluidPhase* fluid_in,
               MFIXReactions* reactions_in,
               bool only_ep_g = false);

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
    amrex::MultiFab* thermodynamic_p_g;
    amrex::MultiFab* thermodynamic_p_go;

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

    // div (ep_g * u)
    amrex::MultiFab* diveu;

    // Solution for MAC projection
    amrex::MultiFab* mac_phi;

    amrex::MultiFab* divtau_o;

    const int level_allocated;

    MFIXEmbeddedBoundaries* m_embedded_boundaries;
    MFIXFluidPhase* fluid;
    MFIXReactions* reactions;

    bool ep_g_only = false;
};

#endif
