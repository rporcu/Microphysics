#ifndef _MFIX_LEVELDATA_H_
#define _MFIX_LEVELDATA_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabFactory.H>
#include <AMReX_MultiFab.H>

class LevelData
{
  public:
    // Default Constructor
    LevelData () = default;

    // Constructor
    LevelData (amrex::BoxArray const& ba,
               amrex::DistributionMapping const& dm,
               const int nghost,
               amrex::FabFactory<amrex::FArrayBox> const& fact);

    // Reset initial values
    void resetValues (const amrex::Real covered_val);

    // Destructor
    ~LevelData ();

    // Cell-centered MultiFabs
    // Void fraction
    amrex::MultiFab* ep_g;
    amrex::MultiFab* ep_go;

    // Gas pressure fraction
    amrex::MultiFab* p_g;
    amrex::MultiFab* p_go;

    // Gas density
    amrex::MultiFab* ro_g;
    amrex::MultiFab* ro_go;

    // Gas molecular weight
    amrex::MultiFab* MW_g;

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

    // Molecular viscosity
    amrex::MultiFab* mu_g;

    // Gas temperature
    amrex::MultiFab* T_g;
    amrex::MultiFab* T_go;
    amrex::MultiFab* T_g_on_eb;

    // Specific heat
    amrex::MultiFab* cp_g;

    // Thermal conductivity
    amrex::MultiFab* k_g;
    amrex::MultiFab* k_g_on_eb;

    // Gas enthalpy
    amrex::MultiFab* h_g;
    amrex::MultiFab* h_go;

    // Gas species
    amrex::MultiFab* X_gk;
    amrex::MultiFab* X_gko;

    // Species diffusion coefficients
    amrex::MultiFab* D_gk;

    // Gas species enthalpy
    amrex::MultiFab* cp_gk;

    // Gas species enthalpy
    amrex::MultiFab* h_gk;

    // Cell-based
    amrex::MultiFab* vort;

    // interphase transfer (drag, heat, etc)
    amrex::MultiFab* txfr;

    // These are multi-component multifabs
    amrex::MultiFab* xslopes_u;
    amrex::MultiFab* yslopes_u;
    amrex::MultiFab* zslopes_u;
    amrex::MultiFab* xslopes_s;
    amrex::MultiFab* yslopes_s;
    amrex::MultiFab* zslopes_s;
    amrex::MultiFab* xslopes_X;
    amrex::MultiFab* yslopes_X;
    amrex::MultiFab* zslopes_X;

    // div (ep_g * u)
    amrex::MultiFab* diveu;

    // Solution for MAC projection
    amrex::MultiFab* mac_phi;

    // RHS for MAC solve
    amrex::MultiFab* mac_rhs;

    // MAC velocities
    amrex::MultiFab* u_mac;
    amrex::MultiFab* v_mac;
    amrex::MultiFab* w_mac;
};

#endif
