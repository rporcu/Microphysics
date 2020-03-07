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
               amrex::FabFactory<amrex::FArrayBox> const& fact,
               const amrex::Real covered_val);

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

    // Cell-based
    amrex::MultiFab* vort;

    amrex::MultiFab* drag;
    
    // These are multi-component multifabs
    amrex::MultiFab* xslopes_u;
    amrex::MultiFab* yslopes_u;
    amrex::MultiFab* zslopes_u;
    amrex::MultiFab* xslopes_s;
    amrex::MultiFab* yslopes_s;
    amrex::MultiFab* zslopes_s;

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

    // used if load_balance_type == "KnapSack"
    amrex::MultiFab* particle_cost;
    amrex::MultiFab* fluid_cost;
};

#endif
