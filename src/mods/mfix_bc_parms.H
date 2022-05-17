#ifndef MFIX_BC_PARMS_H_
#define MFIX_BC_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>

#include <AMReX_Vector.H>
#include <AMReX_LO_BCTYPES.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>

#include <mfix_bc_list.H>
#include <mfix_eb_parms.H>
#include <mfix_regions_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_solids_parms.H>

#include <bits/stdc++.h>

namespace BC
{

  struct BC_t {
    // Constructor
    BC_t (FluidPhase& fluid_in) :
      fluid(fluid_in)
    {}

    // BC type (mi, pi, po, nsw, eb)
    int type;

    const amrex::RealBox* region;

    EB::EB_t eb;

    FLUID_t fluid;

    amrex::Vector<SOLIDS_t> solids;

  };

  // Lists of BCs applied to the domain extent
  extern amrex::Vector<int> bc_xlo, bc_xhi;
  extern amrex::Vector<int> bc_ylo, bc_yhi;
  extern amrex::Vector<int> bc_zlo, bc_zhi;

  // Flag to solve fluid equations
  extern int delp_dir;
  extern amrex::Real delp[3];

  // Flags for DEM domain extent reinforcement.
  extern int domain_bc[6];

  extern std::bitset<6> flow_plane;

  // Vectors storing EB planes for level-set creation
  extern amrex::Vector<amrex::EB2::PlaneIF> flow_planes;
  extern amrex::Vector<amrex::EB2::PlaneIF> wall_planes;

  // Use AMReX LinOpBCType to define the domain extents
  extern std::array<amrex::LinOpBCType,3> ppe_lobc;
  extern std::array<amrex::LinOpBCType,3> ppe_hibc;

  extern std::array<amrex::LinOpBCType,3> diff_vel_lobc;
  extern std::array<amrex::LinOpBCType,3> diff_vel_hibc;

  extern std::array<amrex::LinOpBCType,3> diff_scal_lobc;
  extern std::array<amrex::LinOpBCType,3> diff_scal_hibc;

  extern std::array<amrex::LinOpBCType,3> diff_temperature_lobc;
  extern std::array<amrex::LinOpBCType,3> diff_temperature_hibc;

  extern std::array<amrex::LinOpBCType,3> diff_species_lobc;
  extern std::array<amrex::LinOpBCType,3> diff_species_hibc;

  // Data structure for BC data
  extern amrex::Vector<BC_t> bc;

  void Initialize (amrex::Geometry& geom,
                   const Regions& regions,
                   FluidPhase& fluid,
                   const SolidsPhase& solids);

  void read_bc_density (amrex::ParmParse pp, FLUID_t* fluid);
  void read_bc_velocity (amrex::ParmParse pp, FLUID_t* fluid);
  void read_bc_volflow (amrex::ParmParse pp, FLUID_t* fluid);
  void read_bc_temperature (amrex::ParmParse pp, FLUID_t* fluid);
  void read_bc_species (amrex::ParmParse pp,
                        const amrex::Vector<std::string> species,
                        FLUID_t* fluid);

}

#endif
