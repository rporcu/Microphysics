#ifndef MFIX_BC_PARMS_H_
#define MFIX_BC_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>

#include <AMReX_Vector.H>
#include <AMReX_LO_BCTYPES.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>

#include <MFIX_FLUID_Parms.H>
#include <MFIX_DEM_Parms.H>

namespace BC
{

  struct BC_t {

    // BC type (mi, pi, po, nsw, eb)
    int type;

    const amrex::RealBox* region;

    FLUID::FLUID_t fluid;

    amrex::Vector<DEM::DEM_t> solids;

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

  void Initialize (amrex::Geometry& geom);

}

#endif
