#ifndef MFIX_BC_PARMS_H_
#define MFIX_BC_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <AMReX_Vector.H>
#include <MFIX_FLUID_Parms.H>
#include <MFIX_DEM_Parms.H>

namespace BC
{

  struct BC_t {

    std::string type;

    const amrex::RealBox* region;

    FLUID::FLUID_t fluid;

    amrex::Vector<DEM::DEM_t> solids;

  };

  extern amrex::Vector<BC_t> bc;

  // Flag to solve fluid equations
  extern int delp_dir;

  extern amrex::Real delp[3];

  void Initialize ();

}

#endif
