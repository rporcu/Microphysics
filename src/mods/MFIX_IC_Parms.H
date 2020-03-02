#ifndef MFIX_IC_PARMS_H_
#define MFIX_IC_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <AMReX_Vector.H>
#include <MFIX_FLUID_Parms.H>
#include <MFIX_DEM_Parms.H>

namespace IC
{

  struct IC_t {

    const amrex::RealBox* region;

    FLUID::FLUID_t fluid;

    amrex::Vector<DEM::DEM_t> solids;
    std::string packing;

  };

  extern amrex::Vector<IC_t> ic;

  void Initialize ();

}

#endif
