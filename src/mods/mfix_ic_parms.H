#ifndef MFIX_IC_PARMS_H_
#define MFIX_IC_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <AMReX_Vector.H>
#include <mfix_fluid_parms.H>
#include <mfix_solids_parms.H>

namespace IC
{

  struct IC_t {

    // Constructor
    IC_t (FluidPhase& fluid_in) :
      fluid(fluid_in)
    {}

    const amrex::RealBox* region;

    FLUID_t fluid;

    amrex::Vector<SOLIDS_t> solids;
    std::string packing;

  };

  extern amrex::Vector<IC_t> ic;

  void Initialize (FluidPhase& fluid,
                   const SolidsPhase& solids);

}

#endif
