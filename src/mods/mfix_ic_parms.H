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
    IC_t (const FluidPhase& fluid_phase) :
      fluid(fluid_phase.constraint_type, fluid_phase.nspecies, fluid_phase.parameters)
    {}

    const amrex::RealBox* region;

    FluidPhase::FLUID_t fluid;

    amrex::Vector<SolidsPhase::SOLIDS_t> solids;
    std::string packing;

  };

  extern amrex::Vector<IC_t> ic;

  void Initialize (const FluidPhase& fluid,
                   const SolidsPhase& solids);

}

#endif
