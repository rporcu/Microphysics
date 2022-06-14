#ifndef MFIX_IC_PARMS_H_
#define MFIX_IC_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <AMReX_Vector.H>

#include <mfix_regions_parms.H>
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

    const SOLIDS_t* get_solid (const int phase) const
    {
      for (int n(0); n < solids.size(); ++n)
        if (phase == solids[n].phase)
          return &(solids[n]);

      amrex::Abort("Error: solid type not found");
      return nullptr;
    }
  };

  extern amrex::Vector<IC_t> ic;

  void Initialize (const Regions& regions,
                   FluidPhase& fluid,
                   const SolidsPhase& solids);

}

#endif
