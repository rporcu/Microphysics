#ifndef MFIX_EB_PARMS_H_
#define MFIX_EB_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>

#include <AMReX_Vector.H>
#include <AMReX_LO_BCTYPES.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>

namespace EB
{

  struct EB_t {

    // Temperature to set on EB walls
    // e.g., it can be used to set inhomogeneous Dirichlet BCs for fluid
    // temperature solver
    amrex::Real temperature;

    EB_t ()
      : temperature(0.)
    {}
  };

}

#endif
