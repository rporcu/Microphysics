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

  // Flag to specify if we have Dirichlet BCs on the EB
  extern int fix_temperature;

  // Flag to specify that we have Dirichlet flow though EB face
  extern int has_flow;

  // Flag to compute the area of the EB region
  extern int compute_area;

  struct EB_t {

    // Temperature to set on EB walls
    // e.g., it can be used to set inhomogeneous Dirichlet BCs for fluid
    // temperature solver
    amrex::Real temperature;

    // BC area -- used to convert volumetric flows to velocites
    amrex::Real area;

    // Normal to exclude EB faces. By default all faces within a specified
    // region will
    int has_normal;
    amrex::Vector<amrex::Real> normal;
    amrex::Real normal_tol;

    EB_t ()
      : temperature(0.)
      , area(-1.)
      , has_normal(0)
      , normal(0.)
      , normal_tol(0.)
    {}
  };

}

#endif
