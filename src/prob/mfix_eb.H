#ifndef MFIX_EB_PARMS_H_
#define MFIX_EB_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>

#include <AMReX_Vector.H>
#include <AMReX_LO_BCTYPES.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>

class MFIXEmbeddedBoundaries
{
  public:
    // Constructor
    MFIXEmbeddedBoundaries ()
      : m_fix_temperature(0)
      , m_has_flow(0)
    {}

    int fix_temperature () const { return m_fix_temperature; }
    void set_fix_temperature (const int val) { m_fix_temperature = val; }

    int has_flow () const { return m_has_flow; }
    void set_has_flow (const int val) { m_has_flow = val; }

  private:
    // Flag to specify if we have Dirichlet BCs on the EB
    int m_fix_temperature;

    // Flag to specify that we have Dirichlet flow though EB face
    int m_has_flow;
};


struct EB_t {

  EB_t ()
    : temperature(0.)
    , has_normal(0)
    , normal(0.)
    , normal_tol(0.)
  {}

  // Temperature to set on EB walls
  // e.g., it can be used to set inhomogeneous Dirichlet BCs for fluid
  // temperature solver
  amrex::Real temperature;

  // Normal to exclude EB faces. By default all faces within a specified
  // region will
  int has_normal;
  amrex::Vector<amrex::Real> normal;
  amrex::Real normal_tol;
};

#endif
