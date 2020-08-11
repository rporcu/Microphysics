#ifndef MFIX_PIC_PARMS_H_
#define MFIX_PIC_PARMS_H_

#include <AMReX_REAL.H>

namespace PIC
{

//  enum COLLISIONMODEL {LSD};
//
//  extern COLLISIONMODEL CollisionModel;

  extern int solve;

  // Coefficients for computing the solid stress
  extern amrex::Real Ps;
  extern amrex::Real beta;
  extern amrex::Real ep_cp;

  extern amrex::Real velfac;

  extern amrex::Real small_number;
  extern amrex::Real damping_factor;

  // TODO check on this
  //struct PIC_t {

  //  amrex::Real volfrac;
  //  amrex::Real statwt;
  //  amrex::Real temperature;
  //  amrex::Vector<amrex::Real> velocity;

  //  SPECIES::SPECIES_t species;

  //  INPUT_DIST_t diameter;
  //  INPUT_DIST_t density;
  //};

  void Initialize ();
}

#endif
