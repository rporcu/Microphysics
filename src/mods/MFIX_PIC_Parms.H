#ifndef MFIX_PIC_PARMS_H_
#define MFIX_PIC_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace PIC
{
  constexpr int NMAX = 10;
  extern AMREX_GPU_DEVICE_MANAGED int NPHASE;

//  enum COLLISIONMODEL {LSD};
//
//  extern AMREX_GPU_DEVICE_MANAGED COLLISIONMODEL CollisionModel;

  extern int solve;

  // Coefficients for computing the solid stress
  extern AMREX_GPU_DEVICE_MANAGED amrex::Real Ps;
  extern AMREX_GPU_DEVICE_MANAGED amrex::Real beta;
  extern AMREX_GPU_DEVICE_MANAGED amrex::Real ep_cp;

  extern AMREX_GPU_DEVICE_MANAGED amrex::Real velfac;

  extern AMREX_GPU_DEVICE_MANAGED amrex::Real small_number;
  extern AMREX_GPU_DEVICE_MANAGED amrex::Real damping_factor;

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
