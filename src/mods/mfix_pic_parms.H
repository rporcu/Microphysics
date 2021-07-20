#ifndef MFIX_PIC_PARMS_H_
#define MFIX_PIC_PARMS_H_

#include <AMReX_REAL.H>

#define MFP_DISABLED 0
#define SCALE_MFP_DISABLED 0
#define SCALE_MFP_CHAPMAN  0
#define SCALE_MFP_CUTCHIS  0
#define SCALE_MFP_MUSSER   1

namespace PIC
{
  extern int solve;

  extern int verbose;

  // Coefficients for computing the solid stress
  extern amrex::Real Ps;
  extern amrex::Real beta;
  extern amrex::Real ep_cp;

  extern amrex::Real vel_ref_frame;

  extern amrex::Real small_number;
  extern amrex::Real damping_factor;

  extern amrex::Real damping_factor_wall_normal;
  extern amrex::Real damping_factor_wall_tangent;

  extern amrex::Real advance_vel_p;

  extern int max_iter;

  enum struct InitialStepType {
    Invalid, zero_eps, nth_eps, taylor_approx
  };

  extern InitialStepType initial_step;

  void Initialize ();
}

#endif
