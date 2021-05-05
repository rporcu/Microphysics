#ifndef MFIX_PIC_PARMS_H_
#define MFIX_PIC_PARMS_H_

#include <AMReX_REAL.H>

namespace PIC
{
  extern int solve;

  // Coefficients for computing the solid stress
  extern amrex::Real Ps;
  extern amrex::Real beta;
  extern amrex::Real ep_cp;

  extern amrex::Real velfac;
  extern amrex::Real velfac_alpha;
  extern amrex::Real velfac_intercept;

  extern amrex::Real small_number;
  extern amrex::Real damping_factor;

  extern amrex::Real sigmoidal_offset;

  extern amrex::Real advance_eps_step;

  void Initialize ();
}

#endif
