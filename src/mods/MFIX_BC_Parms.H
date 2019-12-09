#ifndef MFIX_BC_PARMS_H_
#define MFIX_BC_PARMS_H_

#include <AMReX_REAL.H>

namespace BC
{

  // enum DENSITYMODEL {ConstantDensity, IdealGas};

  // Flag to solve fluid equations
  extern int delp_dir;

  extern amrex::Real delp[3];

  void Initialize ();

}

#endif
