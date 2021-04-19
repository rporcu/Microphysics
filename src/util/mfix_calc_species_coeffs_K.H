#ifndef _MFIX_CALC_FLUID_COEFFS_K_H_
#define _MFIX_CALC_FLUID_COEFFS_K_H_

#include <AMReX_REAL.H>

namespace FLUID {

// Compute fluid enthalpy
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real calc_h_g (const amrex::Real cp_g,
                      const amrex::Real T_g,
                      const amrex::Real T_ref,
                      const amrex::Real H_f)
{
  return H_f + cp_g*(T_g - T_ref);
}

// Compute fluid enthalpy
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real calc_T_g (const amrex::Real cp_g,
                      const amrex::Real h_g,
                      const amrex::Real T_ref,
                      const amrex::Real H_f)
{
  return (h_g - H_f)/cp_g + T_ref;
}

}

namespace SOLIDS {

// Compute fluid enthalpy
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real calc_h_s (const amrex::Real cp_s,
                      const amrex::Real T_s,
                      const amrex::Real T_ref,
                      const amrex::Real H_f)
{
  return H_f + cp_s*(T_s - T_ref);
}

// Compute fluid enthalpy
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real calc_T_s (const amrex::Real cp_s,
                      const amrex::Real h_s,
                      const amrex::Real T_ref,
                      const amrex::Real H_f)
{
  return (h_s - H_f)/cp_s + T_ref;
}

}

#endif
