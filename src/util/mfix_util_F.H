#ifndef _MFIX_UTIL_F_H_
#define _MFIX_UTIL_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef BL_SINGLE_PRECISION_PARTICLES
typedef float amrex_particle_real;
#else
typedef double amrex_particle_real;
#endif

#ifdef __cplusplus
extern "C"
{
#endif

  void phifill (
      BL_FORT_FAB_ARG_3D(state),
      const int* dlo, const int* dhi,
      const amrex_real* dx, const amrex_real* glo,
      const amrex_real* time, const int* bc);

#ifdef __cplusplus
}
#endif

#endif
