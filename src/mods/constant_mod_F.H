#ifndef __CONSTANTS_MOD_F_H_
#define __CONSTANTS_MOD_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

void get_gravity (amrex::Real* grav);

#ifdef __cplusplus
}
#endif

#endif
