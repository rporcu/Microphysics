#ifndef __RANDOM_NB_MOD_F_H_
#define __RANDOM_NB_MOD_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

amrex::Real get_random();

#ifdef __cplusplus
}
#endif

#endif
