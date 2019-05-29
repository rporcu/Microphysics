#ifdef __cplusplus
extern "C"
{
#endif

#ifndef _EOS_MOD_F_H_
#define _EOS_MOD_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

amrex::Real sutherland(const amrex::Real& tg);

#endif

#ifdef __cplusplus
}
#endif
