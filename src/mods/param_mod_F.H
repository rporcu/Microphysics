#ifndef __PARAM_MOD_F_H_
#define __PARAM_MOD_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

int get_dim_ic();

int get_dim_bc();

amrex::Real get_my_huge();

_Bool is_defined_db(const amrex::Real& x);

_Bool is_undefined_db(const amrex::Real& x);

amrex::Real get_undefined();

#ifdef __cplusplus
}
#endif

AMREX_GPU_HOST_DEVICE
_Bool is_equal(const amrex::Real& x, const amrex::Real& y);

#endif
