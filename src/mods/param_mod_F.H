#ifdef __cplusplus
extern "C"
{
#endif

#ifndef __PARAM_MOD_F_H_
#define __PARAM_MOD_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

int get_dim_ic();
amrex::Real get_my_huge();

bool is_defined_db(const amrex::Real& x);
bool is_defined_i(const int& x);

bool is_undefined_db(const amrex::Real& x);
bool is_undefined_i(const int& x);

#endif

#ifdef __cplusplus
}
#endif
