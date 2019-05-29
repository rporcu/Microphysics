#ifdef __cplusplus
extern "C"
{
#endif

#ifndef __FLD_CONSTANTS_MOD_F_H_
#define __FLD_CONSTANTS_MOD_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

amrex::Real get_ro_g0();
amrex::Real get_mu_g0();
amrex::Real get_mw_avg();

#endif

#ifdef __cplusplus
}
#endif
