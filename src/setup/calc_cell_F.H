#ifndef __CALC_CELL_F_H_
#define __CALC_CELL_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

void calc_cell_ic(const amrex::Real& dx, const amrex::Real& dy, const amrex::Real& dz,
                  const amrex::Real* ic_lo, const amrex::Real* ic_hi,
                  int& i_w, int& i_e, int& j_s,
                  int& j_n, int& k_b, int& k_t);

#ifdef __cplusplus
}
#endif

#endif
