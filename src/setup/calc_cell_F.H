#ifndef __CALC_CELL_F_H_
#define __CALC_CELL_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

void calc_cell_ic(const Real& dx, const Real& dy, const Real& dz,
                  const Real& x_w, const Real& y_s, const Real& z_b,
                  const Real& x_e, const Real& y_n, const Real& z_t,
                  int& i_w, int& i_e, int& j_s,
                  int& j_n, int& k_b, int& k_t);

#ifdef __cplusplus
}
#endif

#endif
