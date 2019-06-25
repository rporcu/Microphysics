#ifndef __IC_MOD_F_H_
#define __IC_MOD_F_H_

#ifdef __cplusplus
extern "C"
{
#endif

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

int ic_defined(const int& icv);

amrex::Real get_ic_u_g(const int& icv);

amrex::Real get_ic_v_g(const int& icv);

amrex::Real get_ic_w_g(const int& icv);

amrex::Real get_ic_x_w(const int& icv);

amrex::Real get_ic_y_s(const int& icv);

amrex::Real get_ic_z_b(const int& icv);

amrex::Real get_ic_x_e(const int& icv);

amrex::Real get_ic_y_n(const int& icv);

amrex::Real get_ic_z_t(const int& icv);

#ifdef __cplusplus
}
#endif

#endif
