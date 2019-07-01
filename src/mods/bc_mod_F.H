#ifndef __BC_MOD_F_H_
#define __BC_MOD_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

amrex::Real get_bc_u_g(const int& pID);
amrex::Real get_bc_v_g(const int& pID);
amrex::Real get_bc_w_g(const int& pID);

amrex::Real get_bc_t_g(const int& pID);

amrex::Real get_bc_ep_g(const int& pID);

int get_minf();
int get_pinf();
int get_pout();

void get_domain_bc (int* domain_bc_out);

#ifdef __cplusplus
}
#endif

#endif
