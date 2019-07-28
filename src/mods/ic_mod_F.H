#ifndef __IC_MOD_F_H_
#define __IC_MOD_F_H_

#ifdef __cplusplus
extern "C"
{
#endif

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

int ic_defined_cpp(const int& icv);

amrex::Real get_ic_p_g(const int& pID);
amrex::Real get_ic_ep_g(const int& pID);
amrex::Real get_ic_ep_s(const int& pID, const int& solid_type);
void get_ic_pack_type(const int& pID, char c_string[]);
void get_ic_dp_dist(const int& pID, const int& solid_type, char c_string[]);
amrex::Real get_ic_dp_mean(const int& pID, const int& solid_type);
amrex::Real get_ic_dp_std(const int& pID, const int& solid_type);
amrex::Real get_ic_dp_min(const int& pID, const int& solid_type);
amrex::Real get_ic_dp_max(const int& pID, const int& solid_type);
void get_ic_ro_s_dist(const int& pID, const int& solid_type, char c_string[]);
amrex::Real get_ic_ro_s_mean(const int& pID, const int& solid_type);
amrex::Real get_ic_ro_s_std(const int& pID, const int& solid_type);
amrex::Real get_ic_ro_s_min(const int& pID, const int& solid_type);
amrex::Real get_ic_ro_s_max(const int& pID, const int& solid_type);

amrex::Real get_ic_u_g(const int& pID);
amrex::Real get_ic_v_g(const int& pID);
amrex::Real get_ic_w_g(const int& pID);
amrex::Real get_ic_x_w(const int& pID);
amrex::Real get_ic_y_s(const int& pID);
amrex::Real get_ic_z_b(const int& pID);
amrex::Real get_ic_x_e(const int& pID);
amrex::Real get_ic_y_n(const int& pID);
amrex::Real get_ic_z_t(const int& pID);

amrex::Real get_ic_u_s(const int& pID, const int& solid_type);
amrex::Real get_ic_v_s(const int& pID, const int& solid_type);
amrex::Real get_ic_w_s(const int& pID, const int& solid_type);

#ifdef __cplusplus
}
#endif

#endif
