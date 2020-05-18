#ifndef _MFIX_F_H_
#define _MFIX_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

void init_collision (const int* mmax,
                     amrex::Real *min_dp,
                     amrex::Real *min_roi,
                     amrex::Real *max_dp,
                     amrex::Real *max_roi,
                     amrex::Real *sum_dp,
                     amrex::Real *sum_roi,
                     const amrex::Real& tcoll_ratio,
                     amrex::Real* etan,
                     amrex::Real* etan_w,
                     amrex::Real* etat,
                     amrex::Real* etat_w,
                     amrex::Real* neighborhood,
                     amrex::Real* dtsolid);

void mfix_usr3 (amrex::Real* vel_g, const int* ulo, const int* uhi,
                amrex::Real* p_g,   const int* slo, const int* shi,
                amrex::Real* dx, amrex::Real* dy, amrex::Real* dz);

void mfix_collect_fluid (const int* lo, const int* hi,
                         const int* domlo, const int* domhi,
                         amrex::Real* ep_g, const int* slo, const int* shi,
                         amrex::Real* p_g, const int* plo, const int* phi,
                         amrex::Real* vel_g, const int* vlo, const int* vhi,
                         const amrex::Real* vratio, const int* rlo, const int* rhi,
                         const amrex::Real* x_w,
                         const amrex::Real* x_e,
                         const amrex::Real* y_s,
                         const amrex::Real* y_n,
                         const amrex::Real* z_b,
                         const amrex::Real* z_t,
                         const amrex::Real* dx,
                         amrex::Real* sum_ep_g,
                         amrex::Real* sum_p_g,
                         amrex::Real* sum_vol,
                         amrex::Real* sum_velx, amrex::Real* sum_vely, amrex::Real* sum_velz);

#ifdef __cplusplus
}
#endif

#endif
