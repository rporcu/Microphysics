#ifndef _MFIX_DES_F_H_
#define _MFIX_DES_F_H_

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

void set_lsd_collision_coefficients (amrex::Real* mew,
                                     amrex::Real* mew_w,
                                     amrex::Real* kn,
                                     amrex::Real* kn_w,
                                     amrex::Real* kt,
                                     amrex::Real* kt_w,
                                     amrex::Real* en,
                                     amrex::Real* en_w,
                                     amrex::Real* kt_fac,
                                     amrex::Real* kt_fac_w,
                                     amrex::Real* eta_fac,
                                     amrex::Real* eta_w_fac);

void calc_particle_collisions (void* rparticles, const int* nrp,
                               void* gparticles, const int* ngp,
                               const int* nbor_list, const int* size_nl,
                               amrex::Real* tow,
                               amrex::Real* fc,
                               amrex::Real* subdt,
                               int* ncoll);

void set_particle_properties (const int* pstate,
                              amrex::Real* pradius,
                              amrex::Real* pdensity,
                              amrex::Real* pvol,
                              amrex::Real* pmass,
                              amrex::Real* omoi,
                              amrex::Real* omega);

void particle_get_position (void* particles,
                            const int* nrp,
                            amrex::Real* x,
                            amrex::Real* y,
                            amrex::Real* z);

void des_init_time_loop (const  amrex::Real* tstart,
                         const amrex::Real* dt,
                         int* nsubsteps,
                         amrex::Real* subdt);

void call_usr3_des (const int* np,
                    void* particles);

void call_usr2_des (const int* np,
                    void* particles);

#ifdef __cplusplus
}
#endif

#endif
