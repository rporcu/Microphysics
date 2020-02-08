#ifndef _MFIX_EB_F_H_
#define _MFIX_EB_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

void mfix_get_walls (int* bcv,
                     int* exists,
                     void* normal,
                     void* center,
                     int* po_noParOut);

void mfix_get_real_walls(int* bcv,
                         int* exists,
                         void* normal,
                         void* center);

void ls_has_walls(int * has_wall,
                  const amrex::Real* phi,
                  const int* phlo,
                  const int* phhi,
                  const amrex::Real* tol);

void calc_wall_collisions(void* rparticles,
                          const int* ntot,
                          const int* nrp,
                          amrex::Real* tow,
                          amrex::Real* fc,
                          amrex::Real* subdt,
                          const int* valid,
                          const int* vlo,
                          const int* vhi,
                          const amrex::Real* phi, const int* phlo, const int* phhi,
                          const amrex::Real* dx,
                          const int* n_refine);

#ifdef __cplusplus
}
#endif

#endif
