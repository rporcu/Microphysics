#ifndef _MFIX_DIVOP_F_H_
#define _MFIX_DIVOP_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

void compute_divop (
  const int* lo, const int* hi,
  amrex::Real* div, const int* divlo, const int* divhi, 
  const amrex::Real* vel   ,  const int* vello, const int* velhi,
  const amrex::Real* fx   ,  const int* fxlo, const int* fxhi,
  const amrex::Real* fy   ,  const int* fylo, const int* fyhi,
  const amrex::Real* fz   ,  const int* fzlo, const int* fzhi,
  const amrex::Real* ep_g   ,  const int* elo, const int* ehi,
  const amrex::Real* afrac_x,const int*  axlo,const int*  axhi, 
  const amrex::Real* afrac_y,const int*  aylo,const int*  ayhi, 
  const amrex::Real* afrac_z,const int*  azlo,const int*  azhi, 
  const amrex::Real* cent_x,const int*   cxlo,const int*  cxhi, 
  const amrex::Real* cent_y,const int*   cylo,const int*  cyhi, 
  const amrex::Real* cent_z,const int*   czlo,const int*  czhi, 
  const void* flags, const int*  flo, const int*   fhi, 
  const amrex::Real* vfrac, const int* vflo, const int*  vfhi,
  const amrex::Real* bcent, const int*  blo, const int*   bhi,
  const int* domlo, const int* domhi,
  const amrex::Real* dx, const int* ng,
  const amrex::Real* mu=nullptr,
  const int* do_explicit_diffusion=nullptr
);

#ifdef __cplusplus
}
#endif

#endif
