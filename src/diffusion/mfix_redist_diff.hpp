#ifndef _MFIX_DIVOP_DIFF_HPP_
#define _MFIX_DIVOP_DIFF_HPP_

#include <mfix_F.H>
#include <mfix.H>
#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

void 
compute_redist_diff(
               Box& bx,
               MultiFab& divtau,
               MultiFab& ep_g,
               MultiFab& divtau_aux,
               MFIter* mfi,
               const EBCellFlagFab& flags_fab,
               const MultiFab* volfrac,
               const int cyclic_x,
               const int cyclic_y,
               const int cyclic_z,
               Box& domain,
               const Real* dx);
#endif
