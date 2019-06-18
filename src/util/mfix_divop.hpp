#ifndef _MFIX_DIVOP_HPP_
#define _MFIX_DIVOP_HPP_

#include <mfix_F.H>
#include <mfix.H>
#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

void 
compute_divop(Box& bx,
              MultiFab& divergence,
              MultiFab& velocity,
              MultiFab& ep_g,
              MFIter* mfi,
              FArrayBox& fx,
              FArrayBox& fy,
              FArrayBox& fz,
              Array<const MultiCutFab*, AMREX_SPACEDIM>& areafrac,
              Array<const MultiCutFab*, AMREX_SPACEDIM>& facecent,
              const EBCellFlagFab& flags_fab,
              const MultiFab* volfrac,
              const MultiCutFab* bndrycent_fab,
              const int cyclic_x,
              const int cyclic_y,
              const int cyclic_z,
              Box& domain,
              const Real* dx,
              const int* nghost,
              const Array4<Real>* mu = nullptr,
              const int* do_explicit_diffusion = nullptr);

#endif
