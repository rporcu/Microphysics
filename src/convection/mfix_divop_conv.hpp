#ifndef _MFIX_DIVOP_CONV_HPP_
#define _MFIX_DIVOP_CONV_HPP_

#include <mfix.H>

void 
compute_divop_conv(
              Box& bx,
              MultiFab& divergence,
              MultiFab& ep_g,
              int conv_comp, int ncomp,
              MFIter* mfi,
              FArrayBox& fx,
              FArrayBox& fy,
              FArrayBox& fz,
              Array<const MultiCutFab*, AMREX_SPACEDIM>& areafrac,
              Array<const MultiCutFab*, AMREX_SPACEDIM>& facecent,
              const EBCellFlagFab& flags_fab,
              const MultiFab* volfrac,
              const MultiCutFab* bndrycent_fab,
              Box& domain,
              const int cyclic_x,
              const int cyclic_y,
              const int cyclic_z,
              const Real* dx);
#endif
