#ifndef _MFIX_DIVOP_CONV_HPP_
#define _MFIX_DIVOP_CONV_HPP_

#include <mfix.H>


void 
mfix_apply_eb_redistribution ( Box& bx,
                               MultiFab& conv,
                               MultiFab& divc,
                               MultiFab& ep_g,
                               MFIter* mfi,
                               const int conv_comp,
                               const int ncomp,
                               const EBCellFlagFab& flags_fab,
                               const MultiFab* volfrac,
                               Box& domain,
                               const int cyclic_x,
                               const int cyclic_y,
                               const int cyclic_z,
                               const Real* dx);


#endif
