#ifndef _MFIX_SET_GRADP_BCS_HPP_
#define _MFIX_SET_GRADP_BCS_HPP_

#include <mfix_F.H>
#include <mfix.H>
#include <mfix_des_F.H>
#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

void 
set_gradp_bcs (const Box& bx,
               FArrayBox& gp_fab,
               IArrayBox& bct_ilo_fab, 
               IArrayBox& bct_ihi_fab,
               IArrayBox& bct_jlo_fab,
               IArrayBox& bct_jhi_fab,
               IArrayBox& bct_klo_fab,
               IArrayBox& bct_khi_fab,
               Box& domain,
               BcList& bc_list,
               const int* ng);

#endif
