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
               Array4<Real> const& gp,
               Array4<int> const& bct_ilo, 
               Array4<int> const& bct_ihi,
               Array4<int> const& bct_jlo,
               Array4<int> const& bct_jhi,
               Array4<int> const& bct_klo,
               Array4<int> const& bct_khi,
               Box& domain,
               BcList& bc_list,
               const int* ng);

#endif
