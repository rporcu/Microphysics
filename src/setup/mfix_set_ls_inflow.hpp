//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Subroutine: set_ls_inflow                                           !
//!  Purpose: This subroutine does the initial setting of all boundary   !
//!  conditions. The user specifications of the boundary conditions are  !
//!  checked for veracity in various check_data/ routines:               !
//!  (e.g., check_boundary_conditions).                                  !
//!                                                                      !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

#ifndef _MFIX_SET_LS_INFLOW_HPP_
#define _MFIX_SET_LS_INFLOW_HPP_

#include <mfix_F.H>
#include <mfix.H>
#include <mfix_des_F.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

void set_ls_inflow(FArrayBox& ls_phi_fab,
                   const BcList& bc_list,
                   const IArrayBox& bct_ilo_fab,
                   const IArrayBox& bct_ihi_fab,
                   const IArrayBox& bct_jlo_fab,
                   const IArrayBox& bct_jhi_fab,
                   const IArrayBox& bct_klo_fab,
                   const IArrayBox& bct_khi_fab,
                   const Box& domain,
                   const int* ng,
                   const int& nref,
                   const Real* dx);

#endif
