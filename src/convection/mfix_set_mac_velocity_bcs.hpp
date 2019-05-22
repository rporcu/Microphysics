//               
//  This subroutine sets the BCs for velocity components only.
//  
//  Author: Michele Rosso
//  Date: December 20, 2017
// 

#ifndef _MFIX_SET_MAC_VELOCITY_BCS_HPP_
#define _MFIX_SET_MAC_VELOCITY_BCS_HPP_

#include <mfix_F.H>
#include <mfix.H>
#include <mfix_des_F.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

void
set_mac_velocity_bcs(Real* time,
                     BcList& bc_list,
                     const Box& bx,
                     MFIter* mfi,
                     MultiFab& u_g_fab,
                     MultiFab& v_g_fab,
                     MultiFab& w_g_fab,
                     IArrayBox& bct_ilo_fab,
                     IArrayBox& bct_ihi_fab,
                     IArrayBox& bct_jlo_fab,
                     IArrayBox& bct_jhi_fab,
                     IArrayBox& bct_klo_fab,
                     IArrayBox& bct_khi_fab,
                     Box& domain,
                     const int* nghost);

#endif
