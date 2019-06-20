//
//
//  This subroutine sets the BCs for all scalar variables involved in
//  the projection, EXCEPT the pressure and velocity components.
//  To set velocity BCs, use "set_velocity_bcs".
//
//  Author: Michele Rosso
//
//  Date: December 20, 2017
//
//

#ifndef _MFIX_SET_SCALAR_BCS_HPP_
#define _MFIX_SET_SCALAR_BCS_HPP_

#include <mfix_F.H>
#include <mfix.H>
#include <mfix_des_F.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

void set_scalar_bcs(const BcList& bc_list,
                    FArrayBox& ep_g_fab,
                    FArrayBox& ro_g_fab,
                    FArrayBox& mu_g_fab,
                    const IArrayBox& bct_ilo_fab,
                    const IArrayBox& bct_ihi_fab,
                    const IArrayBox& bct_jlo_fab,
                    const IArrayBox& bct_jhi_fab,
                    const IArrayBox& bct_klo_fab,
                    const IArrayBox& bct_khi_fab,
                    const Box& domain,
                    Real* m_bc_ep_g,
                    Real* m_bc_t_g,
                    const int* ng);

#endif
