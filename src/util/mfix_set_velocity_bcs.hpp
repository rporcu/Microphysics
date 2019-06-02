#ifndef _SET_VELOCITY_BCS_HPP_
#define _SET_VELOCITY_BCS_HPP_

#include <mfix_F.H>
#include <mfix.H>
#include <mfix_des_F.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

//              
//  This subroutine sets the BCs for the velocity components only.
//

void
set_velocity_bcs(Real* time,
                 const BcList& bc_list,
                 FArrayBox& vel_fab,
                 const IArrayBox& bct_ilo_fab,
                 const IArrayBox& bct_ihi_fab,
                 const IArrayBox& bct_jlo_fab,
                 const IArrayBox& bct_jhi_fab,
                 const IArrayBox& bct_klo_fab,
                 const IArrayBox& bct_khi_fab,
                 const Box& domain,
                 const int* ng,
                 const int* extrap_dir_bcs);

void
set_vec_bcs(const BcList& bc_list,
            FArrayBox& vec_fab,
            const IArrayBox& bct_ilo_fab,
            const IArrayBox& bct_ihi_fab,
            const IArrayBox& bct_jlo_fab,
            const IArrayBox& bct_jhi_fab,
            const IArrayBox& bct_klo_fab,
            const IArrayBox& bct_khi_fab,
            const Box& domain,
            const int* ng);

#endif
