#ifndef _MFIX_CALC_FLUID_COEFFS_HPP_
#define _MFIX_CALC_FLUID_COEFFS_HPP_

#include <mfix.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

void calc_mu_g (const Box& bx, FArrayBox& T_g_fab, FArrayBox& mu_g_fab);

void calc_D_g (const Box& bx, FArrayBox& temp_fab, FArrayBox& D_g_fab);

void calc_cp_g (const Box& bx, FArrayBox& temp_fab, FArrayBox& cp_g_fab);

void calc_k_g (const Box& bx, FArrayBox& temp_fab, FArrayBox& k_g_fab);

#endif
