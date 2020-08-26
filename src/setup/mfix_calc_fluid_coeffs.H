#ifndef _MFIX_CALC_FLUID_COEFFS_HPP_
#define _MFIX_CALC_FLUID_COEFFS_HPP_

#include <mfix.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

// Compute fluid viscosity
void calc_mu_g (const Box& bx,
                FArrayBox& mu_g_fab);

// Compute fluid species diffusivity
void calc_D_gk (const Box& bx,
               FArrayBox& D_gk_fab);

// Compute fluid specific heat
void calc_cp_g (const Box& bx,
                FArrayBox& cp_g_fab,
                FArrayBox& T_g_fab);

// Compute fluid species specific heat
void calc_cp_gk (const Box& bx,
                 FArrayBox& cp_gk_fab,
                 FArrayBox& T_g_fab);

// Compute fluid enthalpy
void calc_h_g (const Box& bx,
               FArrayBox& h_g_fab,
               FArrayBox& cp_g_fab,
               FArrayBox& T_g_fab);

// Compute fluid species enthalpy
void calc_h_gk (const Box& bx,
                FArrayBox& h_gk_fab,
                FArrayBox& cp_gk_fab,
                FArrayBox& T_g_fab);

// Compute fluid diffusivity
void calc_k_g (const Box& bx,
               FArrayBox& k_g_fab,
               FArrayBox& T_g_fab);

#endif
