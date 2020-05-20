#ifndef _MFIX_CALC_H_G_HPP_
#define _MFIX_CALC_H_G_HPP_

#include <mfix.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

void calc_k_g (const Box& bx, FArrayBox& temp_fab, FArrayBox& k_g_fab);

#endif
