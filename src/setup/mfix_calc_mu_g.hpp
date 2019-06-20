#ifndef __MFIX_CALC_MU_G_HPP_
#define __MFIX_CALC_MU_G_HPP_

#include <mfix_F.H>
#include <mfix.H>
#include <mfix_des_F.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

void calc_mu_g(const Box& bx, FArrayBox& mu_g_fab);

#endif
