#ifndef _MFIX_CALC_CP_G_HPP_
#define _MFIX_CALC_CP_G_HPP_

#include <mfix_F.H>
#include <mfix.H>
#include <mfix_des_F.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

void calc_cp_g (const Box& bx, FArrayBox& cp_g_fab);

#endif
