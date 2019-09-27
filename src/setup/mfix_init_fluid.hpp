#ifndef __MFIX_INIT_FLUID_HPP_
#define __MFIX_INIT_FLUID_HPP_

#include <mfix_F.H>
#include <mfix.H>
#include <mfix_des_F.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

void init_fluid(const Box& sbx,
                const Box& bx,
                const Box& domain,
                const FArrayBox& ep_g,
                FArrayBox& ro_g,
                FArrayBox& trac,
                const FArrayBox& p_g,
                FArrayBox& vel_g,
                FArrayBox& mu_g,
                const Real dx,
                const Real dy,
                const Real dz,
                const Real xlength,
                const Real ylength,
                const Real zlength);

void init_helix(const Box& bx,
                const Box& domain,
                FArrayBox& vel,
                const Real dx,
                const Real dy,
                const Real dz);

void init_periodic_vortices(const Box& bx,
                            const Box& domain,
                            FArrayBox& vel,
                            const Real dx,
                            const Real dy,
                            const Real dz);

void init_fluid_restart(const Box& bx,
                        FArrayBox& mu_g);

void set_ic(const Box& sbx,
            const Box& domain,
            const Real dx,
            const Real dy,
            const Real dz,
            FArrayBox& vel_g);

#endif
