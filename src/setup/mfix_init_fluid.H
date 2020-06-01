#ifndef __MFIX_INIT_FLUID_HPP_
#define __MFIX_INIT_FLUID_HPP_

#include <mfix.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

void init_fluid (const Box& sbx,
                 const Box& bx,
                 const Box& domain,
                 const MFIter& mfi,
                 LevelData& ld,
                 const Real dx,
                 const Real dy,
                 const Real dz,
                 const Real xlength,
                 const Real ylength,
                 const Real zlength,
                 bool test_tracer_conservation,
                 const int advect_enthalpy,
                 const int advect_fluid_species);

void init_periodic_vortices (const Box& bx,
                             const Box& domain,
                             FArrayBox& vel,
                             const Real dx,
                             const Real dy,
                             const Real dz);

void init_periodic_tracer (const Box& bx,
                           const Box& domain,
                           FArrayBox& vel_fab,
                           FArrayBox& trac_fab,
                           const Real dx,
                           const Real dy,
                           const Real dz);

void init_fluid_restart (const Box& bx,
                         const MFIter& mfi,
                         LevelData& ld,
                         const int advect_enthalpy,
                         const int advect_fluid_species);

void set_ic (const Box& sbx,
             const Box& domain,
             const Real dx,
             const Real dy,
             const Real dz,
             FArrayBox& vel_g);

#endif
