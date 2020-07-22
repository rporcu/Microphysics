#ifndef _EOS_MOD_F_H_
#define _EOS_MOD_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#include <cmath>


AMREX_GPU_HOST_DEVICE
amrex::Real idealgas (const amrex::Real& mw,
                      const amrex::Real& pg,
                      const amrex::Real& tg);

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Function: sutherland                                                !
//                                                                      !
//  Purpose: Compute a default value of gas viscosity where gas is      !
//  assumed to be air                                                   !
//                                                                      !
//  Literature/Document References:                                     !
//     Perry, R. H., and Chilton, C. H., Chemical Engineers' Handbook,  !
//        5th Edition, McGraw-Hill Inc., 1973, pp. 248, eqn. 3-133.     !
//     Arnold, J. H., Vapor viscosities and the Sutherland equation,    !
//        Journal of Chemical Physics, 1 (2), 1933, pp. 170-176.        !
//     Sutherland, W., The Viscosity of Gases and Molecular Force,      !
//        Phil. Mag. 5:507-531, 1893.                                   !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

AMREX_GPU_HOST_DEVICE
amrex::Real sutherland (const amrex::Real& tg);

#endif
