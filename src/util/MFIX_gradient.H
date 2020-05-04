#ifndef _MFIX_GRADIENT_H_
#define _MFIX_GRADIENT_H_

#include <AMReX_REAL.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_Geometry.H>

namespace mfix_aux {

void computeGradient (amrex::MultiFab& grad,
                      const amrex::MultiFab& field,
                      const amrex::Geometry& geom);

void EB_computeGradient (amrex::MultiFab& grad,
                         const amrex::MultiFab& field,
                         const amrex::Geometry& geom,
                         bool already_on_centroids);

} // end namespace mfix_aux

#endif
