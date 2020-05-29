#ifndef _MFIX_NORMALIZE_SPECIES_G_H_
#define _MFIX_NORMALIZE_SPECIES_G_H_

#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

void normalize_species_g(const amrex::Vector< amrex::MultiFab* >& X_g);

#endif
