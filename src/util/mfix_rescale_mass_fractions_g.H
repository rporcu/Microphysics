#ifndef _MFIX_RESCALE_MASS_FRACTIONS_G_H_
#define _MFIX_RESCALE_MASS_FRACTIONS_G_H_

#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

void rescale_species(const amrex::Vector< amrex::MultiFab* >& X_g);

#endif
