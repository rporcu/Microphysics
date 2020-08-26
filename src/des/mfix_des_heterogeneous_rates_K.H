#ifndef _MFIX_DES_HETEROGENEOUS_K_H_
#define _MFIX_DES_HETEROGENEOUS_K_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_GpuQualifiers.H>

#include <cmath>

#include <mfix_pc.H>
#include <mfix_reactions_parms.H>

using ChemicalReaction = REACTIONS::ChemicalReaction;

struct ComputeRRateUser
{
    ComputeRRateUser () {};

    /**
     *  \brief Returns: the calculated drag coefficient.
     *
     *  Inputs: // TODO
     */
    AMREX_GPU_EXTERNAL
    AMREX_GPU_HOST_DEVICE
    void operator() (amrex::Real* R_q,
                     const int nreactions,
                     const int* nreactants,
                     const int* nproducts,
                     const int** reactants_id,
                     const amrex::Real** reactants_coeffs,
                     const int** reactants_phases,
                     const int** products_id,
                     const amrex::Real** products_coeffs,
                     const int** products_phases,
                     const int* species_id_s,
                     const amrex::Real* X_sn,
                     const amrex::Real* MW_sn,
                     const int nspecies_s,
                     const amrex::Real ro_s,
                     const amrex::Real ep_s,
                     const int* species_id_g,
                     const amrex::Real* X_gk,
                     const amrex::Real* MW_gk,
                     const int nspecies_g,
                     const amrex::Real ro_g,
                     const amrex::Real ep_g) const;
};

#endif
