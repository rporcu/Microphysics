#ifndef _MFIX_REACTIONS_RATES_K_H_
#define _MFIX_REACTIONS_RATES_K_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Gpu.H>

#include <cmath>

#include <mfix_pc.H>
#include <mfix_fluid_parms.H>
#include <mfix_solids_parms.H>
#include <mfix_reactions_parms.H>


class HeterogeneousRatesUser
{
  public:
    AMREX_GPU_HOST_DEVICE
    HeterogeneousRatesUser () = default;

    /**
     *  \brief Returns: the calculated heterogeneous chemical reactions rates
     *
     *  Inputs: // TODO
     */
    AMREX_GPU_HOST_DEVICE
    void operator() (amrex::Real* R_q,
                     const ReactionsParms& reactions_parms,
                     const SolidsParms& solids_parms,
                     const amrex::Real* X_sn,
                     const amrex::Real ro_s,
                     const amrex::Real ep_s,
                     const amrex::Real T_s,
                     const amrex::RealVect& vel_s,
                     const FluidParms& fluid_parms,
                     const amrex::Real* X_gk,
                     const amrex::Real ro_g,
                     const amrex::Real ep_g,
                     const amrex::Real T_g,
                     const amrex::RealVect& vel_g,
                     const amrex::Real DP,
                     const amrex::Real p_g,
                     const amrex::RunOn& run_on) const;
};

class HomogeneousRatesUser
{
  public:
    AMREX_GPU_HOST_DEVICE
    HomogeneousRatesUser () = default;

    /**
     *  \brief Returns: the calculated homogeneous solids-phase chemical
     *  reaction rates
     *
     *  Inputs: // TODO
     */
    AMREX_GPU_HOST_DEVICE
    void operator() (amrex::Real*, // R_q,
                     const ReactionsParms&, // reactions_parms,
                     const SolidsParms&, // solids_parms,
                     const amrex::Real*, // X_sn,
                     const amrex::Real, // ro_s,
                     const amrex::Real /*ep_s*/) const;

    /**
     *  \brief Returns: the calculated homogeneous fluid-phase chemical
     *  reaction rates
     *
     *  Inputs: // TODO
     */
    AMREX_GPU_HOST_DEVICE
    void operator() (amrex::Real*, // R_q,
                     const ReactionsParms&, // reactions_parms,
                     const FluidParms&, // fluid_parms,
                     const amrex::Real*, // X_gk,
                     const amrex::Real, // ro_g,
                     const amrex::Real /*ep_g*/) const;
};

#endif
