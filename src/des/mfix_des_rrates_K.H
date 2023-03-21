#ifndef _MFIX_REACTIONS_RATES_K_H_
#define _MFIX_REACTIONS_RATES_K_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Gpu.H>

#include <cmath>

//#include <mfix_pc.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_reactions.H>


class HeterogeneousRatesUser
{
  public:
    AMREX_GPU_HOST_DEVICE
    HeterogeneousRatesUser () = default;

    /**
     *  \brief Returns: the calculated heterogeneous chemical reactions rates
     */
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE
#if defined(AMREX_USE_DPCPP)
    SYCL_EXTERNAL
#endif
    void operator() (amrex::Real* R_q,
                     const MFIXReactionsParms& reactions_parms,
                     const MFIXSolidsParms& solids_parms,
                     const amrex::Real* X_sn,
                     const amrex::Real ro_s,
                     const amrex::Real ep_s,
                     const amrex::Real T_s,
                     const amrex::RealVect& vel_s,
                     const MFIXFluidParms& fluid_parms,
                     const amrex::Real* X_gk,
                     const amrex::Real ro_g,
                     const amrex::Real ep_g,
                     const amrex::Real T_g,
                     const amrex::RealVect& vel_g,
                     const amrex::Real DP,
                     const amrex::Real p_g) const;
};

class HomogeneousRatesUser
{
  public:
    AMREX_GPU_HOST_DEVICE
    HomogeneousRatesUser () = default;

    /**
     *  \brief Returns: the calculated homogeneous solids-phase chemical
     *  reaction rates
     */
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE
#if defined(AMREX_USE_DPCPP)
    SYCL_EXTERNAL
#endif
    void operator() (amrex::Real*, // R_q,
                     const MFIXReactionsParms&, // reactions_parms,
                     const MFIXSolidsParms&, // solids_parms,
                     const amrex::Real*, // X_sn,
                     const amrex::Real, // ro_s,
                     const amrex::Real /*ep_s*/) const;

    /**
     *  \brief Returns: the calculated homogeneous fluid-phase chemical
     *  reaction rates
     */
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE
#if defined(AMREX_USE_DPCPP)
    SYCL_EXTERNAL
#endif
    void operator() (amrex::Real*, // R_q,
                     const MFIXReactionsParms&, // reactions_parms,
                     const MFIXFluidParms&, // fluid_parms,
                     const amrex::Real*, // X_gk,
                     const amrex::Real, // ro_g,
                     const amrex::Real /*ep_g*/) const;
};

#endif
