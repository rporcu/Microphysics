#ifndef _MFIX_DES_USR_REACTIONS_RATES_K_H_
#define _MFIX_DES_USR_REACTIONS_RATES_K_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Gpu.H>

#include <cmath>

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
    template <amrex::RunOn run_on, class F1>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (amrex::Real* R_q,
                     F1 accessor,
                     const int idx,
                     const MFIXReactionsParms& reactions_parms,
                     const int /*np*/,
                     const int /*p_id*/,
                     const MFIXSolidsParms& /*solids_parms*/,
                     const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& /*ptile_data*/,
                     const int /*idx_X_sn*/,
                     const amrex::Real /*ro_s*/,
                     const amrex::Real /*ep_s*/,
                     const amrex::Real /*T_s*/,
                     const amrex::RealVect& /*vel_s*/,
                     const MFIXFluidParms& /*fluid_parms*/,
                     amrex::Real* /*X_gk*/,
                     const amrex::Real /*ro_g*/,
                     const amrex::Real /*ep_g*/,
                     const amrex::Real /*T_g*/,
                     const amrex::RealVect& /*vel_g*/,
                     const amrex::Real /*DP*/,
                     const amrex::Real /*p_g*/) const
    {
      // Loop over reactions
      for (int q(0); q < reactions_parms.m_nreactions; q++)
      {
//        int nreact(0);
//        for (int n(0); n < reactions_parms.nreactants[q]; n++)
//        {
//          const int current_species_id = reactions_parms.reactants_id[q][n];
//
//          {
//            const int pos = find(fluid_parms.species_id, fluid_parms.nspecies, current_species_id);
//            if (pos != -1)
//              if (X_gk[pos] > 0)
//                nreact += 1;
//          }
//
//          {
//            const int pos = find(solids_parms.species_id, solids_parms.nspecies, current_species_id);
//            if (pos != -1)
//              if (X_sn[pos] > 0)
//                nreact += 1;
//          }
//        }
//      
//        R_q[accessor(q,idx)] = (nreact == reactions_parms.nreactants[q]) ? 0 : reaction rate;

        R_q[accessor(q,idx)] = 0.;
      }
    }
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
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (amrex::Real* /*R_q*/,
                     const MFIXReactionsParms& /*reactions_parms*/,
                     const int /*np*/,
                     const int /*p_id*/,
                     const MFIXSolidsParms& /*solids_parms*/,
                     const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& /*ptile_data*/,
                     const int /*idx_X_sn*/,
                     const amrex::Real /*ro_s*/,
                     const amrex::Real /*ep_s*/) const
    {
//  // Loop over reactions
//  for (int q(0); q < reactions_parms.nreactions; q++)
//  {
//    int nreact(0);
//    for (int n(0); n < reactions_parms.nreactants[q]; n++)
//    {
//      const int current_species_id = reactions_parms.reactants_id[q][n];
//
//      {
//        const int pos = find(solids_parms.species_id, solids_parms.nspecies, current_species_id);
//        if (pos != -1)
//          if (X_sn[pos] > 0)
//            nreact += 1;
//      }
//    }
//  
//    R_q[q] = (nreact == reactions_parms.nreactants[q]) ? 0 : reaction rate;
//  }
//
    }

    /**
     *  \brief Returns: the calculated homogeneous fluid-phase chemical
     *  reaction rates
     */
    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (amrex::Real* /*R_q*/,
                     const MFIXReactionsParms& /*reactions_parms*/,
                     const int /*np*/,
                     const int /*p_id*/,
                     const MFIXFluidParms& /*fluid_parms*/,
                     amrex::Real* /*X_gk*/,
                     const amrex::Real /*ro_g*/,
                     const amrex::Real /*ep_g*/) const
    {
//  // Loop over reactions
//  for (int q(0); q < reactions_parms.nreactions; q++)
//  {
//    int nreact(0);
//    for (int n(0); n < reactions_parms.nreactants[q]; n++)
//    {
//      const int current_species_id = reactions_parms.reactants_id[q][n];
//
//      {
//        const int pos = find(fluid_parms.species_id, fluid_parms.nspecies, current_species_id);
//        if (pos != -1)
//          if (X_gk[pos] > 0)
//            nreact += 1;
//      }
//    }
//  
//    R_q[q] = (nreact == reactions_parms.nreactants[q]) ? 0 : reaction rate;
//  }
//
    }
};

#endif
