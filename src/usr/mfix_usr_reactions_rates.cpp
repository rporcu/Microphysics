#include <AMReX_GpuQualifiers.H>
#include <AMReX_RealVect.H>

#include <cmath>

#include <mfix_des_rrates_K.H>
#include <mfix_species.H>
#include <mfix_fluid.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_reactions.H>
#include <mfix_algorithm.H>


using namespace amrex;


template <amrex::RunOn run_on>
AMREX_GPU_HOST_DEVICE
void
HeterogeneousRatesUser::operator() (Real* /*R_q*/,
                                    const MFIXReactionsParms& /*reactions_parms*/,
                                    const MFIXSolidsParms& /*solids_parms*/,
                                    const Real* /*X_sn*/,
                                    const Real /*ro_s*/,
                                    const Real /*ep_s*/,
                                    const Real /*T_s*/,
                                    const RealVect& /*vel_s*/,
                                    const MFIXFluidParms& /*fluid_parms*/,
                                    const Real* /*X_gk*/,
                                    const Real /*ro_g*/,
                                    const Real /*ep_g*/,
                                    const Real /*T_g*/,
                                    const RealVect& /*vel_g*/,
                                    const Real /*DP*/,
                                    const Real /*p_g*/) const
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


template <amrex::RunOn run_on>
AMREX_GPU_HOST_DEVICE
void
HomogeneousRatesUser::operator() (Real* /*R_q*/,
                                  const MFIXReactionsParms& /*reactions_parms*/,
                                  const MFIXSolidsParms& /*solids_parms*/,
                                  const Real* /*X_sn*/,
                                  const Real /*ro_s*/,
                                  const Real /*ep_s*/) const
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


template <amrex::RunOn run_on>
AMREX_GPU_HOST_DEVICE
void
HomogeneousRatesUser::operator() (Real* /*R_q*/,
                                  const MFIXReactionsParms& /*reactions_parms*/,
                                  const MFIXFluidParms& /*fluid_parms*/,
                                  const Real* /*X_gk*/,
                                  const Real /*ro_s*/,
                                  const Real /*ep_s*/) const
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


template
AMREX_GPU_HOST_DEVICE
void
HeterogeneousRatesUser::operator()<RunOn::Host>(Real* /*R_q*/,
                                                const MFIXReactionsParms& /*reactions_parms*/,
                                                const MFIXSolidsParms& /*solids_parms*/,
                                                const Real* /*X_sn*/,
                                                const Real /*ro_s*/,
                                                const Real /*ep_s*/,
                                                const Real /*T_s*/,
                                                const RealVect& /*vel_s*/,
                                                const MFIXFluidParms& /*fluid_parms*/,
                                                const Real* /*X_gk*/,
                                                const Real /*ro_g*/,
                                                const Real /*ep_g*/,
                                                const Real /*T_g*/,
                                                const RealVect& /*vel_g*/,
                                                const Real /*DP*/,
                                                const Real /*p_g*/) const;


template
AMREX_GPU_HOST_DEVICE
void
HeterogeneousRatesUser::operator()<RunOn::Gpu>(Real* /*R_q*/,
                                               const MFIXReactionsParms& /*reactions_parms*/,
                                               const MFIXSolidsParms& /*solids_parms*/,
                                               const Real* /*X_sn*/,
                                               const Real /*ro_s*/,
                                               const Real /*ep_s*/,
                                               const Real /*T_s*/,
                                               const RealVect& /*vel_s*/,
                                               const MFIXFluidParms& /*fluid_parms*/,
                                               const Real* /*X_gk*/,
                                               const Real /*ro_g*/,
                                               const Real /*ep_g*/,
                                               const Real /*T_g*/,
                                               const RealVect& /*vel_g*/,
                                               const Real /*DP*/,
                                               const Real /*p_g*/) const;


template
AMREX_GPU_HOST_DEVICE
void
HomogeneousRatesUser::operator()<RunOn::Host>(Real* /*R_q*/,
                                              const MFIXReactionsParms& /*reactions_parms*/,
                                              const MFIXSolidsParms& /*solids_parms*/,
                                              const Real* /*X_sn*/,
                                              const Real /*ro_s*/,
                                              const Real /*ep_s*/) const;


template
AMREX_GPU_HOST_DEVICE
void
HomogeneousRatesUser::operator()<RunOn::Gpu>(Real* /*R_q*/,
                                             const MFIXReactionsParms& /*reactions_parms*/,
                                             const MFIXSolidsParms& /*solids_parms*/,
                                             const Real* /*X_sn*/,
                                             const Real /*ro_s*/,
                                             const Real /*ep_s*/) const;


template
AMREX_GPU_HOST_DEVICE
void
HomogeneousRatesUser::operator()<RunOn::Host>(Real* /*R_q*/,
                                              const MFIXReactionsParms& /*reactions_parms*/,
                                              const MFIXFluidParms& /*fluid_parms*/,
                                              const Real* /*X_gk*/,
                                              const Real /*ro_s*/,
                                              const Real /*ep_s*/) const;


template
AMREX_GPU_HOST_DEVICE
void
HomogeneousRatesUser::operator()<RunOn::Gpu>(Real* /*R_q*/,
                                             const MFIXReactionsParms& /*reactions_parms*/,
                                             const MFIXFluidParms& /*fluid_parms*/,
                                             const Real* /*X_gk*/,
                                             const Real /*ro_s*/,
                                             const Real /*ep_s*/) const;
