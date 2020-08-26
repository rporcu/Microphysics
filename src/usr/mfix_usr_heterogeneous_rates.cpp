#include <AMReX_GpuQualifiers.H>

#include <cmath>

#include <mfix_des_heterogeneous_rates_K.H>
#include <mfix_species_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_algorithm.H>


AMREX_GPU_HOST_DEVICE
void
ComputeRRateUser::operator() (amrex::Real* R_q,
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
                              const amrex::Real ep_g) const
{
  // Loop over reactions
  for (int q(0); q < nreactions; q++)
  {
    amrex::Real k = 0.1;

    // Initially set it as constant
    R_q[q] = k;

  //  // Loop over reactants
  //  for (int n(0); n < nreactants[q]; n++)
  //  {
  //    const int current_species_id = reactants_id[q][n];

  //    // Current species concentration
  //    amrex::Real concentration(0.);

  //    if (reactants_phases[q][n] == CHEMICALPHASE::Solid)
  //    {
  //      const int pos = MFIXfind(species_id_s, nspecies_s, current_species_id);

  //      concentration = (ep_s*ro_s*X_sn[pos])/ (MW_sn[pos]);
  //    }
  //    else if (reactants_phases[q][n] == CHEMICALPHASE::Fluid)
  //    {
  //      const int pos = MFIXfind(species_id_s, nspecies_s, current_species_id);

  //      concentration = (ep_g*ro_g*X_gk[pos])/ (MW_gk[pos]);
  //    }

  //    R_q[q] *= concentration;
  //  }
  }

}
