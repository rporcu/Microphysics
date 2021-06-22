#include <AMReX_GpuQualifiers.H>

#include <cmath>

#include <mfix_des_heterogeneous_rates_K.H>
#include <mfix_species_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_algorithm.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE
void
ComputeRRateUser::operator() (Real* R_q,
                              const int nreactions,
                              const int* nreactants,
                              const int* /*nproducts*/,
                              const int** reactants_id,
                              const Real** /*reactants_coeffs*/,
                              const int** /*reactants_phases*/,
                              const int** /*products_id*/,
                              const Real** /*products_coeffs*/,
                              const int** /*products_phases*/,
                              const int* species_id_s,
                              const Real* X_sn,
                              const Real* /*MW_sn*/,
                              const int nspecies_s,
                              const Real /*ro_s*/,
                              const Real /*ep_s*/,
                              const int* species_id_g,
                              const Real* X_gk,
                              const Real* /*MW_gk*/,
                              const int nspecies_g,
                              const Real /*ro_g*/,
                              const Real /*ep_g*/) const
{
  // Loop over reactions
  for (int q(0); q < nreactions; q++)
  {

    int nreact(0);
    for (int n(0); n < nreactants[q]; n++)
    {
      const int current_species_id = reactants_id[q][n];

      {
        const int pos = MFIXfind(species_id_g, nspecies_g, current_species_id);
        if (pos != -1)
          if (X_gk[pos] > 0)
            nreact += 1;
      }

      {
        const int pos = MFIXfind(species_id_s, nspecies_s, current_species_id);
        if (pos != -1)
          if (X_sn[pos] > 0)
            nreact += 1;
      }
    }
  R_q[q] = (nreact == nreactants[q]) ? 1.e-3 : 0.0;
  }

}
