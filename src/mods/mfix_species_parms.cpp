#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <mfix_species_parms.H>


namespace SPECIES
{
  int DiffusivityModel = DIFFUSIVITYMODEL::Invalid;

  int SpecificHeatModel = SPECIFICHEATMODEL::Invalid;

  // Flag to solve species equations
  int solve;

  // Total number of species
  int nspecies(0);

  // Species names
  std::vector<std::string> species(0);

  // Species unique identifying code (at the moment = their index in the input
  // entries)
  std::vector<int> species_id(0);

  // Specified species molecular weight
  std::vector<amrex::Real> MW_k0(0);

  // Specified species diffusion coefficients
  std::vector<amrex::Real> D_k0(0);

  // Flag to solve enthalpy species equations
  int solve_enthalpy(0);

  // Specified constant species specific heat
  std::vector<amrex::Real> cp_k0(0);


  void Initialize ()
  {
    amrex::ParmParse pp("species");

    if (pp.contains("solve"))
    {
      pp.getarr("solve", species);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.size() > 0,
          "No input provided for species.solve");

      // Disable the species solver if the species are defined as "None" (case
      // insensitive) or 0
      if (amrex::toLower(species[0]).compare("none") == 0 or
          (species[0]).compare("0") == 0)
      {
        solve = 0;
      }
      else
      {
        solve = 1;
        nspecies = species.size();

        species_id.resize(nspecies);
        for (int n(0); n < nspecies; n++) {
          species_id[n] = n;
        }

        MW_k0.resize(nspecies);
        D_k0.resize(nspecies);

        // Get species temperature inputs -----------------------------//
        amrex::ParmParse ppMFIX("mfix");
        int advect_enthalpy(0);
        ppMFIX.query("advect_enthalpy", advect_enthalpy);

        if (advect_enthalpy == 1) {
          solve_enthalpy = 1;
          cp_k0.resize(nspecies);
        }
      }

      if(solve)
      {
        // Get molecular weights input --------------------------------//
        for (int n(0); n < nspecies; n++) {
          std::string name = "species." + species[n];
          amrex::ParmParse ppSpecies(name.c_str());
          int exists = ppSpecies.query("molecular_weight", MW_k0[n]);

          if (not exists) {
            if ( amrex::ParallelDescriptor::IOProcessor() )
              amrex::Warning(species[n] + "_MW not provided. Assuming " +
                             species[n] + "_MW = 0");
          }
        }

        // Get diffusivity model input --------------------------------//
        std::string diffusivity_model;
        pp.get("diffusivity", diffusivity_model);

        if (amrex::toLower(diffusivity_model).compare("constant") == 0)
        {
          DiffusivityModel = DIFFUSIVITYMODEL::Constant;

          for (int n(0); n < nspecies; n++) {
            std::string name = "species." + species[n];
            amrex::ParmParse ppSpecies(name.c_str());
            ppSpecies.get("diffusivity.constant", D_k0[n]);
          }
        }
        else {
          amrex::Abort("Unknown species mass diffusivity model!");
        }

        if (solve_enthalpy)
        {
          // Get specific heat model input ------------------------//
          std::string specific_heat_model;
          pp.query("specific_heat", specific_heat_model);

          if (amrex::toLower(specific_heat_model).compare("constant") == 0)
          {
            SpecificHeatModel = SPECIFICHEATMODEL::Constant;

            for (int n(0); n < nspecies; n++) {
              std::string name = "species." + species[n];
              amrex::ParmParse ppSpecies(name.c_str());
              ppSpecies.get("specific_heat.constant", cp_k0[n]);
            }
          }
          else {
            SpecificHeatModel = SPECIFICHEATMODEL::Constant;

            if ( amrex::ParallelDescriptor::IOProcessor() )
              amrex::Warning("Species specific heat model not provided."
                             " Assuming constant model with cp_gk = 0");
          }

        }
      }
    }
  }

}
