#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <mfix_species_parms.H>
#include <mfix_fluid_parms.H>


namespace SPECIES
{
  int DiffusivityModel = DIFFUSIVITYMODEL::Invalid;

  int SpecificHeatModel = SPECIFICHEATMODEL::Invalid;

  // Flag to solve species equations
  int solve(0);

  // Total number of species
  int nspecies(0);

  // Species names
  amrex::Vector<std::string> species(0);

  // Species unique identifying code (at the moment = their index in the input
  // entries)
  amrex::Vector<int> species_id(0);

  // Specified species molecular weight
  amrex::Vector<amrex::Real> MW_k0(0);

  // Specified species diffusion coefficients
  amrex::Vector<amrex::Real> D_k0(0);

  // Flag to solve enthalpy species equations
  int solve_enthalpy(0);

  // Specified constant species specific heat
  amrex::Vector<amrex::Real> cp_k0(0);

  // Enthalpy of formation
  amrex::Vector<amrex::Real> H_fk0(0);


  void Initialize ()
  {
    amrex::ParmParse pp("species");

    if (pp.contains("solve")) {

      pp.getarr("solve", species);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.size() > 0,
                                       "No input provided for species.solve");

      // Disable the species solver if the species are defined as "None" (case
      // insensitive) or 0
      if (amrex::toLower(species[0]).compare("none") == 0) {
        solve = 0;

      } else {
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

          std::string specific_heat_model;
          pp.query("specific_heat", specific_heat_model);

          if (amrex::toLower(specific_heat_model).compare("constant") == 0) {
            SpecificHeatModel = SPECIFICHEATMODEL::Constant;
            cp_k0.resize(nspecies);
          } else if (amrex::toLower(specific_heat_model).compare("nasa7-poly") == 0) {
            SpecificHeatModel = SPECIFICHEATMODEL::NASA7Polynomials;
            cp_k0.resize(nspecies*12);
          } else {
            amrex::Abort("Don't know this specific heat model!");
          }

          H_fk0.resize(nspecies);
        }
      }

      if(solve) {
        // Get molecular weights input --------------------------------//
        for (int n(0); n < nspecies; n++) {
          std::string name = "species." + species[n];
          amrex::ParmParse ppSpecies(name.c_str());

          if (!ppSpecies.query("molecular_weight", MW_k0[n])) {

            if (amrex::ParallelDescriptor::IOProcessor()) {
              std::string message = "Input not provided. Assuming MW_" + species[n] + " = 0";
              amrex::Warning(message.c_str());
            }
          }
        }

        // Get diffusivity model input --------------------------------//
        std::string diffusivity_model;
        pp.get("diffusivity", diffusivity_model);

        if (amrex::toLower(diffusivity_model).compare("constant") == 0) {
          DiffusivityModel = DIFFUSIVITYMODEL::Constant;

          for (int n(0); n < nspecies; n++) {
            std::string name = "species." + species[n];
            amrex::ParmParse ppSpecies(name.c_str());
            int value_is_present = ppSpecies.query("diffusivity.constant", D_k0[n]);

            if (!value_is_present) {
              std::string message = "Assuming diffusivity D_" + species[n] + " = 0";
              amrex::Warning(message.c_str());
            }
          }

        } else {
          amrex::Abort("Unknown species mass diffusivity model!");
        }

        if (solve_enthalpy) {
          // Get specific heat model input ------------------------//
          if (SpecificHeatModel == SPECIFICHEATMODEL::Constant) {

            for (int n(0); n < nspecies; n++) {
              std::string name = "species." + species[n];
              amrex::ParmParse ppSpecies(name.c_str());
              ppSpecies.get("specific_heat.constant", cp_k0[n]);

              if(!ppSpecies.query("enthalpy_of_formation", H_fk0[n])) {

                if (amrex::ParallelDescriptor::IOProcessor()) {
                  std::string message = "Input not provided. Assuming Hf_" + species[n] + " = 0";
                  amrex::Warning(message.c_str());
                }
              }
            }

          } else if (SpecificHeatModel == SPECIFICHEATMODEL::NASA7Polynomials) {

            for (int n(0); n < nspecies; n++) {
              // Non-Normalization coefficient
              AMREX_ALWAYS_ASSERT_WITH_MESSAGE(MW_k0[n] > 0., "Wrong molecular weight inputs");
              const amrex::Real coeff = FluidPhase::R / MW_k0[n];

              std::string name = "species." + species[n];
              amrex::ParmParse ppSpecies(name.c_str());

              for (int i(0); i < 6; ++i) {
                amrex::Vector<amrex::Real> aa(0);
                std::string field = "specific_heat.NASA7.a"+std::to_string(i);
                ppSpecies.getarr(field.c_str(), aa);

                if (aa.size() == 1) {

                  cp_k0[n*12 + i] = aa[0] * coeff;
                  cp_k0[n*12 + i+6] = 0.;

                } else if (aa.size() == 2) {

                  cp_k0[n*12 + i] = aa[0] * coeff;
                  cp_k0[n*12 + i+6] = aa[1] * coeff;

                } else {

                  amrex::Abort("Input error");
                }
              }

            }

          } else {
            amrex::Abort("Unkown specific heat model");

            // Get enthalpy of formation model input ------------------------//
            for (int n(0); n < nspecies; n++) {
              std::string name = "species." + species[n];
              amrex::ParmParse ppSpecies(name.c_str());

              if(!ppSpecies.query("enthalpy_of_formation", H_fk0[n])) {

                if (amrex::ParallelDescriptor::IOProcessor()) {
                  std::string message = "Input not provided. Assuming Hf_" + species[n] + " = 0";
                  amrex::Warning(message.c_str());
                }
              }
            }
          }

        }
      }
    }
  }

}
