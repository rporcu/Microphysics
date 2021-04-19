#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <mfix_solids_parms.H>
#include <mfix_species_parms.H>
#include <AMReX_ParmParse.H>

namespace SOLIDS
{
  int NTYPES;

  // Fluid molecular weight model
  int MolecularWeightModel = MOLECULARWEIGHTMODEL::Invalid;

  // Fluid specific heat model
  int SpecificHeatModel = SPECIFICHEATMODEL::Invalid;

  // Fluid thermal conductivity model
  int ThermalConductivityModel = THERMALCONDUCTIVITYMODEL::Invalid;

  // Names of solids used in IC/BC setups
  amrex::Vector<std::string> names;

  // Specified constant solid phase molecular weight
  amrex::Real MW_s0(0);

  // Specified constant specific heat
  amrex::Vector<amrex::Real> cp_s0(0);

  // Specified constant solid enthalpy of formation
  amrex::Real H_f0(0);

  // Specified constant solid reference temperature
  amrex::Real T_ref(0);

  // Flag to solve species solids equations
  int solve_species(0);

  // Species unique identifying code
  amrex::Vector<int> species_id;

  // Particle species
  amrex::Vector<std::string> species;

  // Total number of dem species
  int nspecies;

  // Flag to set solids as a mixture of solids species
  int is_a_mixture(0);

  // Specified constant gas phase species molecular weight
  amrex::Vector<amrex::Real> MW_sn0(0);

  // Specified constant specific heat
  amrex::Vector<amrex::Real> cp_sn0(0);

  // Specified constant specific heat
  amrex::Vector<amrex::Real> H_fn0(0);

  amrex::Real enthalpy_source(0);


  void Initialize ()
  {
    int check_energy  = 1;

    {
      amrex::ParmParse pp_mfix("mfix");
      pp_mfix.query("advect_enthalpy", check_energy);

      {
        amrex::ParmParse pp_mfix_solids("mfix.particles");
        pp_mfix_solids.query("enthalpy_source", enthalpy_source);
      }
    }

    amrex::ParmParse pp("solids");

    pp.queryarr("types", names);

    NTYPES = names.size();

    for(int solid(0); solid < NTYPES; ++solid)
    {
      amrex::ParmParse ppSolid(names[solid].c_str());

      if(check_energy)
      {
        cp_s0.resize(NTYPES);

        if (ppSolid.contains("specific_heat"))
        {
          // Read in the specific heat model we are using
          std::string specific_heat_model;
          ppSolid.query("specific_heat", specific_heat_model );

          // Set up ParmParse to read in specific heat inputs
          std::string cp_str = names[solid]+".specific_heat";
          amrex::ParmParse pSOLIDS_CP(cp_str.c_str());

          if (amrex::toLower(specific_heat_model).compare("constant") == 0)
          {
            amrex::Real cp0_in(-1.);
            SpecificHeatModel = SPECIFICHEATMODEL::Constant;
            pSOLIDS_CP.get("constant", cp0_in);

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cp0_in > 0.0,
                "Invalid DEM constant specific heat.");

            cp_s0[solid] = cp0_in;
          }
          else if (amrex::toLower(specific_heat_model).compare("nasa9-poly") == 0)
          {
            SpecificHeatModel = SPECIFICHEATMODEL::NASA9Polynomials;
            amrex::Abort("Not yet implemented.");
          }
          else
          {
            amrex::Abort("Unknown fluid specific heat model!");
          }
        } // end specific heat

        // Get molecular weight inputs ---------------------------//
        std::string molecular_weight_model;
        ppSolid.query("molecular_weight", molecular_weight_model);

        if (amrex::toLower(molecular_weight_model).compare("constant") == 0)
        {
          MolecularWeightModel = MOLECULARWEIGHTMODEL::Constant;
          ppSolid.get("molecular_weight.constant", MW_s0);
        }
        else if(amrex::toLower(molecular_weight_model).compare("mixture") == 0)
        {
          MolecularWeightModel = MOLECULARWEIGHTMODEL::Mixture;
        }
        else
        {
          MolecularWeightModel = MOLECULARWEIGHTMODEL::Constant;

          if ( amrex::ParallelDescriptor::IOProcessor() )
            amrex::Warning("Solid molecular weight model not provided. "
              "Assuming constant model with MW_s = 0");
        }

        // Query the reference temperature
        pp.query("reference_temperature", T_ref);

        // Query the enthalpy_of_formation
        ppSolid.query("enthalpy_of_formation", H_f0);
      } // check_energy

      // Solids species inputs
      if (ppSolid.contains("species") ||
          MolecularWeightModel == MOLECULARWEIGHTMODEL::Mixture)
      {
        ppSolid.getarr("species", species);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.size() > 0,
            "No input provided for solids.species");

        // Disable the species solver if the species are defined as "None"
        // (caseinsensitive) or 0
        if (amrex::toLower(species[0]).compare("none") == 0 ||
            (species[0]).compare("0") == 0)
        {
          solve_species = 0;
          nspecies = 0;
        }
        else {
          solve_species = 1;
          nspecies = species.size();

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nspecies <= SPECIES::nspecies,
              "Solids species number is higher than total species number");

          species_id.resize(nspecies);
          MW_sn0.resize(nspecies);

          if (check_energy) {
            cp_sn0.resize(nspecies);
            H_fn0.resize(nspecies);
          }

          for (int n(0); n < nspecies; n++) {
            auto it = std::find(SPECIES::species.begin(), SPECIES::species.end(),
                species[n]);

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != SPECIES::species.end(),
                "Solid species missing in input");

            const auto pos = std::distance(SPECIES::species.begin(), it);

            species_id[n] = SPECIES::species_id[pos];
            MW_sn0[n] = SPECIES::MW_k0[pos];

            if (check_energy) {
              cp_sn0[n]  = SPECIES::cp_k0[pos];
              H_fn0[n]  = SPECIES::H_fk0[pos];
            }
          }
        }
      }

      // Flag to determine if we want to solve the solid as a mixture
      is_a_mixture = solve_species &&
        (MolecularWeightModel == MOLECULARWEIGHTMODEL::Mixture);

    }
  } // Initialize()

}
