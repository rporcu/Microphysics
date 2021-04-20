#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>


namespace FLUID
{
  // Fluid density model
  int DensityModel = DENSITYMODEL::Invalid;

  // Fluid molecular weight model
  int MolecularWeightModel = MOLECULARWEIGHTMODEL::Invalid;

  // Fluid viscosity model
  int ViscosityModel = VISCOSITYMODEL::Invalid;

  // Fluid specific heat model
  int SpecificHeatModel = SPECIFICHEATMODEL::Invalid;

  // Fluid thermal conductivity model
  int ThermalConductivityModel = THERMALCONDUCTIVITYMODEL::Invalid;

  // Flag to solve fluid equations
  int solve;

  // Specified constant gas density
  amrex::Real ro_g0(0);

  // Specified constant gas molecular weight
  amrex::Real MW_g0(0);

  // Specified constant tracer value
  amrex::Real trac_0(0);

  // Specified constant gas viscosity
  amrex::Real mu_g0(0);

  // Average molecular weight of gas
  amrex::Real mw_avg(0);

  // Flag to solve enthalpy fluid equation
  int solve_enthalpy(0);

  // Specified constant gas temperature
  amrex::Real T_g0(273.15);

  // Specified constant gas specific heat
  amrex::Real cp_g0(1);

  // Specified constant gas enthalpy of formation
  amrex::Real H_f0(0);

  // Specified constant gas reference temperature
  // TODO TODO TODO
  amrex::Real T_ref(0);

  // Specified constant gas phase thermal conductivity coefficient
  amrex::Real k_g0(0);

  // Flag to solve species fluid equations
  int solve_species(0);

  // Fluid phase species names
  amrex::Vector<std::string> species;

  // Species unique identifying code
  amrex::Vector<int> species_id;

  // Total number of fluid species
  int nspecies(0);

  // Specified constant gas phase species molecular weight
  amrex::Vector<amrex::Real> MW_gk0(0);

  // Specified constant gas phase species diffusion coefficients
  amrex::Vector<amrex::Real> D_gk0(0);

  // Flag to understand if fluid is a mixture of fluid species
  int is_a_mixture(0);

  // Specified constant gas phase species specific heat
  amrex::Vector<amrex::Real> cp_gk0(0);

  // Enthalpy of formation
  amrex::Vector<amrex::Real> H_fk0(0);

  // Name to later reference when building inputs for IC/BC regions
  std::string name;


  void Initialize ()
  {
    amrex::ParmParse pp("fluid");

    std::vector<std::string> fluid_name;
    pp.queryarr("solve", fluid_name);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid_name.size() == 1,
       "Fluid solver not specified. fluid.sove = ? ");

    // Disable the fluid solver if the fluid is defined as "None"
    if (amrex::toLower(fluid_name[0]).compare("none") == 0 || fluid_name[0] == "0" ) {
      solve = 0;
    } else {
      solve = 1;
      name = fluid_name[0];
    }

    if (solve)
    {
      amrex::ParmParse ppFluid(name.c_str());

      // Get density inputs ------------------------------------//
      std::string density_model;
      ppFluid.get("density", density_model);

      if (amrex::toLower(density_model).compare("constant") == 0)
      {
        DensityModel = DENSITYMODEL::Constant;
        ppFluid.get("density.constant", ro_g0);
      }
      else if(amrex::toLower(density_model).compare("idealgas") == 0)
      {
        DensityModel = DENSITYMODEL::IdealGas;
        amrex::Abort("Not yet implemented.");
      }
      else
      {
        amrex::Abort("Unknown fluid density model!");
      }

      // Get molecular weight inputs ---------------------------//
      std::string molecular_weight_model;
      ppFluid.query("molecular_weight", molecular_weight_model);

      if (amrex::toLower(molecular_weight_model).compare("constant") == 0)
      {
        MolecularWeightModel = MOLECULARWEIGHTMODEL::Constant;
        ppFluid.get("molecular_weight.constant", MW_g0);
      }
      else if(amrex::toLower(molecular_weight_model).compare("mixture") == 0)
      {
        MolecularWeightModel = MOLECULARWEIGHTMODEL::Mixture;
      }
      else
      {
        MolecularWeightModel = MOLECULARWEIGHTMODEL::Constant;

        if ( amrex::ParallelDescriptor::IOProcessor() )
          amrex::Warning("Fluid molecular weight model not provided. "
            "Assuming constant model with MW_g = 0");
      }

      // Get viscosity inputs ----------------------------------//
      std::string viscosity_model;
      ppFluid.get("viscosity", viscosity_model );

      if (amrex::toLower(viscosity_model).compare("constant") == 0)
      {
        ViscosityModel = VISCOSITYMODEL::Constant;
        ppFluid.get("viscosity.constant", mu_g0);
      }
      else if (amrex::toLower(viscosity_model).compare("sutherland") == 0)
      {
        ViscosityModel = VISCOSITYMODEL::Sutherland;
        amrex::Abort("Not yet implemented.");
      }
      else
      {
        amrex::Abort("Unknown fluid viscosity model!");
      }


      // Get fluid temperature inputs ----------------------------------//
      amrex::ParmParse ppMFIX("mfix");

      int advect_enthalpy(0);
      ppMFIX.query("advect_enthalpy", advect_enthalpy);

      if (advect_enthalpy == 1) {
        solve_enthalpy = 1;

        if (MolecularWeightModel != MOLECULARWEIGHTMODEL::Mixture)
        {
          // Get specific heat inputs ------------------------------------//
          std::string specific_heat_model;
          ppFluid.get("specific_heat", specific_heat_model);

          if (amrex::toLower(specific_heat_model).compare("constant") == 0)
          {
            SpecificHeatModel = SPECIFICHEATMODEL::Constant;
            ppFluid.get("specific_heat.constant", cp_g0);
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

          // Query the enthalpy_of_formation
          ppFluid.query("enthalpy_of_formation", H_f0);
        }

        // Query the reference temperature
        ppFluid.query("reference_temperature", T_ref);

        // Get thermal conductivity inputs -----------------------------//
        std::string thermal_conductivity_model;
        ppFluid.get("thermal_conductivity", thermal_conductivity_model);

        if (amrex::toLower(thermal_conductivity_model).compare("constant") == 0)
        {
          ThermalConductivityModel = THERMALCONDUCTIVITYMODEL::Constant;
          ppFluid.get("thermal_conductivity.constant", k_g0);
        }
        else 
        {
          amrex::Abort("Unknown fluid thermal conductivity model!");
        }
      }

      int advect_fluid_species(0);
      ppMFIX.query("advect_fluid_species", advect_fluid_species);

      if (advect_fluid_species) {
        // Fluid species inputs
        if (ppFluid.contains("species") ||
            MolecularWeightModel == MOLECULARWEIGHTMODEL::Mixture)
        {
          ppFluid.getarr("species", species);

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.size() > 0, 
              "No input provided for fluid.species");

          // Disable the species solver if the species are defined as "None" (case
          // insensitive) or 0
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
                "Fluid species_g number is higher than species number");

            species_id.resize(nspecies);
            MW_gk0.resize(nspecies);
            D_gk0.resize(nspecies);

            if (solve_enthalpy) {
              cp_gk0.resize(nspecies);
              H_fk0.resize(nspecies);
            }

            for (int n(0); n < nspecies; n++) {
              auto it = std::find(SPECIES::species.begin(), SPECIES::species.end(),
                  species[n]);

              AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != SPECIES::species.end(),
                  "Fluid species missing in input");

              const auto pos = std::distance(SPECIES::species.begin(), it);

              species_id[n] = SPECIES::species_id[pos];
              MW_gk0[n] = SPECIES::MW_k0[pos];
              D_gk0[n] = SPECIES::D_k0[pos];

              if (solve_enthalpy) {
                cp_gk0[n] = SPECIES::cp_k0[pos];
                H_fk0[n]  = SPECIES::H_fk0[pos];
              }
            }
          }
        }
      }

      // Flag to determine if we want to solve the fluid as a mixture
      is_a_mixture = solve_species &&
        (FLUID::MolecularWeightModel == FLUID::MOLECULARWEIGHTMODEL::Mixture);

      pp.query("T_g0", T_g0);
      pp.query("trac0", trac_0);
      pp.query("mw_avg", mw_avg);

    }
  } // Initialize()
}
