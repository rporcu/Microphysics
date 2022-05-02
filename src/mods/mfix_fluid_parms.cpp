#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_algorithm.H>


using namespace amrex;

FluidPhase::FluidPhase()
  : ViscosityModel(VISCOSITYMODEL::Invalid)
  , SpecificHeatModel(SPECIFICHEATMODEL::Invalid)
  , ThermalConductivityModel(THERMALCONDUCTIVITYMODEL::Invalid)
  , solve(0)
  , solve_density(0)
  , solve_tracer(0)
  , trac_0(0)
  , mu_g0(0)
  , mw_avg(0)
  , solve_enthalpy(0)
  , T_g0(273.15)
  , T_ref(0)
  , k_g0(0)
  , thermodynamic_pressure(-1.0)
  , thermodynamic_pressure_defined(0)
  , solve_species(0)
  , species_names(0)
  , species_IDs(0)
  , d_species_IDs(0)
  , nspecies(0)
  , MW_gk0(0)
  , d_MW_gk0(0)
  , D_g0(0)
  , is_a_mixture(0)
  , H_fk0(0)
  , d_H_fk0(0)
  , cp_gk0(0)
  , d_cp_gk0(0)
  , stoich_coeffs(0)
  , d_stoich_coeffs(0)
  , name(std::string())
  , parameters(nullptr)
  , is_initialized(0)
{}

FluidPhase::~FluidPhase()
{
  if (parameters != nullptr)
    delete parameters;
}

void
FluidPhase::Initialize (const Species& species,
                        const Reactions& reactions)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.is_initialized,
      "Species not initialized. Can't initialize fluid phase before species initialization");

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(reactions.is_initialized,
      "Reactions not initialized. Can't initialize fluid phase before reactions initialization");

  // Flag for initialzation
  is_initialized = 1;

  amrex::ParmParse pp("fluid");

  std::vector<std::string> fluid_name;
  pp.queryarr("solve", fluid_name);

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid_name.size() == 1,
     "Fluid solver not specified. fluid.sove = ? ");

  // Disable the fluid solver if the fluid is defined as "None"
  if (amrex::toLower(fluid_name[0]).compare("none") == 0) {
    solve = 0;
  } else {
    solve = 1;
    name = fluid_name[0];
  }

  if (solve) {

    pp.query("T_g0",  T_g0);
    pp.query("trac0",  trac_0);
    pp.query("mw_avg", mw_avg);
    thermodynamic_pressure_defined = pp.query("thermodynamic_pressure",
                                              thermodynamic_pressure);

    // Constraint type
    {
      ParmParse ppMFIX("mfix");

      std::string constraint_str = "IncompressibleFluid";
      ppMFIX.query("constraint_type", constraint_str);
      constraint_str = amrex::toLower(constraint_str);

      if (constraint_str.compare("incompressiblefluid") == 0) {
        constraint_type = ConstraintType::IncompressibleFluid;
      }
      else if (constraint_str.compare("idealgasopensystem") == 0) {
        constraint_type = ConstraintType::IdealGasOpenSystem;
      }
      else if (constraint_str.compare("idealgasclosedsystem") == 0) {
        constraint_type = ConstraintType::IdealGasClosedSystem;
      }
      else {
        amrex::Abort("Don't know this constraint type!");
      }
    }

    if (constraint_type == ConstraintType::IncompressibleFluid &&
        thermodynamic_pressure_defined) {

      amrex::Warning("When the incompressible fluid constraint is selected, "
          "the fluid thermodynamic pressure input will be ignored");
    }

    if (constraint_type == ConstraintType::IdealGasClosedSystem &&
        (!thermodynamic_pressure_defined)) {

      amrex::Abort("When the idealgas closedsystem constraint is selected, "
          "the fluid thermodynamic pressure input must be provided");
    }

    amrex::ParmParse ppFluid(name.c_str());

    // Get viscosity inputs ----------------------------------//
    std::string viscosity_model;
    ppFluid.get("viscosity", viscosity_model);

    if (amrex::toLower(viscosity_model).compare("constant") == 0) {
      ViscosityModel = VISCOSITYMODEL::Constant;
      ppFluid.get("viscosity.constant", mu_g0);

    } else if (amrex::toLower(viscosity_model).compare("sutherland") == 0) {
      ViscosityModel = VISCOSITYMODEL::Sutherland;
      amrex::Abort("Not yet implemented.");

    } else {
      amrex::Abort("Unknown fluid viscosity model!");
    }


    // Get fluid temperature inputs ----------------------------------//
    amrex::ParmParse ppMFIX("mfix");

    ppMFIX.query("advect_density", solve_density);

    ppMFIX.query("advect_enthalpy", solve_enthalpy);

    ppMFIX.query("advect_tracer", solve_tracer);

    // Query mfix solve fluid species
    ppMFIX.query("solve_species", solve_species);

    // Query fluid species
    ppFluid.queryarr("species", species_names);

    if (!solve_species) {
      species_names.clear();
      nspecies = 0;
    } else if (amrex::toLower(species_names[0]).compare("none") == 0) {
      species_names.clear();
      solve_species = 0;
      nspecies = 0;
    } else {
      solve_species = 1;
      nspecies = species_names.size();

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species_names.size() > 0, 
                                       "No input provided for fluid.names_names");

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nspecies <= species.nspecies,
          "Fluid species number is higher than total species number");
    }

    if (solve_enthalpy == 1) {
      solve_enthalpy = 1;

      // Query the reference temperature
      ppFluid.query("reference_temperature", T_ref);

      // Get thermal conductivity inputs -----------------------------//
      std::string thermal_conductivity_model;
      ppFluid.get("thermal_conductivity", thermal_conductivity_model);

      if (amrex::toLower(thermal_conductivity_model).compare("constant") == 0) {
        ThermalConductivityModel = THERMALCONDUCTIVITYMODEL::Constant;
        ppFluid.get("thermal_conductivity.constant", k_g0);

      } else {
        amrex::Abort("Unknown fluid thermal conductivity model!");
      }
    }

    if (solve_enthalpy && !solve_species) {

      H_fk0.resize(1);

      // Get specific heat inputs ------------------------------------//
      std::string specific_heat_model;
      ppFluid.get("specific_heat", specific_heat_model);

      if (amrex::toLower(specific_heat_model).compare("constant") == 0) {
        SpecificHeatModel = SPECIFICHEATMODEL::Constant;
        cp_gk0.resize(1);
        ppFluid.get("specific_heat.constant", cp_gk0[0]);

        // Query the enthalpy_of_formation
        ppFluid.query("enthalpy_of_formation", H_fk0[0]);

      } else if (amrex::toLower(specific_heat_model).compare("nasa7-poly") == 0) {

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(MW_gk0[0] > 0., "Wrong fluid molecular weight");
        const amrex::Real coeff = FluidPhase::R / MW_gk0[0];

        SpecificHeatModel = SPECIFICHEATMODEL::NASA7Polynomials;
        cp_gk0.resize(12);

        for (int i(0); i < 6; ++i) {
          amrex::Vector<amrex::Real> aa(0);
          std::string field = "specific_heat.NASA7.a"+std::to_string(i);
          ppFluid.getarr(field.c_str(), aa);

          if (aa.size() == 1) {

            cp_gk0[i] = aa[0] * coeff;
            cp_gk0[i+6] = 0.;

          } else if (aa.size() == 2) {

            cp_gk0[i] = aa[0] * coeff;
            cp_gk0[i+6] = aa[1] * coeff;

          } else {

            amrex::Abort("Input error");
          }
        }

      } else if (amrex::toLower(specific_heat_model).compare("nasa9-poly") == 0) {
        SpecificHeatModel = SPECIFICHEATMODEL::NASA9Polynomials;
        amrex::Abort("Not yet implemented.");

      } else {
        amrex::Abort("Unknown fluid specific heat model!");
      }
    }

    if (solve_species) {

      species_IDs.resize(nspecies);
      MW_gk0.resize(nspecies);
      D_g0 = species.D_0;

      if (solve_enthalpy) {

        if (nspecies > 1) {
          std::string specific_heat_model;
          ppFluid.get("specific_heat", specific_heat_model);

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::toLower(specific_heat_model).compare("mixture") == 0,
              "When solving fluid enthalpy and species, fluid specific heat model must be a mixture");
        }

        SpecificHeatModel = species.SpecificHeatModel;

        if (SpecificHeatModel == Species::SPECIFICHEATMODEL::Constant) {
          cp_gk0.resize(nspecies);

        } else if (SpecificHeatModel == Species::SPECIFICHEATMODEL::NASA7Polynomials) {
          cp_gk0.resize(nspecies*12);
        }

        H_fk0.resize(nspecies);
      }

      for (int n(0); n < nspecies; n++) {
        auto it = std::find(species.names.begin(), species.names.end(), species_names[n]);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != species.names.end(),
                                         "Fluid species missing in input");

        const auto pos = std::distance(species.names.begin(), it);

        species_IDs[n] = species.IDs[pos];
        MW_gk0[n] = species.MW_k0[pos];

        if (solve_enthalpy) {

          if (SpecificHeatModel == Species::SPECIFICHEATMODEL::Constant) {
            cp_gk0[n] = species.cp_k0[pos];

          } else if (SpecificHeatModel == Species::SPECIFICHEATMODEL::NASA7Polynomials) {
            std::copy(&species.cp_k0[pos*12], &species.cp_k0[pos*12] + 12, &cp_gk0[n*12]);
          }

          H_fk0[n]  = species.H_fk0[pos];
        }
      }

    } else {

      if (ppFluid.contains("molecular_weight")) {
        MW_gk0.resize(1);
        ppFluid.get("molecular_weight", MW_gk0[0]);
      }

      if (solve_enthalpy) {
        H_fk0.resize(1);

        // Get specific heat model input ------------------------//
        std::string specific_heat_model;
        ppFluid.query("specific_heat", specific_heat_model);

        if (amrex::toLower(specific_heat_model).compare("constant") == 0) {
          SpecificHeatModel = SPECIFICHEATMODEL::Constant;
          cp_gk0.resize(1);
          ppFluid.get("specific_heat.constant", cp_gk0[0]);

          ppFluid.query("enthalpy_of_formation", H_fk0[0]);

        } else if (amrex::toLower(specific_heat_model).compare("nasa7-poly") == 0) {

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(MW_gk0.size() > 0, "No fluid molecular weight input");
          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(MW_gk0[0] > 0., "Wrong fluid molecular weight");

          const amrex::Real coeff = FluidPhase::R / MW_gk0[0];

          SpecificHeatModel = SPECIFICHEATMODEL::NASA7Polynomials;
          cp_gk0.resize(12);

          for (int i(0); i < 6; ++i) {
            amrex::Vector<amrex::Real> aa(0);
            std::string field = "specific_heat.NASA7.a"+std::to_string(i);
            ppFluid.getarr(field.c_str(), aa);

            if (aa.size() == 1) {

              cp_gk0[i] = aa[0] * coeff;
              cp_gk0[i+6] = 0.;

            } else if (aa.size() == 2) {

              cp_gk0[i] = aa[0] * coeff;
              cp_gk0[i+6] = aa[1] * coeff;

            } else {

              amrex::Abort("Input error");
            }
          }

        } else {
          amrex::Abort("Don't know this specific heat model!");
        }
      }
    }

    // Flag to determine if we want to solve the fluid as a mixture
    is_a_mixture = static_cast<int>(nspecies > 1);

    // Create the stoichiometric table for the fluid phase, that is associate
    // the total stoichiometric coefficient for each fluid species in each
    // reaction
    if (reactions.solve) {

      const int nreactions = reactions.nreactions;

      constexpr int Fluid = CHEMICALPHASE::Fluid;
      constexpr int Heterogeneous = REACTIONTYPE::Heterogeneous;

      // Allocate space for necessary data
      stoich_coeffs.resize(nspecies*nreactions, 0.);

      for (int n_g(0); n_g < nspecies; n_g++) {
        // Get the ID of the current species n_g
        const int species_id = species_IDs[n_g];

        // Loop over reactions to compute each contribution
        for (int q(0); q < nreactions; q++) {

          ChemicalReaction* chem_reaction = reactions.get(q);
          const auto& phases = chem_reaction->get_phases();

          const auto& reactants_IDs = chem_reaction->get_reactants_ids();
          const auto& reactants_phases = chem_reaction->get_reactants_phases();
          const auto& reactants_coeffs = chem_reaction->get_reactants_coeffs();

          const auto& products_IDs = chem_reaction->get_products_ids();
          const auto& products_phases = chem_reaction->get_products_phases();
          const auto& products_coeffs = chem_reaction->get_products_coeffs();

          // Do something only if reaction is heterogeneous and contains
          // a solid compound
          if (chem_reaction->get_type() == Heterogeneous &&
              std::find(phases.begin(), phases.end(), Fluid) != phases.end()) {

            // Add reactant contribution (if any)
            {
              for (size_t pos(0); pos < reactants_IDs.size(); ++pos) {
                if (species_id == reactants_IDs[pos] && reactants_phases[pos] == Fluid) {
                  stoich_coeffs[n_g*nreactions+q] += reactants_coeffs[pos];
                }
              }
            }

            // Add products contribution (if any)
            {
              for (size_t pos(0); pos < products_IDs.size(); ++pos) {
                if (species_id == products_IDs[pos] && products_phases[pos] == Fluid) {
                  stoich_coeffs[n_g*nreactions+q] += products_coeffs[pos];
                }
              }
            }
          }
        }
      }
    }
  }

  // Check on inputs in case of Ideal Gas EOS
  if (constraint_type == ConstraintType::IdealGasOpenSystem ||
      constraint_type == ConstraintType::IdealGasClosedSystem) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(MW_gk0.size() > 0, "Inputs error: fluid molecular_weight not provided");

    for (size_t i(0); i < MW_gk0.size(); ++i) {
      int abort(0);
      if (MW_gk0[i] < 1.e-15) {
        Print() << "Invalid molecular weight for species " << species_names[i] << "\n";
      }

      if (abort)
        amrex::Abort("Inputs error");
    }
  }

  if (solve_species) {
    d_species_IDs.resize(species_IDs.size());
    Gpu::copyAsync(Gpu::hostToDevice, species_IDs.begin(), species_IDs.end(), d_species_IDs.begin());
  }
  const int* p_h_species_IDs = solve_species ? species_IDs.data() : nullptr;
  const int* p_d_species_IDs = solve_species ? d_species_IDs.data() : nullptr;

  const int MW_gk0_provided = static_cast<int>(MW_gk0.size() > 0);
  if (MW_gk0_provided) {
    d_MW_gk0.resize(MW_gk0.size());
    Gpu::copyAsync(Gpu::hostToDevice, MW_gk0.begin(), MW_gk0.end(), d_MW_gk0.begin());
  }
  const Real* p_h_MW_gk0 = MW_gk0_provided ? MW_gk0.data() : nullptr;
  const Real* p_d_MW_gk0 = MW_gk0_provided ? d_MW_gk0.data() : nullptr;

  if (solve_enthalpy) {
    d_H_fk0.resize(H_fk0.size());
    Gpu::copyAsync(Gpu::hostToDevice, H_fk0.begin(), H_fk0.end(), d_H_fk0.begin());
  }
  const Real* p_h_H_fk0 = solve_enthalpy ? H_fk0.data() : nullptr;
  const Real* p_d_H_fk0 = solve_enthalpy ? d_H_fk0.data() : nullptr;

  if (solve_enthalpy) {
    d_cp_gk0.resize(cp_gk0.size());
    Gpu::copyAsync(Gpu::hostToDevice, cp_gk0.begin(), cp_gk0.end(), d_cp_gk0.begin());
  }
  const Real* p_h_cp_gk0 = solve_enthalpy ? cp_gk0.data() : nullptr;
  const Real* p_d_cp_gk0 = solve_enthalpy ? d_cp_gk0.data() : nullptr;

  if (reactions.solve) {
    d_stoich_coeffs.resize(stoich_coeffs.size());
    Gpu::copyAsync(Gpu::hostToDevice, stoich_coeffs.begin(), stoich_coeffs.end(), d_stoich_coeffs.begin());
  }
  const Real* p_h_stoich_coeffs = reactions.solve ? stoich_coeffs.data() : nullptr;
  const Real* p_d_stoich_coeffs = reactions.solve ? d_stoich_coeffs.data() : nullptr;

  int ncoefficients = 0;
  if (SpecificHeatModel == SPECIFICHEATMODEL::Constant)
    ncoefficients = 1;
  else if (SpecificHeatModel == SPECIFICHEATMODEL::NASA7Polynomials)
    ncoefficients = 6;

  parameters = new FluidParms(T_ref, mu_g0, k_g0, nspecies, p_h_species_IDs,
                              p_d_species_IDs, p_h_MW_gk0, p_d_MW_gk0, D_g0,
                              ncoefficients, p_h_cp_gk0, p_d_cp_gk0, p_h_H_fk0,
                              p_d_H_fk0, reactions.nreactions, p_h_stoich_coeffs,
                              p_d_stoich_coeffs, SpecificHeatModel);
}
