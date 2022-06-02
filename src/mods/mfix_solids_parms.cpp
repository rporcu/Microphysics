#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>

#include <mfix_solids_parms.H>
#include <mfix_species_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_algorithm.H>


using namespace amrex;

SolidsPhase::SolidsPhase()
  : NTYPES(0)
  , SpecificHeatModel(SPECIFICHEATMODEL::Invalid)
  , ThermalConductivityModel(THERMALCONDUCTIVITYMODEL::Invalid)
  , names(0)
  , solve_mass(1)
  , solve_momentum(1)
  , T_ref(0)
  , enthalpy_source(0)
  , solve_enthalpy(1)
  , solve_species(0)
  , species_names(0)
  , species_IDs(0)
  , d_species_IDs(0)
  , nspecies(0)
  , MW_sn0(0)
  , d_MW_sn0(0)
  , is_a_mixture(0)
  , H_fn0(0)
  , d_H_fn0(0)
  , cp_sn0(0)
  , d_cp_sn0(0)
  , stoich_coeffs(0)
  , d_stoich_coeffs(0)
  , kp_sn0(0)
  , d_kp_sn0(0)
  , flpc(0.4)
  , rough(2.E-8)
  , do_pfp_cond(0)
  , parameters(nullptr)
  , is_initialized(0)
{}


SolidsPhase::~SolidsPhase()
{
  if (parameters != nullptr)
    delete parameters;
}


void
SolidsPhase::Initialize (const Species& species,
                         const Reactions& reactions)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.is_initialized,
      "Species not initialized. Can't initialize solids phase before species initialization");

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(reactions.is_initialized,
      "Reactions not initialized. Can't initialize solids phase before reactions initialization");

  // Flag for initialization
  is_initialized = 1;

  int solve_energy = 0;

  {
    amrex::ParmParse pp_mfix("mfix");
    pp_mfix.query("advect_enthalpy", solve_energy);

    {
      amrex::ParmParse pp_mfix_solids("mfix.particles");
      pp_mfix_solids.query("enthalpy_source", enthalpy_source);
      pp_mfix_solids.query("update_mass", solve_mass);
      pp_mfix_solids.query("update_momentum", solve_momentum);
      pp_mfix_solids.query("update_enthalpy", solve_enthalpy);
    }
  }

  // Flag for pfp conduction
  int kp_type_read(0), kp_read(0), flpc_read(0), rough_read(0);
  int do_conduction(0);

  amrex::ParmParse pp("solids");

  pp.queryarr("types", names);

  if (names.size() > 0 && amrex::toLower(names[0]).compare("none") != 0) {

    NTYPES = names.size();

    if(solve_energy) {
      // Query the reference temperature
      pp.query("reference_temperature", T_ref);
    }

    // Currently we do not support species when NTYPES > 1
    if (NTYPES > 1) {
      amrex::Abort("Not yet implemented");

      MW_sn0.resize(NTYPES);

      if(solve_energy) {

        // Read in the specific heat model we are using
        std::string specific_heat_model;
        pp.get("specific_heat", specific_heat_model);

        if (amrex::toLower(specific_heat_model).compare("constant") == 0) {
          SpecificHeatModel = SPECIFICHEATMODEL::Constant;
          cp_sn0.resize(NTYPES);
        } else if (amrex::toLower(specific_heat_model).compare("nasa7-poly") == 0) {
          SpecificHeatModel = SPECIFICHEATMODEL::NASA7Polynomials;
          cp_sn0.resize(NTYPES*12);
        } else if (amrex::toLower(specific_heat_model).compare("nasa9-poly") == 0) {
          SpecificHeatModel = SPECIFICHEATMODEL::NASA9Polynomials;
          amrex::Abort("Not yet implemented.");
        } else {
          amrex::Abort("Unknown fluid specific heat model!");
        }

        H_fn0.resize(NTYPES);

        // Read in the conductivity model we are using
        std::string conductivity_model;
        kp_type_read = pp.query("thermal_conductivity", conductivity_model);

        if (amrex::toLower(conductivity_model).compare("constant") == 0) {
          ThermalConductivityModel = THERMALCONDUCTIVITYMODEL::Constant;
          kp_sn0.resize(NTYPES);
        } /*else {
          amrex::Abort("Unknown particle conductivity model!");
          }*/

      } // solve energy

      for(int solid(0); solid < NTYPES; ++solid) {

        amrex::ParmParse ppSolid(names[solid].c_str());

        ppSolid.query("molecular_weight", MW_sn0[solid]);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(MW_sn0[solid] > 0, "Invalid DEM molecular weight.");

        if(solve_energy) {

          // Set up ParmParse to read in specific heat inputs
          std::string cp_str = names[solid]+".specific_heat";
          amrex::ParmParse pSOLIDS_CP(cp_str.c_str());

          if (SpecificHeatModel == SPECIFICHEATMODEL::Constant) {
            pSOLIDS_CP.get("constant", cp_sn0[solid]);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cp_sn0[solid] > 0, "Invalid DEM constant specific heat.");

            // Query the enthalpy_of_formation
            ppSolid.query("enthalpy_of_formation", H_fn0[solid]);

          } else if (SpecificHeatModel == SPECIFICHEATMODEL::NASA7Polynomials) {

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cp_sn0[solid] > 0., "Wrong DEM molecular weight");
            const amrex::Real coeff = FluidParms::R / cp_sn0[solid];

            for (int i(0); i < 6; ++i) {
              amrex::Vector<amrex::Real> aa(0);
              std::string field = "NASA7.a"+std::to_string(i);
              pSOLIDS_CP.getarr(field.c_str(), aa);

              if (aa.size() == 1) {

                cp_sn0[solid*12 + i] = aa[0] * coeff;
                cp_sn0[solid*12 + i+6] = 0.;

              } else if (aa.size() == 2) {

                cp_sn0[solid*12 + i] = aa[0] * coeff;
                cp_sn0[solid*12 + i+6] = aa[1] * coeff;

              } else {

                amrex::Abort("Input error");
              }
            }
          }

          // Set up ParmParse to read in conductivity for each solids phase
          std::string kp_str = names[solid]+".thermal_conductivity";
          amrex::ParmParse pSOLIDS_KP(kp_str.c_str());

          if (ThermalConductivityModel == THERMALCONDUCTIVITYMODEL::Constant) {
            kp_read = pSOLIDS_KP.query("constant", kp_sn0[solid]);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kp_sn0[solid] > 0, "Invalid DEM thermal conductivity.");
          } /*else {
            amrex::Abort("Input error");
            }*/

        } // solve_energy
      } // for solid

      // Flag to determine if we want to solve the solid as a mixture
      is_a_mixture = false;

    } else if (NTYPES == 1) {

      amrex::ParmParse ppSolid(names[0].c_str());

      // Species
      solve_species = ppSolid.queryarr("species", species_names);

      if (!solve_species) {
        species_names.clear();
        nspecies = 0;
      } else {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species_names.size() > 0,
                                         "No input provided for solids.species");
      }

      // Solids species inputs
      if (solve_species) {

        // Disable the species solver if the species are defined as "None"
        // (caseinsensitive) or 0
        if (amrex::toLower(species_names[0]).compare("none") == 0) {
          solve_species = 0;
          nspecies = 0;

        } else {
          solve_species = 1;
          nspecies = species_names.size();

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nspecies <= species.nspecies,
              "Solids species number is higher than total species number");

          species_IDs.resize(nspecies);
          MW_sn0.resize(nspecies);

          if (solve_energy) {

            if (species.SpecificHeatModel == Species::SPECIFICHEATMODEL::Constant) {
              SpecificHeatModel = SPECIFICHEATMODEL::Constant;
              cp_sn0.resize(nspecies);
            } else if (species.SpecificHeatModel == Species::SPECIFICHEATMODEL::NASA7Polynomials) {
              SpecificHeatModel = SPECIFICHEATMODEL::NASA7Polynomials;
              cp_sn0.resize(nspecies*12);
            }

            H_fn0.resize(nspecies);
          }

          for (int n(0); n < nspecies; n++) {
            auto it = std::find(species.names.begin(), species.names.end(), species_names[n]);

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != species.names.end(),
                                             "Solid species missing in input");

            const auto pos = std::distance(species.names.begin(), it);

            species_IDs[n] = species.IDs[pos];
            MW_sn0[n] = species.MW_k0[pos];

            if (solve_energy) {
              if (SpecificHeatModel == SPECIFICHEATMODEL::Constant) {
                cp_sn0[n] = species.cp_k0[pos];
              } else if (SpecificHeatModel == SPECIFICHEATMODEL::NASA7Polynomials) {
                std::copy(&species.cp_k0[pos*12], &species.cp_k0[pos*12] + 12, &cp_sn0[n*12]);
              }

              H_fn0[n] = species.H_fk0[pos];
            }
          }
        }
      } else {
        species_IDs.resize(1, 0);
        MW_sn0.resize(1);

        ppSolid.query("molecular_weight", MW_sn0[0]);

        if (solve_energy) {
          H_fn0.resize(1);

          // Get specific heat model input ------------------------//
          std::string specific_heat_model;
          ppSolid.query("specific_heat", specific_heat_model);

          if (amrex::toLower(specific_heat_model).compare("constant") == 0) {
            SpecificHeatModel = SPECIFICHEATMODEL::Constant;
            cp_sn0.resize(1);
            ppSolid.get("specific_heat.constant", cp_sn0[0]);

            // Get enthalpy of formation model input ------------------------//
            ppSolid.query("enthalpy_of_formation", H_fn0[0]);

          } else if (amrex::toLower(specific_heat_model).compare("nasa7-poly") == 0) {
            SpecificHeatModel = SPECIFICHEATMODEL::NASA7Polynomials;
            cp_sn0.resize(12);

            for (int i(0); i < 6; ++i) {

              AMREX_ALWAYS_ASSERT_WITH_MESSAGE(MW_sn0[0] > 0., "Wrong solids molecular weight");
              const amrex::Real coeff = FluidParms::R / MW_sn0[0];

              amrex::Vector<amrex::Real> aa(0);
              std::string field = "specific_heat.NASA7.a"+std::to_string(i);
              ppSolid.getarr(field.c_str(), aa);

              if (aa.size() == 1) {

                cp_sn0[i] = aa[0] * coeff;
                cp_sn0[i+6] = 0.;

              } else if (aa.size() == 2) {

                cp_sn0[i] = aa[0] * coeff;
                cp_sn0[i+6] = aa[1] * coeff;

              } else {

                amrex::Abort("Input error");
              }
            }

          } else {
            amrex::Abort("Don't know this specific heat model!");
          }
        }
      }

      // Flag to determine if we want to solve the solid as a mixture
      is_a_mixture = static_cast<int>(nspecies > 1);

      // Create the stoichiometric table for the fluid phase, that is associate
      // the total stoichiometric coefficient for each fluid species in each
      // reaction
      if (reactions.solve) {

        const int nreactions = reactions.nreactions;

        constexpr int Solid = CHEMICALPHASE::Solid;
        constexpr int Heterogeneous = REACTIONTYPE::Heterogeneous;

        // Allocate space for necessary data
        stoich_coeffs.resize(nspecies*nreactions, 0.);

        for (int n_s(0); n_s < nspecies; n_s++) {
          // Get the ID of the current species n_s
          const int species_id = species_IDs[n_s];

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
                std::find(phases.begin(), phases.end(), Solid) != phases.end()) {

              // Add reactant contribution (if any)
              {
                for (int pos(0); pos < reactants_IDs.size(); ++pos) {
                  if (species_id == reactants_IDs[pos] && reactants_phases[pos] == Solid) {
                    stoich_coeffs[n_s*nreactions+q] += reactants_coeffs[pos];
                  }
                }
              }

              // Add products contribution (if any)
              {
                for (int pos(0); pos < products_IDs.size(); ++pos) {
                  if (species_id == products_IDs[pos] && products_phases[pos] == Solid) {
                    stoich_coeffs[n_s*nreactions+q] += products_coeffs[pos];
                  }
                }
              }
            }
          }
        }
      }

      // Do not support multiple species for solids conductivities!
      // Always reads in the conductivity model from first solids phase
      std::string conductivity_model;
      kp_type_read = ppSolid.query("thermal_conductivity", conductivity_model);

      if (amrex::toLower(conductivity_model).compare("constant") == 0) {
          ThermalConductivityModel = THERMALCONDUCTIVITYMODEL::Constant;
          kp_sn0.resize(1);
          kp_read = ppSolid.query("thermal_conductivity.constant", kp_sn0[0]);
          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kp_sn0[0] > 0, "Invalid DEM thermal conductivity.");
      } /*else {
          amrex::Abort("Unknown particle conductivity model!");
          }*/

  } // ntypes == 1

    // FLPC and roughness are always dimension 1
    if (solve_energy) {
        amrex::ParmParse ppSolid(names[0].c_str());
        flpc_read  = ppSolid.query("flpc", flpc);
        rough_read = ppSolid.query("min_conduction_dist", rough);
    }

    if(kp_type_read && kp_read && flpc_read && rough_read) do_conduction = 1;

    d_species_IDs.resize(species_IDs.size());
    Gpu::copyAsync(Gpu::hostToDevice, species_IDs.begin(), species_IDs.end(), d_species_IDs.begin());
    const int* p_h_species_IDs = solve_species ? species_IDs.data() : nullptr;
    const int* p_d_species_IDs = solve_species ? d_species_IDs.data() : nullptr;

    d_MW_sn0.resize(MW_sn0.size());
    Gpu::copyAsync(Gpu::hostToDevice, MW_sn0.begin(), MW_sn0.end(), d_MW_sn0.begin());
    const Real* p_h_MW_sn0 = solve_species ? MW_sn0.data() : nullptr;
    const Real* p_d_MW_sn0 = solve_species ? d_MW_sn0.data() : nullptr;

    if (solve_energy) {
      d_H_fn0.resize(H_fn0.size());
      Gpu::copyAsync(Gpu::hostToDevice, H_fn0.begin(), H_fn0.end(), d_H_fn0.begin());
    }
    const Real* p_h_H_fn0 = solve_energy ? H_fn0.data() : nullptr;
    const Real* p_d_H_fn0 = solve_energy ? d_H_fn0.data() : nullptr;

    if (solve_energy) {
      d_cp_sn0.resize(cp_sn0.size());
      Gpu::copyAsync(Gpu::hostToDevice, cp_sn0.begin(), cp_sn0.end(), d_cp_sn0.begin());
    }
    const Real* p_h_cp_sn0 = solve_energy ? cp_sn0.data() : nullptr;
    const Real* p_d_cp_sn0 = solve_energy ? d_cp_sn0.data() : nullptr;

    if (reactions.solve) {
      d_stoich_coeffs.resize(stoich_coeffs.size());
      Gpu::copyAsync(Gpu::hostToDevice, stoich_coeffs.begin(), stoich_coeffs.end(), d_stoich_coeffs.begin());
    }
    const Real* p_h_stoich_coeffs = reactions.solve ? stoich_coeffs.data() : nullptr;
    const Real* p_d_stoich_coeffs = reactions.solve ? d_stoich_coeffs.data() : nullptr;

    if (solve_energy) {
      d_kp_sn0.resize(kp_sn0.size());
      Gpu::copyAsync(Gpu::hostToDevice, kp_sn0.begin(), kp_sn0.end(), d_kp_sn0.begin());
    }
    const Real* p_h_kp_sn0 = solve_energy ? kp_sn0.data() : nullptr;
    const Real* p_d_kp_sn0 = solve_energy ? d_kp_sn0.data() : nullptr;

    int ncoefficients = 0;
    if (SpecificHeatModel == SPECIFICHEATMODEL::Constant)
      ncoefficients = 1;
    else if (SpecificHeatModel == SPECIFICHEATMODEL::NASA7Polynomials)
      ncoefficients = 6;

    parameters = new SolidsParms(T_ref, nspecies, p_h_species_IDs, p_d_species_IDs,
                                 p_h_MW_sn0, p_d_MW_sn0, ncoefficients,
                                 p_h_cp_sn0, p_d_cp_sn0, p_h_H_fn0, p_d_H_fn0,
                                 reactions.nreactions, p_h_stoich_coeffs, p_d_stoich_coeffs,
                                 p_h_kp_sn0, p_d_kp_sn0, flpc, rough, do_conduction,
                                 SpecificHeatModel, ThermalConductivityModel);

  } else {
    parameters = new SolidsParms();
  }

}
