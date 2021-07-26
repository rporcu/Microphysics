#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <mfix_solids_parms.H>
#include <mfix_species_parms.H>
#include <AMReX_ParmParse.H>


using namespace amrex;

SolidsPhase::SolidsPhase()
  : NTYPES(0)
  , SpecificHeatModel(SPECIFICHEATMODEL::Invalid)
  , ThermalConductivityModel(THERMALCONDUCTIVITYMODEL::Invalid)
  , EnthalpyOfFormationModel(ENTHALPYOFFORMATIONMODEL::Invalid)
  , names(0)
  //, T_ref(298.15)
  , T_ref(0)
  , solve_species(0)
  , species(0)
  , species_id(0)
  , d_species_id(0)
  , nspecies(0)
  , MW_sn0(0)
  , d_MW_sn0(0)
  , is_a_mixture(0)
  , H_fn0(0)
  , d_H_fn0(0)
  , cp_sn0(0)
  , d_cp_sn0(0)
  , parameters(nullptr)
  , enthalpy_source(0)
  , update_mass(1)
  , update_momentum(1)
  , update_enthalpy(1)
{}


SolidsPhase::~SolidsPhase()
{
  if (parameters != nullptr)
    delete parameters;
}


void
SolidsPhase::Initialize ()
{
  int solve_energy = 0;

  {
    amrex::ParmParse pp_mfix("mfix");
    pp_mfix.query("advect_enthalpy", solve_energy);

    {
      amrex::ParmParse pp_mfix_solids("mfix.particles");
      pp_mfix_solids.query("enthalpy_source", enthalpy_source);
      pp_mfix_solids.query("update_mass", update_mass);
      pp_mfix_solids.query("update_momentum", update_momentum);
      pp_mfix_solids.query("update_enthalpy", update_enthalpy);
    }
  }

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
          cp_sn0.resize(NTYPES*5);
        } else if (amrex::toLower(specific_heat_model).compare("nasa9-poly") == 0) {
          SpecificHeatModel = SPECIFICHEATMODEL::NASA9Polynomials;
          amrex::Abort("Not yet implemented.");
        } else {
          amrex::Abort("Unknown fluid specific heat model!");
        }

        H_fn0.resize(NTYPES);
      }

      for(int solid(0); solid < NTYPES; ++solid) {

        amrex::ParmParse ppSolid(names[solid].c_str());

        if(solve_energy) {

          // Set up ParmParse to read in specific heat inputs
          std::string cp_str = names[solid]+".specific_heat";
          amrex::ParmParse pSOLIDS_CP(cp_str.c_str());

          if (SpecificHeatModel == SPECIFICHEATMODEL::Constant) {
            pSOLIDS_CP.get("constant", cp_sn0[solid]);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cp_sn0[solid] > 0, "Invalid DEM constant specific heat.");

          } else if (SpecificHeatModel == SPECIFICHEATMODEL::NASA7Polynomials) {
            pSOLIDS_CP.get("NASA7.a1", cp_sn0[solid*5 + 0]);
            pSOLIDS_CP.get("NASA7.a2", cp_sn0[solid*5 + 1]);
            pSOLIDS_CP.get("NASA7.a3", cp_sn0[solid*5 + 2]);
            pSOLIDS_CP.get("NASA7.a4", cp_sn0[solid*5 + 3]);
            pSOLIDS_CP.get("NASA7.a5", cp_sn0[solid*5 + 4]);
          } 

          // Query the enthalpy_of_formation
          ppSolid.query("enthalpy_of_formation", H_fn0[solid]);

        } // solve_energy
      }

      // Flag to determine if we want to solve the solid as a mixture
      is_a_mixture = false;

    } else if (NTYPES == 1) {

      amrex::ParmParse ppSolid(names[0].c_str());

      // Species
      solve_species = ppSolid.queryarr("species", species);

      if (!solve_species) {
        species.resize(1);
        species[0] = "None";
      }

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.size() > 0,
                                       "No input provided for solids.species");

      // Disable the species solver if the species are defined as "None"
      // (caseinsensitive) or 0
      if (amrex::toLower(species[0]).compare("none") == 0) {
        solve_species = 0;
        nspecies = 1; // anonymous species

      } else {
        solve_species = 1;
        nspecies = species.size();

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nspecies <= SPECIES::nspecies,
            "Solids species number is higher than total species number");
      }

      // Solids species inputs
      if (solve_species) {
        species_id.resize(nspecies);
        MW_sn0.resize(nspecies);

        if (solve_energy) {

          if (SPECIES::SpecificHeatModel == SPECIES::SPECIFICHEATMODEL::Constant) {
            cp_sn0.resize(nspecies);
          } else if (SPECIES::SpecificHeatModel == SPECIES::SPECIFICHEATMODEL::NASA7Polynomials) {
            cp_sn0.resize(nspecies*5);
          }

          H_fn0.resize(nspecies);
        }

        for (int n(0); n < nspecies; n++) {
          auto it = std::find(SPECIES::species.begin(), SPECIES::species.end(), species[n]);

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != SPECIES::species.end(),
                                           "Solid species missing in input");

          const auto pos = std::distance(SPECIES::species.begin(), it);

          species_id[n] = SPECIES::species_id[pos];
          MW_sn0[n] = SPECIES::MW_k0[pos];

          if (solve_energy) {
            if (SPECIES::SpecificHeatModel == SPECIES::SPECIFICHEATMODEL::Constant) {
              cp_sn0[n] = SPECIES::cp_k0[pos];
            } else if (SPECIES::SpecificHeatModel == SPECIES::SPECIFICHEATMODEL::NASA7Polynomials) {
              std::copy(&SPECIES::cp_k0[pos], &SPECIES::cp_k0[pos] + 5, &cp_sn0[n]);
            }

            H_fn0[n] = SPECIES::H_fk0[pos];
          }
        }
      } else {
        species_id.resize(1, 0);
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

          } else if (amrex::toLower(specific_heat_model).compare("nasa7-poly") == 0) {
            SpecificHeatModel = SPECIFICHEATMODEL::NASA7Polynomials;
            cp_sn0.resize(5);
            ppSolid.get("specific_heat.NASA7.a1", cp_sn0[0]);
            ppSolid.get("specific_heat.NASA7.a2", cp_sn0[1]);
            ppSolid.get("specific_heat.NASA7.a3", cp_sn0[2]);
            ppSolid.get("specific_heat.NASA7.a4", cp_sn0[3]);
            ppSolid.get("specific_heat.NASA7.a5", cp_sn0[4]);
          } else {
            amrex::Abort("Don't know this specific heat model!");
          }

          // Get enthalpy of formation model input ------------------------//
          std::string enthalpy_of_formation_model("constant");
          pp.query("enthalpy_of_formation", enthalpy_of_formation_model);

          if (amrex::toLower(enthalpy_of_formation_model).compare("constant") == 0) {

            EnthalpyOfFormationModel = ENTHALPYOFFORMATIONMODEL::Constant;
            ppSolid.query("enthalpy_of_formation.constant", H_fn0[0]);
          }
        }
      }

      // Flag to determine if we want to solve the solid as a mixture
      is_a_mixture = static_cast<int>(nspecies > 1);
    }

    d_species_id.resize(species_id.size());
    Gpu::copyAsync(Gpu::hostToDevice, species_id.begin(), species_id.end(), d_species_id.begin());
    const int* p_h_species_id = solve_species ? species_id.data() : nullptr;
    const int* p_d_species_id = solve_species ? d_species_id.data() : nullptr;

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

    int ncoefficients = 0;
    if (SpecificHeatModel == SPECIFICHEATMODEL::Constant)
      ncoefficients = 1;
    else if (SpecificHeatModel == SPECIFICHEATMODEL::NASA7Polynomials)
      ncoefficients = 5;

    parameters = new SolidsParms(T_ref, nspecies, p_h_species_id, p_d_species_id,
                                 p_h_MW_sn0, p_d_MW_sn0, ncoefficients,
                                 p_h_cp_sn0, p_d_cp_sn0, p_h_H_fn0, p_d_H_fn0);

  } else {
    parameters = new SolidsParms();
  }

}
