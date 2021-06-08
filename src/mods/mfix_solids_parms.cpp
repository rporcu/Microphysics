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
  , names(0)
  , MW_s0(0)
  , cp_s0(0)
  , H_f0(0)
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
{}

SolidsPhase::~SolidsPhase()
{
  if (parameters != nullptr)
    delete parameters;
}

void
SolidsPhase::Initialize ()
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

      // Query the reference temperature
      pp.query("reference_temperature", T_ref);

      // Query the enthalpy_of_formation
      ppSolid.query("enthalpy_of_formation", H_f0);
    } // check_energy

    // Solids species inputs
    if (ppSolid.contains("species")) {
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
            if (SPECIES::SpecificHeatModel == SPECIES::SPECIFICHEATMODEL::Invalid)
              cp_sn0[n] = cp_s0[solid];
            else
              cp_sn0[n] = SPECIES::cp_k0[pos];

            H_fn0[n]  = SPECIES::H_fk0[pos];
          }
        }
      }
    }

    // Flag to determine if we want to solve the solid as a mixture
    is_a_mixture = static_cast<int>(nspecies > 1);
  }

  d_cp_s0.resize(cp_s0.size());
  Gpu::copyAsync(Gpu::hostToDevice, cp_s0.begin(), cp_s0.end(), d_cp_s0.begin());
  Real* p_h_cp_s0 = cp_s0.data();
  Real* p_d_cp_s0 = d_cp_s0.data();

  d_species_id.resize(species_id.size());
  Gpu::copyAsync(Gpu::hostToDevice, species_id.begin(), species_id.end(), d_species_id.begin());

  d_MW_sn0.resize(MW_sn0.size());
  Gpu::copyAsync(Gpu::hostToDevice, MW_sn0.begin(), MW_sn0.end(), d_MW_sn0.begin());

  d_H_fn0.resize(H_fn0.size());
  Gpu::copyAsync(Gpu::hostToDevice, H_fn0.begin(), H_fn0.end(), d_H_fn0.begin());

  d_cp_sn0.resize(cp_sn0.size());
  Gpu::copyAsync(Gpu::hostToDevice, cp_sn0.begin(), cp_sn0.end(), d_cp_sn0.begin());
  Real* p_h_cp_sn0 = cp_sn0.data();
  Real* p_d_cp_sn0 = d_cp_sn0.data();

  parameters = new SolidParms(p_h_cp_s0, p_d_cp_s0, p_h_cp_sn0, p_d_cp_sn0);
}
