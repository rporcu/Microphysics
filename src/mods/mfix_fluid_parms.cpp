#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>


using namespace amrex;

FluidPhase::FluidPhase()
  : DensityModel(DENSITYMODEL::Invalid)
  , ViscosityModel(VISCOSITYMODEL::Invalid)
  , SpecificHeatModel(SPECIFICHEATMODEL::Invalid)
  , ThermalConductivityModel(THERMALCONDUCTIVITYMODEL::Invalid)
  , solve(0)
  , solve_density(0)
  , ro_g0(0)
  , solve_tracer(0)
  , trac_0(0)
  , mu_g0(0)
  , mw_avg(0)
  , solve_enthalpy(0)
  , T_g0(273.15)
  //, T_ref(298.15)
  , T_ref(0)
  , k_g0(0)
  , solve_species(0)
  , species(0)
  , species_id(0)
  , d_species_id(0)
  , nspecies(0)
  , MW_gk0(0)
  , d_MW_gk0(0)
  , D_g0(0)
  , is_a_mixture(0)
  , H_fk0(0)
  , d_H_fk0(0)
  , cp_gk0(0)
  , d_cp_gk0(0)
  , name(std::string())
  , parameters(nullptr)
{}

FluidPhase::~FluidPhase()
{
  if (parameters != nullptr)
    delete parameters;
}

void
FluidPhase::Initialize ()
{
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

    amrex::ParmParse ppFluid(name.c_str());

    // Get density inputs ------------------------------------//
    std::string density_model;
    ppFluid.get("density", density_model);

    if (amrex::toLower(density_model).compare("constant") == 0) {
      DensityModel = DENSITYMODEL::Constant;
      ppFluid.get("density.constant", ro_g0);

    } else if(amrex::toLower(density_model).compare("idealgas") == 0) {
      DensityModel = DENSITYMODEL::IdealGas;
      amrex::Abort("Not yet implemented.");

    } else {
      amrex::Abort("Unknown fluid density model!");
    }

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

    int advect_enthalpy(0);
    ppMFIX.query("advect_enthalpy", advect_enthalpy);

    solve_species = ppFluid.queryarr("species", species);

    if (!solve_species) {
      species.clear();
      solve_species = 0;
      nspecies = 0;
    } else if (amrex::toLower(species[0]).compare("none") == 0) {
      species.clear();
      solve_species = 0;
      nspecies = 0;
    } else {
      solve_species = 1;
      nspecies = species.size();

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.size() > 0, 
                                       "No input provided for fluid.species");

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nspecies <= SPECIES::nspecies,
          "Fluid species number is higher than total species number");
    }

    if (advect_enthalpy == 1) {
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

    MW_gk0.resize(1);
    ppFluid.query("molecular_weight", MW_gk0[0]);

    if (advect_enthalpy && !solve_species) {

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

      species_id.resize(nspecies);
      MW_gk0.resize(nspecies);
      D_g0 = SPECIES::D_0;

      if (solve_enthalpy) {

        SpecificHeatModel = SPECIES::SpecificHeatModel;

        if (SpecificHeatModel == SPECIES::SPECIFICHEATMODEL::Constant) {
          cp_gk0.resize(nspecies);

        } else if (SpecificHeatModel == SPECIES::SPECIFICHEATMODEL::NASA7Polynomials) {
          cp_gk0.resize(nspecies*12);
        }

        H_fk0.resize(nspecies);
      }

      for (int n(0); n < nspecies; n++) {
        auto it = std::find(SPECIES::species.begin(), SPECIES::species.end(), species[n]);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != SPECIES::species.end(),
                                         "Fluid species missing in input");

        const auto pos = std::distance(SPECIES::species.begin(), it);

        species_id[n] = SPECIES::species_id[pos];
        MW_gk0[n] = SPECIES::MW_k0[pos];

        if (solve_enthalpy) {

          if (SpecificHeatModel == SPECIES::SPECIFICHEATMODEL::Constant) {
            cp_gk0[n] = SPECIES::cp_k0[pos];

          } else if (SpecificHeatModel == SPECIES::SPECIFICHEATMODEL::NASA7Polynomials) {
            std::copy(&SPECIES::cp_k0[pos*12], &SPECIES::cp_k0[pos*12] + 12, &cp_gk0[n*12]);
          }

          H_fk0[n]  = SPECIES::H_fk0[pos];
        }
      }
    } else {
      MW_gk0.resize(1);

      ppFluid.query("molecular_weight", MW_gk0[0]);

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
  }

  d_species_id.resize(species_id.size());
  Gpu::copyAsync(Gpu::hostToDevice, species_id.begin(), species_id.end(), d_species_id.begin());
  const int* p_h_species_id = species_id.data();
  const int* p_d_species_id = d_species_id.data();

  d_MW_gk0.resize(MW_gk0.size());
  Gpu::copyAsync(Gpu::hostToDevice, MW_gk0.begin(), MW_gk0.end(), d_MW_gk0.begin());
  const Real* p_h_MW_gk0 = MW_gk0.data();
  const Real* p_d_MW_gk0 = d_MW_gk0.data();

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

  int ncoefficients = 0;
  if (SpecificHeatModel == SPECIFICHEATMODEL::Constant)
    ncoefficients = 1;
  else if (SpecificHeatModel == SPECIFICHEATMODEL::NASA7Polynomials)
    ncoefficients = 6;

  parameters = new FluidParms(T_ref, mu_g0, k_g0, nspecies, p_h_species_id,
                              p_d_species_id, p_h_MW_gk0, p_d_MW_gk0, D_g0,
                              ncoefficients, p_h_cp_gk0, p_d_cp_gk0, p_h_H_fk0,
                              p_d_H_fk0, SpecificHeatModel);
}
