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
  , ro_g0(0)
  , MW_g0(0)
  , trac_0(0)
  , mu_g0(0)
  , mw_avg(0)
  , solve_enthalpy(0)
  , T_g0(273.15)
  , H_f0(0)
  , T_ref(0)
  , cp_g0(1)
  , k_g0(0)
  , solve_species(0)
  , species(0)
  , species_id(0)
  , d_species_id(0)
  , nspecies(0)
  , MW_gk0(0)
  , d_MW_gk0(0)
  , D_gk0(0)
  , d_D_gk0(0)
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
      if (ppFluid.contains("species")) {
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

              if (SPECIES::SpecificHeatModel == SPECIES::SPECIFICHEATMODEL::Invalid)
                cp_gk0[n] = cp_g0;
              else
                cp_gk0[n] = SPECIES::cp_k0[pos];

              H_fk0[n]  = SPECIES::H_fk0[pos];
            }
          }
        }
      }
    }

    // Flag to determine if we want to solve the fluid as a mixture
    is_a_mixture = static_cast<int>(nspecies > 1);

    pp.query("T_g0",  T_g0 );
    pp.query("trac0",  trac_0 );
    pp.query("mw_avg", mw_avg );
  }

  d_species_id.resize(species_id.size());
  Gpu::copyAsync(Gpu::hostToDevice, species_id.begin(), species_id.end(), d_species_id.begin());

  d_MW_gk0.resize(MW_gk0.size());
  Gpu::copyAsync(Gpu::hostToDevice, MW_gk0.begin(), MW_gk0.end(), d_MW_gk0.begin());
  Real* p_MW_gk0 = d_MW_gk0.data();

  d_D_gk0.resize(D_gk0.size());
  Gpu::copyAsync(Gpu::hostToDevice, D_gk0.begin(), D_gk0.end(), d_D_gk0.begin());
  Real* p_D_gk0 = d_D_gk0.data();

  d_H_fk0.resize(H_fk0.size());
  Gpu::copyAsync(Gpu::hostToDevice, H_fk0.begin(), H_fk0.end(), d_H_fk0.begin());

  d_cp_gk0.resize(cp_gk0.size());
  Gpu::copyAsync(Gpu::hostToDevice, cp_gk0.begin(), cp_gk0.end(), d_cp_gk0.begin());
  Real* p_cp_gk0 = d_cp_gk0.data();

  parameters = new FluidParms(MW_g0, p_MW_gk0, mu_g0, k_g0, p_D_gk0, cp_g0, p_cp_gk0);
}
