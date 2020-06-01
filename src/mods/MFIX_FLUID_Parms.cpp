#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <MFIX_FLUID_Parms.H>
#include <MFIX_SPECIES_Parms.H>


namespace FLUID
{

  int solve;

  // Specified constant gas temperature
  amrex::Real T_g0 = 273.15;

  // Specified constant gas specific heat
  amrex::Real cp_g0 = 1.0;

  // Gas phase thermal conductivity coefficient
  amrex::Real k_g0 = 0.0;

  // Specified constant gas density
  amrex::Real ro_g0;

  // Specified constant tracer value
  amrex::Real trac_0 = 0.0;

  // Specified constant gas viscosity
  amrex::Real mu_g0;

  // Average molecular weight of gas
  amrex::Real mw_avg;

  // Flag to solve species fluid equations
  int solve_species = 0;

  // Fluid phase species names
  std::vector<std::string> species_g;

  // Total number of fluid species
  int nspecies_g = 0;

  // Fluid phase species diffusion coefficients
  std::vector<amrex::Real> D_g0;

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
    if (fluid_name[0] == "None" ||
        fluid_name[0] == "none" ||
        fluid_name[0] == "NONE" ||
        fluid_name[0] == "0" ) {
      solve = 0;
    } else {
      solve = 1;
      name = fluid_name[0];
    }


    if( solve ) {

      amrex::ParmParse ppFluid(name.c_str());

      // Get density inputs ------------------------------------//
      std::string density_model;
      ppFluid.get("density", density_model );

      if (density_model == "constant")
      {
        //DENSITYMODEL DensityModel = ConstantDensity;
        ppFluid.get("density.constant",   ro_g0 );
      }
      else if(density_model == "idealgas")
      {
        //DENSITYMODEL DensityModel = IdealGas;
        amrex::Abort("Not yet implemented.");
      }
      else
      {
        amrex::Abort("Unknown fluid density model!");
      }


      // Get viscosity inputs ----------------------------------//
      std::string viscosity_model;
      ppFluid.get("viscosity", viscosity_model );

      if (viscosity_model == "constant")
      {
        //VISCOSITYMODEL ViscosityModel = ConstantViscosity;
        ppFluid.get("viscosity.constant", mu_g0 );
      }
      else if (viscosity_model == "sutherland")
      {
        //VISCOSITYMODEL ViscosityModel = Sutherland;
        amrex::Abort("Not yet implemented.");
        // TODO: get const, tg_ref, mu_ref
      }
      else
      {
        amrex::Abort("Unknown fluid viscosity model!");
      }


      if (ppFluid.contains("specific_heat")) {
        // Get specific heat inputs ------------------------------------//
        std::string specific_heat_model;
        ppFluid.query("specific_heat", specific_heat_model );

        if (specific_heat_model == "constant")
        {
          //SPECIFICHEATMODEL SpecificHeatModel = ConstantSpecificHeat;
          // If this value is not set in the inputs file, the default is 1.0
          ppFluid.get("specific_heat.constant",   cp_g0 );
        }
        else if (specific_heat_model == "nasa9-poly")
        {
          //SPECIFICHEATMODEL SpecificHeatModel = NASA9_Polynomial;
          amrex::Abort("Not yet implemented.");
          // TODO: get Tlow, Thigh, coefficients
        }
        else 
        {
          amrex::Abort("Unknown fluid specific heat model!");
        }
      }

      if (ppFluid.contains("thermal_conductivity")) {
        // Get specific heat inputs ------------------------------------//
        std::string thermal_conductivity_model;
        ppFluid.query("thermal_conductivity", thermal_conductivity_model );

        if (thermal_conductivity_model == "constant")
        {
          //THERMALCONDUCTIVITYMODEL ThermalConductivityModel = ConstantThermalConductivity;
          ppFluid.get("thermal_conductivity.constant", k_g0 );
        }
        else 
        {
          amrex::Abort("Unknown fluid thermal conductivity model!");
        }
      }


      // Fluid species inputs
      if (ppFluid.contains("species"))
      {
        ppFluid.getarr("species", species_g);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species_g.size() > 0, "No input provided for fluid.species");

        // Disable the species solver if the species are defined as "None" (case
        // insensitive) or 0
        if (amrex::toLower(species_g[0]).compare("none") == 0 or
            (species_g[0]).compare("0") == 0)
        {
          solve_species = 0;
          nspecies_g = 0;
        }
        else {
          solve_species = 1;
          nspecies_g = species_g.size();

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nspecies_g <= SPECIES::nspecies,
              "Fluid species_g number is higher than species_g number");

          D_g0.resize(nspecies_g);

          for (int n(0); n < nspecies_g; n++) {
            std::string one_specie_g = species_g[n];
            std::vector<std::string>::iterator it;

            it = std::find(SPECIES::species.begin(), SPECIES::species.end(), one_specie_g);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != SPECIES::species.end(),
                "Fluid specie " + one_specie_g + " is not included in provided species");

            D_g0[n] = SPECIES::D_0[std::distance(SPECIES::species.begin(), it)];
          }
        }
      }


      pp.query("T_g0",  T_g0 );
      pp.query("trac0",  trac_0 );
      pp.query("mw_avg", mw_avg );

    }
  }
}
