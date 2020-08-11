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

  // Particle species
  amrex::Vector<std::string> species;

  // Total number of dem species
  int nspecies;

  // Specified constant specific heat
  amrex::Vector<amrex::Real> cp_p0;

  void Initialize ()
  {

    int check_energy  = 0;
    // int check_species = 0;
    {
      amrex::ParmParse pp_mfix("mfix");
      pp_mfix.query("advect_enthalpy", check_energy );
    }

    amrex::ParmParse pp("solids");

    pp.queryarr("types", names);

    NTYPES = names.size();

    for(int solid(0); solid<NTYPES; ++solid) {

      amrex::ParmParse ppSOLIDS(names[solid].c_str());

      if(check_energy){

        if (ppSOLIDS.contains("specific_heat")) {

          cp_p0.resize(NTYPES);

          // Read in the specific heat model we are using
          std::string specific_heat_model;
          ppSOLIDS.query("specific_heat", specific_heat_model );

          // Set up ParmParse to read in specific heat inputs
          std::string cp_str = names[solid]+".specific_heat";
          amrex::ParmParse pSOLIDS_CP(cp_str.c_str());

          if (amrex::toLower(specific_heat_model).compare("constant") == 0)
          {
            amrex::Real cp0_in(-1.);
            SpecificHeatModel = SPECIFICHEATMODEL::Constant;
            pSOLIDS_CP.get("constant", cp0_in);

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cp0_in > 0.0, "Invalid DEM constant specific heat.");

            cp_p0[solid] = cp0_in;

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

      } // check_energy


///////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Read species inputs -----------------------------------------//
        // int species = 0;
        // ppSOLIDS.query("species", species);
        // if (species > 0)
        // {
        //       ppSOLIDS.getarr("species.names", species_dem);
        //       nspecies_dem = species_dem.size();
        // }
///////////////////////////////////////////////////////////////////////////////////////////////////////////
    }

  };

}
