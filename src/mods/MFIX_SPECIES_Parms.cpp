#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <MFIX_SPECIES_Parms.H>


namespace SPECIES
{
  // Total number of species
  int nspecies(0);

  // Species names
  std::vector<std::string> species(0);

  // Species diffusion coefficients
  std::vector<amrex::Real> D_0(0);


  void Initialize ()
  {
    amrex::ParmParse pp("species");

    if (pp.contains("solve"))
    {
      pp.getarr("solve", species);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.size() > 0, "No input provided for species.solve");

      // Disable the species solver if the species are defined as "None" (case
      // insensitive) or 0
      if (amrex::toLower(species[0]).compare("none") == 0 or
          (species[0]).compare("0") == 0) {
        nspecies = 0;
      }
      else {
        nspecies = species.size();
        D_0.resize(nspecies);
      }


      if( nspecies > 0 ) {

        // Get diffusivity model input --------------------------------//
        std::string diffusivity_model;
        pp.get("diffusivity", diffusivity_model);

        if (amrex::toLower(diffusivity_model).compare("constant") == 0) {

          //SPECIESDIFFUSIVITYMODEL SpeciesDiffusivityModel = ConstantSpeciesDiffusivity;
          for (int n(0); n < nspecies; n++) {
            std::string name = "species." + species[n];
            amrex::ParmParse ppSpecies(name.c_str());
            ppSpecies.get("diffusivity", D_0[n]);
          }
        }
        else {
          amrex::Abort("Unknown fluid species mass diffusivity model!");
        }
      }
    }
  }

}
