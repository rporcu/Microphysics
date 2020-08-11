#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>

#include <mfix_ic_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_regions_parms.H>
#include <mfix_species_parms.H>

using namespace amrex;

namespace IC
{

  amrex::Vector<IC_t> ic;

  void Initialize ()
  {

    amrex::ParmParse pp("ic");

    std::vector<std::string> regions;
    pp.queryarr("regions", regions);

    // Query if advect_enthalpy so we check if ic temperature inputs are correct
    amrex::ParmParse ppMFIX("mfix");
    int advect_enthalpy(0);

    ppMFIX.query("advect_enthalpy", advect_enthalpy);

    // Loop over ICs
    for(size_t icv=0; icv < regions.size(); icv++){

      amrex::Real volfrac_total(0.0);

      IC_t new_ic;

      // Set the region for the initial condition.
      new_ic.region = REGIONS::getRegion(regions[icv]);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE( new_ic.region != NULL, "Invalid ic region!");

      // Get fluid data.

      if(FLUID::solve) {

        std::string field = "ic."+regions[icv]+"."+FLUID::name;
        amrex::ParmParse ppFluid(field.c_str());

        ppFluid.get("volfrac", new_ic.fluid.volfrac);
        volfrac_total += new_ic.fluid.volfrac;

        ppFluid.getarr("velocity", new_ic.fluid.velocity, 0, 3);

        if (advect_enthalpy) {
          ppFluid.get("temperature", new_ic.fluid.temperature); 
        }

        new_ic.fluid.pressure_defined = ppFluid.query("pressure", new_ic.fluid.pressure);

        if (FLUID::solve_species) {

          std::string species_field = field+".species";
          amrex::ParmParse ppSpecies(species_field.c_str());

          for (int n(0); n < FLUID::nspecies_g; n++) {
            std::string fluid_specie = FLUID::species_g[n];
            ppSpecies.query(fluid_specie.c_str(), new_ic.fluid.species.mass_fractions[n]);
          }
        }
      }

      if (DEM::solve or PIC::solve)
      {
        // If we initialize particles with particle generator
        if (FLUID::solve and new_ic.fluid.volfrac < 1.0)
        {
          // Get the list of solids used in defining the IC region
          std::vector<std::string> solids_types;
          {
            std::string field = "ic."+regions[icv];
            amrex::ParmParse ppSolid(field.c_str());
            ppSolid.getarr("solids", solids_types);
            ppSolid.get("packing", new_ic.packing);
          }

          for(size_t lcs(0); lcs < solids_types.size(); ++ lcs) {

            SOLIDS::SOLIDS_t new_solid;

            std::string field = "ic."+regions[icv]+"."+solids_types[lcs];
            amrex::ParmParse ppSolid(field.c_str());

            new_solid.name = solids_types[lcs];

            ppSolid.get("volfrac", new_solid.volfrac);
            volfrac_total += new_ic.fluid.volfrac;

            ppSolid.getarr("velocity", new_solid.velocity, 0, 3);

            if (advect_enthalpy) {
              ppSolid.query("temperature", new_solid.temperature); 
            }

            new_solid.statwt = 1.0;
            ppSolid.query("statwt", new_solid.statwt);

            // Get information about diameter distribution.
            ppSolid.get("diameter", new_solid.diameter.distribution);

            std::string dp_field = "ic."+regions[icv]+"."+solids_types[lcs]+".diameter";
            amrex::ParmParse ppSolidDp(dp_field.c_str());

            if( new_solid.diameter.distribution == "constant") {
              ppSolidDp.get("constant", new_solid.diameter.mean);

            } else { // This could probably be an else-if to better catch errors
              ppSolidDp.get("mean", new_solid.diameter.mean);
              ppSolidDp.get("std" , new_solid.diameter.std);
              ppSolidDp.get("min" , new_solid.diameter.min);
              ppSolidDp.get("max" , new_solid.diameter.max);
            }

            // Get information about density distribution.
            ppSolid.get("density", new_solid.density.distribution);

            std::string roh_field = "ic."+regions[icv]+"."+solids_types[lcs]+".density";
            amrex::ParmParse ppSolidRho(roh_field.c_str());

            if( new_solid.diameter.distribution == "constant") {
              ppSolidRho.get("constant", new_solid.density.mean);

            } else { // This could probably be an else-if to better catch errors
              ppSolidRho.get("mean", new_solid.density.mean);
              ppSolidRho.get("std" , new_solid.density.std );
              ppSolidRho.get("min" , new_solid.density.min );
              ppSolidRho.get("max" , new_solid.density.max );
            }

            if (DEM::solve_species /*TODO or PIC::solve_species*/) {

              std::string species_field = field+".species";
              amrex::ParmParse ppSpecies(species_field.c_str());

              for (int n(0); n < DEM::nspecies_dem; n++) {
                std::string dem_specie = DEM::species_dem[n];
                ppSpecies.query(dem_specie.c_str(), new_solid.species.mass_fractions[n]);
              }
            }

            new_ic.solids.push_back(new_solid);
          }
        }
        // If we initialize particles through particle_input.dat
        else {
          // Get the list of solids used in defining the IC region
          std::vector<std::string> solids_types(0);
          {
            std::string field = "ic."+regions[icv];
            amrex::ParmParse ppSolid(field.c_str());
            ppSolid.queryarr("solids", solids_types);
          }

          for(size_t lcs(0); lcs < solids_types.size(); ++ lcs) {

            SOLIDS::SOLIDS_t new_solid;

            std::string field = "ic."+regions[icv]+"."+solids_types[lcs];
            amrex::ParmParse ppSolid(field.c_str());

            new_solid.name = solids_types[lcs];

            if (advect_enthalpy) {
              ppSolid.get("temperature", new_solid.temperature); 
            }

            if (DEM::solve_species /*TODO or PIC::solve_species*/) {

              std::string species_field = field+".species";
              amrex::ParmParse ppSpecies(species_field.c_str());

              for (int n(0); n < DEM::nspecies_dem; n++) {
                std::string dem_specie = DEM::species_dem[n];
                ppSpecies.query(dem_specie.c_str(), new_solid.species.mass_fractions[n]);
              }
            }

            new_ic.solids.push_back(new_solid);
          }
        }
      }

      IC::ic.push_back(new_ic);

    }


#if defined(AMREX_DEBUG)
    //Dump out what we read for debugging!
    for (int icv(0); icv<ic.size(); icv++){

      amrex::Print() << std::endl << "Summarizing IC regions:\n" << std::endl;
      amrex::Print() << " IC: " << icv << std::endl << std::endl <<
        "      lo: " << ic[icv].region->lo(0) << "  "
                     << ic[icv].region->lo(1) << "  "
                     << ic[icv].region->lo(2) << std::endl <<
        "      hi: " << ic[icv].region->hi(0) << "  "
                     << ic[icv].region->hi(1) << "  "
                     << ic[icv].region->hi(2) << std::endl;

      if(FLUID::solve){

        amrex::Print() << std::endl;
        amrex::Print() << "   Fluid:     volfrac: " << ic[icv].fluid.volfrac     << std::endl;
        amrex::Print() << "             pressure: " << ic[icv].fluid.pressure    << std::endl;
        amrex::Print() << "          temperature: " << ic[icv].fluid.temperature << std::endl;
        amrex::Print() << "             velocity: ";
        amrex::Print() << ic[icv].fluid.velocity[0] << "  ";
        amrex::Print() << ic[icv].fluid.velocity[1] << "  ";
        amrex::Print() << ic[icv].fluid.velocity[2] << std::endl;
        amrex::Print() << std::endl;
      }


      if((DEM::solve or PIC::solve) and ic[icv].solids.size()>0){

        amrex::Print() << "       Solids packing: " << ic[icv].packing << std::endl;

        for(int lcs(0); lcs<ic[icv].solids.size(); ++lcs){
          amrex::Print() << std::endl;
          amrex::Print() << "   Solid:       index: " << lcs << std::endl;
          amrex::Print() << "              volfrac: " << ic[icv].solids[lcs].volfrac     << std::endl;
          amrex::Print() << "          temperature: " << ic[icv].solids[lcs].temperature << std::endl;
          amrex::Print() << "             velocity: ";
          amrex::Print() << ic[icv].solids[lcs].velocity[0] << "  ";
          amrex::Print() << ic[icv].solids[lcs].velocity[1] << "  ";
          amrex::Print() << ic[icv].solids[lcs].velocity[2] << std::endl;

          amrex::Print() << "diameter distribution: " << ic[icv].solids[lcs].diameter.distribution << std::endl;
          if( ic[icv].solids[lcs].density.distribution == "constant" ){
            amrex::Print() << "                value: " << ic[icv].solids[lcs].diameter.mean << std::endl;
          } else {
            amrex::Print() << "                 mean: " << ic[icv].solids[lcs].diameter.mean << std::endl;
            amrex::Print() << "                  std: " << ic[icv].solids[lcs].diameter.std  << std::endl;
            amrex::Print() << "                  min: " << ic[icv].solids[lcs].diameter.min  << std::endl;
            amrex::Print() << "                  max: " << ic[icv].solids[lcs].diameter.max  << std::endl;
          }

          amrex::Print() << " desnity distribution: " << ic[icv].solids[lcs].density.distribution << std::endl;
          if( ic[icv].solids[lcs].density.distribution == "constant" ){
            amrex::Print() << "                value: " << ic[icv].solids[lcs].density.mean << std::endl;
          } else {
            amrex::Print() << "                 mean: " << ic[icv].solids[lcs].density.mean << std::endl;
            amrex::Print() << "                  std: " << ic[icv].solids[lcs].density.std  << std::endl;
            amrex::Print() << "                  min: " << ic[icv].solids[lcs].density.min  << std::endl;
            amrex::Print() << "                  max: " << ic[icv].solids[lcs].density.max  << std::endl;
          }

          amrex::Print() << std::endl;

        }

      }

    }
#endif
  }
}
