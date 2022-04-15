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

  void Initialize (const FluidPhase& fluid,
                   const SolidsPhase& solids)
  {

    amrex::ParmParse pp("ic");

    std::vector<std::string> regions;
    pp.queryarr("regions", regions);

    // Query if advect_enthalpy so we check if ic temperature inputs are correct
    amrex::ParmParse ppMFIX("mfix");
    int advect_enthalpy(0);

    ppMFIX.query("advect_enthalpy", advect_enthalpy);

    // Loop over ICs
    for(size_t icv=0; icv < regions.size(); icv++) {

      amrex::Real volfrac_total(0.0);

      IC_t new_ic(fluid);

      // Set the region for the initial condition.
      new_ic.region = REGIONS::getRegion(regions[icv]);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE( new_ic.region != NULL, "Invalid ic region!");

      // Get fluid data.

      if(fluid.solve) {

        std::string field = "ic."+regions[icv]+"."+fluid.name;
        amrex::ParmParse ppFluid(field.c_str());

        ppFluid.get("volfrac", new_ic.fluid.volfrac);
        volfrac_total += new_ic.fluid.volfrac;

        ppFluid.getarr("velocity", new_ic.fluid.velocity, 0, 3);

        if (fluid.constraint_type == ConstraintType::IncompressibleFluid) {
          ppFluid.get("density", new_ic.fluid.density);
        }

        if (fluid.solve_enthalpy && fluid.constraint_type == ConstraintType::IncompressibleFluid) {
          ppFluid.get("temperature", new_ic.fluid.temperature); 
        }

        new_ic.fluid.pressure_defined = ppFluid.query("pressure", new_ic.fluid.pressure);

        if (fluid.solve_species) {

          const int nspecies_g = fluid.nspecies;
          new_ic.fluid.species.resize(nspecies_g);

          std::string species_field = field+".species";
          amrex::ParmParse ppSpecies(species_field.c_str());

          // Auxiliary variable to check that species sum up to 1
          Real total_mass_fraction(0);

          for (int n(0); n < nspecies_g; n++) {
            // Get the name of the fluid species we want to get the IC
            std::string fluid_species = fluid.species[n];
            // Get the IC mass fraction for the current species
            ppSpecies.get(fluid_species.c_str(), new_ic.fluid.species[n].mass_fraction);
            total_mass_fraction += new_ic.fluid.species[n].mass_fraction;
          }

          // Sanity check that the input species mass fractions sum up to 1
          if (!(Math::abs(total_mass_fraction-1) < 1.e-15)) {
            std::string message = "Error: species ICs mass fractions in region "
              + regions[icv] + " sum up to " + std::to_string(total_mass_fraction) + "\n";

            amrex::Abort(message);
          }
        }

        if (fluid.constraint_type == ConstraintType::IdealGasOpenSystem ||
            fluid.constraint_type == ConstraintType::IdealGasClosedSystem) {

          // Get density
          const int density_defined = ppFluid.query("density", new_ic.fluid.density);
          new_ic.fluid.density_defined = density_defined;

          // Get temperature
          const int temperature_defined = ppFluid.query("temperature", new_ic.fluid.temperature);
          new_ic.fluid.temperature_defined = temperature_defined;

          // Get thermodynamic pressure
          const int thermodynamic_p_defined = ppFluid.query("thermodynamic_pressure",
                                                            new_ic.fluid.thermodynamic_pressure);
          new_ic.fluid.thermodynamic_pressure_defined = thermodynamic_p_defined;

          const int sum_defined = density_defined + temperature_defined + thermodynamic_p_defined;

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sum_defined == 2,
              "Initial conditions inputs must provide exactly two quantities"
              " among fluid temperature, density, and thermodynamic pressure");

          auto& fluid_parms = *fluid.parameters;
          Real MW_g(0.);

          for (int n(0); n < fluid.nspecies; n++) {
            MW_g += new_ic.fluid.species[n].mass_fraction / fluid_parms.get_MW_gk<RunOn::Host>(n);
          }

          MW_g = 1./MW_g;

          if (!thermodynamic_p_defined) {

            new_ic.fluid.thermodynamic_pressure = (new_ic.fluid.density *
                fluid_parms.R * new_ic.fluid.temperature) / MW_g;

          } else if (!temperature_defined) {

            new_ic.fluid.temperature = (new_ic.fluid.thermodynamic_pressure *
                MW_g) / (fluid_parms.R * new_ic.fluid.density);

          } else if (!density_defined) {

            new_ic.fluid.density = (new_ic.fluid.thermodynamic_pressure *
                MW_g) / (fluid_parms.R * new_ic.fluid.temperature);

          } else {
            amrex::Abort("How did we arrive here?");
          }
        }
      }

      if (DEM::solve || PIC::solve) {

        // If we initialize particles with particle generator
        if (fluid.solve && new_ic.fluid.volfrac < 1.0) {

          // Get the list of solids used in defining the IC region
          std::vector<std::string> solids_types;
          {
            std::string field = "ic."+regions[icv];
            amrex::ParmParse ppSolid(field.c_str());
            ppSolid.getarr("solids", solids_types);
            ppSolid.get("packing", new_ic.packing);
          }

          for(size_t lcs(0); lcs < solids_types.size(); ++ lcs) {

            SolidsPhase::SOLIDS_t new_solid;

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

            if (solids.solve_species)
            {
              std::string species_field = field+".species";
              amrex::ParmParse ppSpecies(species_field.c_str());

              // TODO: check this when nb of solids > 1
              new_solid.species.resize(solids.nspecies);

              amrex::Real total_mass_fraction(0);

              for (int n(0); n < solids.nspecies; n++) {
                std::string current_species = solids.species[n];
                ppSpecies.get(current_species.c_str(), new_solid.species[n].mass_fraction);

                total_mass_fraction += new_solid.species[n].mass_fraction;
              }

              // Sanity check that the input species mass fractions sum up to 1
              if (!(amrex::Math::abs(total_mass_fraction-1) < 1.e-15)) {
                std::string message = "Error: SOLID type " + solids_types[lcs]
                  + " species ICs mass fractions in region " + regions[icv]
                  + " sum up to " + std::to_string(total_mass_fraction) + "\n";

                amrex::Abort(message);
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

            SolidsPhase::SOLIDS_t new_solid;

            std::string field = "ic."+regions[icv]+"."+solids_types[lcs];
            amrex::ParmParse ppSolid(field.c_str());

            new_solid.name = solids_types[lcs];

            if (advect_enthalpy) {
              ppSolid.get("temperature", new_solid.temperature); 
            }

            if (solids.solve_species) {

              std::string species_field = field+".species";
              amrex::ParmParse ppSpecies(species_field.c_str());

              // TODO: check this when nb of solids > 1
              new_solid.species.resize(solids.nspecies);

              for (int n(0); n < solids.nspecies; n++) {
                std::string current_species = solids.species[n];
                ppSpecies.query(current_species.c_str(), new_solid.species[n].mass_fraction);
              }
            }

            new_ic.solids.push_back(new_solid);
          }
        }
      }

      IC::ic.push_back(new_ic);

    }


//    // This is a check that the initial thermodynamic pressure is uniform in the
//    // whole domain
//    if (fluid.constraint_type == ConstraintType::IdealGasOpenSystem ||
//        fluid.constraint_type == ConstraintType::IdealGasClosedSystem) {
//
//      for (int icv(1); icv < IC::ic.size(); ++icv) {
//        const Real diff = std::abs(IC::ic[0].fluid.pressure -
//                                   IC::ic[icv].fluid.pressure);
//
//        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(diff < 1.e-15,
//            "ICs for thermodynamic pressure are not uniform in space");
//      }
//    }


#if 0
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

      if(fluid.solve){

        amrex::Print() << std::endl;
        amrex::Print() << "   Fluid:     volfrac: " << ic[icv].fluid.volfrac     << std::endl;
        amrex::Print() << "              density: " << ic[icv].fluid.density     << std::endl;
        amrex::Print() << "termodynamic_pressure: " << ic[icv].fluid.thermodynamic_pressure    << std::endl;
        amrex::Print() << "             pressure: " << ic[icv].fluid.pressure    << std::endl;
        amrex::Print() << "          temperature: " << ic[icv].fluid.temperature << std::endl;
        amrex::Print() << "             velocity: ";
        amrex::Print() << ic[icv].fluid.velocity[0] << "  ";
        amrex::Print() << ic[icv].fluid.velocity[1] << "  ";
        amrex::Print() << ic[icv].fluid.velocity[2] << std::endl;
        amrex::Print() << std::endl;
      }


      if((DEM::solve || PIC::solve) && ic[icv].solids.size()>0){

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
