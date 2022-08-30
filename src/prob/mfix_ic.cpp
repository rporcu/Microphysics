#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>

#include <mfix_ic.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_regions.H>
#include <mfix_species.H>


using namespace amrex;


void
MFIXInitialConditions::Initialize (const MFIXRegions& regions,
                                   MFIXFluidPhase& fluid,
                                   const MFIXSolidsPhase& solids,
                                   MFIXDEM& dem,
                                   MFIXPIC& pic)
{
  amrex::ParmParse pp("ic");

  std::vector<std::string> input_regions;
  pp.queryarr("regions", input_regions);

  // Loop over ICs
  for (size_t icv=0; icv < input_regions.size(); icv++) {

    amrex::Real volfrac_total(0.0);

    IC_t new_ic(fluid);

    // Set the region for the initial condition.
    new_ic.region = regions.getRegion(input_regions[icv]);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(new_ic.region != nullptr, "Invalid ic region!");

    // Get fluid data.
    if (fluid.solve()) {

      std::string field = "ic."+input_regions[icv]+"."+fluid.names(0);
      amrex::ParmParse ppFluid(field.c_str());

      ppFluid.get("volfrac", new_ic.fluid.volfrac);
      volfrac_total += new_ic.fluid.volfrac;

      ppFluid.getarr("velocity", new_ic.fluid.velocity, 0, 3);

      if (fluid.solve_species()) {

        const int nspecies_g = fluid.nspecies();
        new_ic.fluid.species.resize(nspecies_g);

        std::string species_field = field+".species";
        amrex::ParmParse ppSpecies(species_field.c_str());

        // Auxiliary variable to check that species sum up to 1
        Real total_mass_fraction(0);

        for (int n(0); n < nspecies_g; n++) {
          // Get the name of the fluid species we want to get the IC
          std::string fluid_species = fluid.species_names(n);
          // Get the IC mass fraction for the current species
          ppSpecies.get(fluid_species.c_str(), new_ic.fluid.species[n].mass_fraction);
          total_mass_fraction += new_ic.fluid.species[n].mass_fraction;
        }

        // Sanity check that the input species mass fractions sum up to 1
        if (!(Math::abs(total_mass_fraction-1) < 1.e-15)) {
          std::string message = "Error: species ICs mass fractions in region "
            + input_regions[icv] + " sum up to " + std::to_string(total_mass_fraction) + "\n";

          amrex::Abort(message);
        }
      }

      // Get density
      new_ic.fluid.density_defined = ppFluid.query("density", new_ic.fluid.density);

      // Get temperature
      new_ic.fluid.temperature_defined = ppFluid.query("temperature", new_ic.fluid.temperature);

      // Get pressure
      new_ic.fluid.pressure_defined = ppFluid.query("pressure", new_ic.fluid.pressure);

      if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IncompressibleFluid) {

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(new_ic.fluid.density_defined,
            "Error: when the incompressible fluid constraint is selected, "
            "fluid initial conditions must include density");
      } else {

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!new_ic.fluid.density_defined,
            "Error: when the ideal gas EOS constraint is selected, "
            "fluid initial conditions must NOT include density");

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.p_therm_defined(),
            "Error: fluid thermodynamic pressure must be defined before "
            "computing IC density out of the ideal gas EOS");

        const auto& fluid_parms = fluid.parameters();
        Real MW_g(0.);

        if (fluid.isMixture()) {
          for (int n(0); n < fluid.nspecies(); n++) {
            MW_g += new_ic.fluid.species[n].mass_fraction / fluid_parms.get_MW_gk<RunOn::Host>(n);
          }

          MW_g = 1./MW_g;
        } else {
          MW_g = fluid_parms.get_MW_g<RunOn::Host>();
        }

        new_ic.fluid.density = (fluid.thermodynamic_pressure() * MW_g) /
                               (MFIXFluidPhase::R * new_ic.fluid.temperature);
      }
    }

    if (dem.solve() || pic.solve()) {

      // If we initialize particles with particle generator
      if (fluid.solve() && new_ic.fluid.volfrac < 1.0) {

        // Get the list of solids used in defining the IC region
        std::vector<std::string> solids_names;
        {
          std::string field = "ic."+input_regions[icv];
          amrex::ParmParse ppSolid(field.c_str());
          ppSolid.getarr("solids", solids_names);
          ppSolid.get("packing", new_ic.packing);
        }

        for (size_t lcs(0); lcs < solids_names.size(); ++lcs) {

          SOLIDS_t new_solid;

          std::string field = "ic."+input_regions[icv]+"."+solids_names[lcs];
          amrex::ParmParse ppSolid(field.c_str());

          new_solid.name = solids_names[lcs];
          new_solid.phase = solids.name_to_phase(solids_names[lcs]);

          ppSolid.get("volfrac", new_solid.volfrac);
          volfrac_total += new_ic.fluid.volfrac;

          ppSolid.getarr("velocity", new_solid.velocity, 0, 3);

          if (fluid.solve_enthalpy() || solids.solve_enthalpy()) {
            ppSolid.query("temperature", new_solid.temperature);
          }

          new_solid.statwt = 1.0;
          ppSolid.query("statwt", new_solid.statwt);

          { // Get information about diameter distribution.
            int err = new_solid.diameter.set(field, "diameter");
            if(err) {
              std::string message = " Error: IC region " + input_regions[icv]
                + " has invalid diameter parameters for solids " + solids_names[lcs] + "\n";
              amrex::Print() << message;
              amrex::Abort(message);
            }
          }

          { // Get information about density distribution
            int err = new_solid.density.set(field, "density");
            if(err) {
              std::string message = " Error: BC region " + input_regions[icv]
                + " has invalid density parameters for solids " + solids_names[lcs] + "\n";
              amrex::Print() << message;
              amrex::Abort(message);
            }
          }

          if (solids.solve_species())
          {
            std::string species_field = field+".species";
            amrex::ParmParse ppSpecies(species_field.c_str());

            new_solid.species.resize(solids.nspecies());

            amrex::Real total_mass_fraction(0);

            for (int n(0); n < solids.nspecies(); n++) {
              std::string current_species = solids.species_names(n);
              ppSpecies.get(current_species.c_str(), new_solid.species[n].mass_fraction);

              total_mass_fraction += new_solid.species[n].mass_fraction;
            }

            // Sanity check that the input species mass fractions sum up to 1
            if (!(amrex::Math::abs(total_mass_fraction-1) < 1.e-15)) {
              std::string message = "Error: SOLID type " + solids_names[lcs]
                + " species ICs mass fractions in region " + input_regions[icv]
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
        std::vector<std::string> solids_names(0);
        {
          std::string field = "ic."+input_regions[icv];
          amrex::ParmParse ppSolid(field.c_str());
          ppSolid.queryarr("solids", solids_names);
        }

        for(size_t lcs(0); lcs < solids_names.size(); ++ lcs) {

          SOLIDS_t new_solid;

          std::string field = "ic."+input_regions[icv]+"."+solids_names[lcs];
          amrex::ParmParse ppSolid(field.c_str());

          new_solid.name = solids_names[lcs];
          new_solid.phase = solids.name_to_phase(solids_names[lcs]);

          if (fluid.solve_enthalpy() || solids.solve_enthalpy()) {
            ppSolid.get("temperature", new_solid.temperature);
          }

          if (solids.solve_species()) {

            std::string species_field = field+".species";
            amrex::ParmParse ppSpecies(species_field.c_str());

            new_solid.species.resize(solids.nspecies());

            for (int n(0); n < solids.nspecies(); n++) {
              std::string current_species = solids.species_names(n);
              ppSpecies.query(current_species.c_str(), new_solid.species[n].mass_fraction);
            }
          }

          new_ic.solids.push_back(new_solid);
        }
      }
    }

    m_ic.push_back(new_ic);
  }


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

    if(fluid.solve()){

      amrex::Print() << std::endl;
      amrex::Print() << "   Fluid:     volfrac: " << ic[icv].fluid.volfrac     << std::endl;
      amrex::Print() << "              density: " << ic[icv].fluid.density     << std::endl;
      amrex::Print() << "             pressure: " << ic[icv].fluid.pressure    << std::endl;
      amrex::Print() << "          temperature: " << ic[icv].fluid.temperature << std::endl;
      amrex::Print() << "             velocity: ";
      amrex::Print() << ic[icv].fluid.velocity[0] << "  ";
      amrex::Print() << ic[icv].fluid.velocity[1] << "  ";
      amrex::Print() << ic[icv].fluid.velocity[2] << std::endl;
      amrex::Print() << std::endl;
    }


    if((dem.solve() || m_pic.solve()) && ic[icv].solids.size()>0){

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

        amrex::Print() << " density distribution: " << ic[icv].solids[lcs].density.distribution << std::endl;
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
