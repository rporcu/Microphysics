#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>

#include <MFIX_BC_Parms.H>
#include <MFIX_REGIONS_Parms.H>

namespace BC
{

  // Direction of pressure drop (0:x, 1:y, 2:z)
  int delp_dir = -1;

  // Specified constant gas density
  amrex::Real delp[3];

  amrex::Vector<BC_t> bc;

  void Initialize ()
  {

    // Initialize the periodic pressure drop to zero in all directions
    for (int dir=0; dir < 3; dir ++ )
      delp[dir] = 0.0;

    amrex::ParmParse pp("bc");

    amrex::Real delp_in(0.0);
    pp.query("delp", delp_in);

    if( delp_in != 0.0) {

      pp.query("delp_dir", delp_dir);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE( 0 <= delp_dir && delp_dir <= 2,
         "Direction of pressure drop (bc.delp_dir) must be specified.");

      delp[delp_dir] = delp_in;

    }

    std::vector<std::string> regions;
    pp.queryarr("regions", regions);

    // Loop over BCs
    for(int bcv=0; bcv<regions.size(); bcv++){

      amrex::Real volfrac_total(0.0);

      BC_t new_bc;

      // Set the region for the initial condition.
      new_bc.region = REGIONS::getRegion(regions[bcv]);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE( new_bc.region != NULL, "Invalid bc region!");


      // Get the BC type (MI/PI/PO...)
      pp.get(regions[bcv].c_str(),new_bc.type);

      // Get fluid data.

      if(FLUID::solve) {

        std::string field = "bc."+regions[bcv]+"."+FLUID::name;
        amrex::ParmParse ppFluid(field.c_str());

        ppFluid.query("volfrac", new_bc.fluid.volfrac);
        volfrac_total += new_bc.fluid.volfrac;

        ppFluid.queryarr("velocity", new_bc.fluid.velocity, 0, 3);
        ppFluid.query("temperature", new_bc.fluid.temperature);

        new_bc.fluid.pressure_defined = ppFluid.query("pressure", new_bc.fluid.pressure);
      }

      if(DEM::solve) {

        // Get the list of solids used in defining the BC region
        std::vector<std::string> solids_types;
        {
          std::string field = "bc."+regions[bcv];
          amrex::ParmParse ppSolid(field.c_str());
          ppSolid.queryarr("solids", solids_types);
        }

        for(int lcs(0); lcs < solids_types.size(); ++ lcs){

          DEM::DEM_t new_solid;

          std::string field = "bc."+regions[bcv]+"."+solids_types[lcs];
          amrex::ParmParse ppSolid(field.c_str());

          amrex::Real temp_volfrac(0.0);
          ppSolid.get("volfrac", new_solid.volfrac);
          volfrac_total += new_bc.fluid.volfrac;

          ppSolid.getarr("velocity", new_solid.velocity, 0, 3);
          ppSolid.query("temperature", new_solid.temperature);

          // Get information about diameter distribution.
          ppSolid.get("diameter", new_solid.diameter.distribution);

          std::string dp_field = "bc."+regions[bcv]+"."+solids_types[lcs]+".diameter";
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

          std::string roh_field = "bc."+regions[bcv]+"."+solids_types[lcs]+".density";
          amrex::ParmParse ppSolidRho(roh_field.c_str());

          if( new_solid.diameter.distribution == "constant") {
            ppSolidRho.get("constant", new_solid.density.mean);

          } else { // This could probably be an else-if to better catch errors
            ppSolidRho.get("mean", new_solid.density.mean);
            ppSolidRho.get("std" , new_solid.density.std );
            ppSolidRho.get("min" , new_solid.density.min );
            ppSolidRho.get("max" , new_solid.density.max );
          }

          new_bc.solids.push_back(new_solid);
        }
      }

      BC::bc.push_back(new_bc);

    }


    //Dump out what we read for debugging!
    for (int bcv(0); bcv<bc.size(); bcv++){


      amrex::Print() << std::endl << "Summarizing BC regions:\n" << std::endl;
      amrex::Print() << " BC: " << bcv << "    Type: " << bc[bcv].type << std::endl << std::endl <<
        "      lo: " << bc[bcv].region->lo(0) << "  "
                     << bc[bcv].region->lo(1) << "  "
                     << bc[bcv].region->lo(2) << std::endl <<
        "      hi: " << bc[bcv].region->hi(0) << "  "
                     << bc[bcv].region->hi(1) << "  "
                     << bc[bcv].region->hi(2) << std::endl;

      if(FLUID::solve){

        amrex::Print() << std::endl;
        amrex::Print() << "   Fluid:     volfrac: " << bc[bcv].fluid.volfrac     << std::endl;
        amrex::Print() << "             pressure: " << bc[bcv].fluid.pressure    << std::endl;
        amrex::Print() << "          temperature: " << bc[bcv].fluid.temperature << std::endl;
        if( bc[bcv].fluid.velocity.size() > 0 ) {
          amrex::Print() << "             velocity: ";
          amrex::Print() << bc[bcv].fluid.velocity[0] << "  ";
          amrex::Print() << bc[bcv].fluid.velocity[1] << "  ";
          amrex::Print() << bc[bcv].fluid.velocity[2] << std::endl;
        }
        amrex::Print() << std::endl;
      }


      if(DEM::solve){

        for(int lcs(0); lcs<bc[bcv].solids.size(); ++lcs){
          amrex::Print() << std::endl;
          amrex::Print() << "   Solid:       index: " << lcs << std::endl;
          amrex::Print() << "              volfrac: " << bc[bcv].solids[lcs].volfrac     << std::endl;
          amrex::Print() << "          temperature: " << bc[bcv].solids[lcs].temperature << std::endl;
          if(bc[bcv].solids[lcs].velocity.size() > 0){
            amrex::Print() << "             velocity: ";
            amrex::Print() << bc[bcv].solids[lcs].velocity[0] << "  ";
            amrex::Print() << bc[bcv].solids[lcs].velocity[1] << "  ";
            amrex::Print() << bc[bcv].solids[lcs].velocity[2] << std::endl;
          }

          amrex::Print() << "diameter distribution: " << bc[bcv].solids[lcs].diameter.distribution << std::endl;
          if( bc[bcv].solids[lcs].density.distribution == "constant" ){
            amrex::Print() << "                value: " << bc[bcv].solids[lcs].diameter.mean << std::endl;
          } else {
            amrex::Print() << "                 mean: " << bc[bcv].solids[lcs].diameter.mean << std::endl;
            amrex::Print() << "                  std: " << bc[bcv].solids[lcs].diameter.std  << std::endl;
            amrex::Print() << "                  min: " << bc[bcv].solids[lcs].diameter.min  << std::endl;
            amrex::Print() << "                  max: " << bc[bcv].solids[lcs].diameter.max  << std::endl;
          }

          amrex::Print() << " desnity distribution: " << bc[bcv].solids[lcs].density.distribution << std::endl;
          if( bc[bcv].solids[lcs].density.distribution == "constant" ){
            amrex::Print() << "                value: " << bc[bcv].solids[lcs].density.mean << std::endl;
          } else {
            amrex::Print() << "                 mean: " << bc[bcv].solids[lcs].density.mean << std::endl;
            amrex::Print() << "                  std: " << bc[bcv].solids[lcs].density.std  << std::endl;
            amrex::Print() << "                  min: " << bc[bcv].solids[lcs].density.min  << std::endl;
            amrex::Print() << "                  max: " << bc[bcv].solids[lcs].density.max  << std::endl;
          }

          amrex::Print() << std::endl;

        }

      }

    }  // END loop over BCs to print output

    exit(0);

  }// END Initialize

}
