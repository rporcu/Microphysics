#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>

#include <MFIX_IC_Parms.H>
#include <MFIX_REGIONS_Parms.H>

namespace IC
{

  amrex::Vector<IC_t> ic;

  void Initialize ()
  {

    amrex::ParmParse pp("ic");

    std::vector<std::string> regions;
    pp.queryarr("regions", regions);

    // Loop over ICs
    for(int icv=0; icv<regions.size(); icv++){

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
        ppFluid.query("temperature", new_ic.fluid.temperature);

        new_ic.fluid.pressure_defined = ppFluid.query("pressure", new_ic.fluid.pressure);
      }

      if(DEM::solve && new_ic.fluid.volfrac < 1.0) {

        // Get the list of solids used in defining the IC region
        std::vector<std::string> solids_types;
        {
          std::string field = "ic."+regions[icv];
          amrex::ParmParse ppSolid(field.c_str());
          ppSolid.getarr("solids", solids_types);
          ppSolid.get("packing", new_ic.packing);
        }

        for(int lcs(0); lcs < solids_types.size(); ++ lcs){

          DEM::DEM_t new_solid;

          std::string field = "ic."+regions[icv]+"."+solids_types[lcs];
          amrex::ParmParse ppSolid(field.c_str());

          //amrex::Real temp_volfrac(0.0);
          ppSolid.get("volfrac", new_solid.volfrac);
          volfrac_total += new_ic.fluid.volfrac;

          ppSolid.getarr("velocity", new_solid.velocity, 0, 3);
          ppSolid.query("temperature", new_solid.temperature);

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

          new_ic.solids.push_back(new_solid);
        }
      }

      IC::ic.push_back(new_ic);

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


      if(DEM::solve && ic[icv].solids.size()>0){

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
