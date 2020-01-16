#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_Geometry.H>

#include <AMReX_ParmParse.H>

#include <MFIX_BC_Parms.H>
#include <MFIX_REGIONS_Parms.H>
#include <MFIX_BcList.H>

namespace BC
{

  // Direction of pressure drop (0:x, 1:y, 2:z)
  int delp_dir = -1;
  amrex::Real delp[3];

  int domain_bc[6];

  // EB planes for level-set creation
  amrex::Vector<amrex::EB2::PlaneIF> flow_planes;
  amrex::Vector<amrex::EB2::PlaneIF> wall_planes;

  // Lists of BCs applied to the domain extent
  amrex::Vector<int> bc_xlo, bc_xhi;
  amrex::Vector<int> bc_ylo, bc_yhi;
  amrex::Vector<int> bc_zlo, bc_zhi;

  std::array<amrex::LinOpBCType,3> ppe_lobc;
  std::array<amrex::LinOpBCType,3> ppe_hibc;

  std::array<amrex::LinOpBCType,3> diff_vel_lobc;
  std::array<amrex::LinOpBCType,3> diff_vel_hibc;

  std::array<amrex::LinOpBCType,3> diff_scal_lobc;
  std::array<amrex::LinOpBCType,3> diff_scal_hibc;


  // Data structure storing individual BC information
  amrex::Vector<BC_t> bc;


  void Initialize (amrex::Geometry& geom)
  {

    BcList bc_mask;

    // Integer ids for BCs
    const int pout_ = bc_mask.get_pout();
    const int pinf_ = bc_mask.get_pinf();
    const int minf_ = bc_mask.get_minf();
    const int nsw_  = bc_mask.get_nsw();

    // Set flag to keep particles from leaving unless periodic.
    for (int dir(0); dir<3; ++dir) {
      if ( geom.isPeriodic(dir)) {
        domain_bc[2*dir  ] = 0;
        domain_bc[2*dir+1] = 0;
      } else {
        domain_bc[2*dir  ] = 1;
        domain_bc[2*dir+1] = 1;
      }
    }


    // Default all sides of the domain to Neumann
    for (int dir=0; dir < 3; dir ++ ){
      if( geom.isPeriodic(dir)){

        ppe_lobc[dir] = amrex::LinOpBCType::Periodic;
        ppe_hibc[dir] = amrex::LinOpBCType::Periodic;

        diff_vel_lobc[dir] = amrex::LinOpBCType::Periodic;
        diff_vel_hibc[dir] = amrex::LinOpBCType::Periodic;

        diff_scal_lobc[dir] = amrex::LinOpBCType::Periodic;
        diff_scal_hibc[dir] = amrex::LinOpBCType::Periodic;

      } else {

        ppe_lobc[dir] = amrex::LinOpBCType::Neumann;
        ppe_hibc[dir] = amrex::LinOpBCType::Neumann;

        diff_vel_lobc[dir] = amrex::LinOpBCType::Dirichlet;
        diff_vel_hibc[dir] = amrex::LinOpBCType::Dirichlet;

        diff_scal_lobc[dir] = amrex::LinOpBCType::Dirichlet;
        diff_scal_hibc[dir] = amrex::LinOpBCType::Dirichlet;

      }
    }


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

    const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

    // Flag to see if particles 'see' a wall at pressure outflows
    int po_noParOut = 0; // default behavior for PO's -- letting particles exit the domain
    pp.query("po_no_par_out", po_noParOut);

    // Get the list of region names used to define BCs
    std::vector<std::string> regions;
    pp.queryarr("regions", regions);

    // Loop over BCs
    for(int bcv=0; bcv<regions.size(); bcv++){

      amrex::Real volfrac_total(0.0);

      BC_t new_bc;


      // Set the region for the initial condition.
      new_bc.region = REGIONS::getRegion(regions[bcv]);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE( new_bc.region != NULL, "Invalid bc region!");


      // Get the BC type (MI/PI/PO/NSW)
      std::string bc_type;
      pp.get(regions[bcv].c_str(),bc_type);

      // Convert the input string into the integers
      if( bc_type == "mi") {
        new_bc.type = minf_;
      } else if( bc_type == "pi") {
        new_bc.type = pinf_;
      } else if( bc_type == "po") {
        new_bc.type = pout_;
      } else if (bc_type == "nsw") {
        new_bc.type = nsw_;
      } else {
        amrex::Abort("Unknown BC type inputs file. Fix it.");
      }

      const auto plo = geom.ProbLoArray();
      const auto phi = geom.ProbHiArray();

      amrex::RealArray normal = {0.0, 0.0, 0.0};
      amrex::RealArray point;
      point[0] = new_bc.region->lo(0) + 0.5*(new_bc.region->hi(0)-new_bc.region->lo(0));
      point[1] = new_bc.region->lo(1) + 0.5*(new_bc.region->hi(1)-new_bc.region->lo(1));
      point[2] = new_bc.region->lo(2) + 0.5*(new_bc.region->hi(2)-new_bc.region->lo(2));

      std::string field_bcRegion = "bc."+regions[bcv];
      amrex::ParmParse ppRegion(field_bcRegion.c_str());

      int dir_int = -1;

      if(new_bc.type == nsw_ ){

        // Walls need a normal because they might not be on a domain extent.

        amrex::Vector<amrex::Real> normal_in(3);
        ppRegion.getarr("normal", normal_in, 0, 3);

        normal={normal_in[0], normal_in[1], normal_in[2]};

      } else {

        // This covers mass inflows (mi), pressure inflows (pi), and
        // pressure outflows (po). These all must align with a domain
        // extent.

        int sum_same_loc(0);
        for (int dir(0); dir<3; ++dir){

          int same_loc = std::abs(new_bc.region->lo(dir) - new_bc.region->hi(dir)) < tolerance ? 1 : 0;

          if (same_loc ){

            if ( std::abs(new_bc.region->lo(dir) - plo[dir] ) < tolerance ){

              point[dir] = plo[dir]+1.0e-15;
              normal[dir] =  1.0;

              dir_int = 2*dir;

              if( new_bc.type == pinf_ || new_bc.type == pout_){
                ppe_lobc[dir] = amrex::LinOpBCType::Dirichlet;
                diff_vel_lobc[dir] = amrex::LinOpBCType::Neumann;
                diff_scal_lobc[dir] = amrex::LinOpBCType::Neumann;
              }


            } else if ( std::abs( new_bc.region->hi(dir) - phi[dir] ) < tolerance ){

              point[dir] = phi[dir]-1.0e-15;
              normal[dir] = -1.0;

              dir_int = 2*dir+1;

              if( new_bc.type == pinf_ || new_bc.type == pout_){
                ppe_hibc[dir] = amrex::LinOpBCType::Dirichlet;
                diff_vel_hibc[dir] = amrex::LinOpBCType::Neumann;
                diff_scal_hibc[dir] = amrex::LinOpBCType::Neumann;
              }

            } else {
              amrex::Print() << "Flow BCs must be located on domain extents!" << std::endl;
              amrex::Print() << "BC Name: " << regions[bcv] <<  std::endl;
              amrex::Print() << "  Invalid direction: " << dir << std::endl;
              amrex::Abort("Fix the inputs file!");
            }
          }

          sum_same_loc += same_loc;
        }

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( sum_same_loc == 1, "Invalid bc region!");

      }


      // Store the BC ID for quick look-up when setting BC types
      if( dir_int == 0) {
        bc_xlo.push_back(bcv);
      } else if( dir_int == 1) {
        bc_xhi.push_back(bcv);
      } else if( dir_int == 2){
        bc_ylo.push_back(bcv);
      } else if( dir_int == 3){
        bc_yhi.push_back(bcv);
      } else if( dir_int == 4){
        bc_zlo.push_back(bcv);
      } else if( dir_int == 5){
        bc_zhi.push_back(bcv);
      }


      // Enforce the boundary for pressure outflows if specified.
      if( !po_noParOut && (new_bc.type == pinf_ || new_bc.type == pout_) ){
        domain_bc[dir_int ] = 0;
      }


      // flow_planes are used to create 'walls' in the level-set so that
      // particles don't fall out inflows (or outflows if desired).
      if( new_bc.type == minf_) {
        BC::flow_planes.emplace_back(point, normal, false);

      } else if( new_bc.type == pinf_ ) {
        BC::flow_planes.emplace_back(point, normal, false);

      } else if( new_bc.type == pout_ && !po_noParOut) {
        BC::flow_planes.emplace_back(point, normal, false);

      } else if (new_bc.type == nsw_) {
        BC::wall_planes.emplace_back(point, normal, false);

      }


      // Get fluid data.
      if(FLUID::solve) {

        std::string field = "bc."+regions[bcv]+"."+FLUID::name;
        amrex::ParmParse ppFluid(field.c_str());

        // Mass inflows need fluid velocity and volume fraction.
        if( new_bc.type == minf_) {
          ppFluid.get("volfrac", new_bc.fluid.volfrac);
          volfrac_total += new_bc.fluid.volfrac;
          ppFluid.getarr("velocity", new_bc.fluid.velocity, 0, 3);
        }

        // Read in fluid pressure
        new_bc.fluid.pressure_defined =
          ppFluid.query("pressure", new_bc.fluid.pressure);
        if(( new_bc.type == pinf_ || new_bc.type == pout_) &&
           !new_bc.fluid.pressure_defined) {
          amrex::Print() << "Pressure BCs must have pressure defined!" << std::endl;
          amrex::Abort("Fix the inputs file!");
        }

        ppFluid.query("temperature", new_bc.fluid.temperature);

      }

      if(DEM::solve) {

        // Get the list of solids used in defining the BC region
        std::vector<std::string> solids_types;
        {
          ppRegion.queryarr("solids", solids_types);
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

#if 0
    //Dump out what we read for debugging!
    for (int bcv(0); bcv<bc.size(); bcv++){


      amrex::Print() << std::endl << "Summarizing BC regions:\n" << std::endl;
      amrex::Print() << " BC: " << bcv << "    Type: " << bc[bcv].type << "\n\n" <<
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
#endif
  }// END Initialize

}
