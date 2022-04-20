#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_Geometry.H>

#include <AMReX_ParmParse.H>

#include <mfix_bc_parms.H>
#include <mfix_eb_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_regions_parms.H>
#include <mfix_species_parms.H>
#include <mfix_bc_list.H>

namespace BC
{

  // Direction of pressure drop (0:x, 1:y, 2:z)
  int delp_dir = -1;
  amrex::Real delp[3];

  int domain_bc[6];

  // EB planes for level-set creation
  amrex::Vector<amrex::EB2::PlaneIF> flow_planes;
  amrex::Vector<amrex::EB2::PlaneIF> wall_planes;

  std::bitset<6> flow_plane(std::string("000000"));

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

  std::array<amrex::LinOpBCType,3> diff_temperature_lobc;
  std::array<amrex::LinOpBCType,3> diff_temperature_hibc;

  std::array<amrex::LinOpBCType,3> diff_species_lobc;
  std::array<amrex::LinOpBCType,3> diff_species_hibc;


  // Data structure storing individual BC information
  amrex::Vector<BC_t> bc;


  void Initialize (amrex::Geometry& geom,
                   const FluidPhase& fluid,
                   const SolidsPhase& solids)
  {
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
    for (int dir=0; dir < 3; dir ++ ) {
      if (geom.isPeriodic(dir)) {

        ppe_lobc[dir] = amrex::LinOpBCType::Periodic;
        ppe_hibc[dir] = amrex::LinOpBCType::Periodic;

        diff_vel_lobc[dir] = amrex::LinOpBCType::Periodic;
        diff_vel_hibc[dir] = amrex::LinOpBCType::Periodic;

        diff_scal_lobc[dir] = amrex::LinOpBCType::Periodic;
        diff_scal_hibc[dir] = amrex::LinOpBCType::Periodic;

        diff_temperature_lobc[dir] = amrex::LinOpBCType::Periodic;
        diff_temperature_hibc[dir] = amrex::LinOpBCType::Periodic;

        diff_species_lobc[dir] = amrex::LinOpBCType::Periodic;
        diff_species_hibc[dir] = amrex::LinOpBCType::Periodic;

      } else {

        ppe_lobc[dir] = amrex::LinOpBCType::Neumann;
        ppe_hibc[dir] = amrex::LinOpBCType::Neumann;

        diff_vel_lobc[dir] = amrex::LinOpBCType::Dirichlet;
        diff_vel_hibc[dir] = amrex::LinOpBCType::Dirichlet;

        diff_scal_lobc[dir] = amrex::LinOpBCType::Dirichlet;
        diff_scal_hibc[dir] = amrex::LinOpBCType::Dirichlet;

        diff_temperature_lobc[dir] = amrex::LinOpBCType::Dirichlet;
        diff_temperature_hibc[dir] = amrex::LinOpBCType::Dirichlet;

        diff_species_lobc[dir] = amrex::LinOpBCType::Dirichlet;
        diff_species_hibc[dir] = amrex::LinOpBCType::Dirichlet;

      }
    }


    // Initialize the periodic pressure drop to zero in all directions
    for (int dir=0; dir < 3; dir ++)
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
    for(size_t bcv=0; bcv < regions.size(); bcv++){

      amrex::Real volfrac_total(0.0);

      BC_t new_bc(fluid);

      // Set the region for the initial condition.
      new_bc.region = REGIONS::getRegion(regions[bcv]);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE( new_bc.region != NULL, "Invalid bc region!");


      // Get the BC type (MI/PI/PO/NSW/EB)
      std::string bc_type;
      pp.get(regions[bcv].c_str(), bc_type);

      // Convert the input string into the integers
      if( bc_type == "mi") {
        new_bc.type = BCList::minf;
      } else if( bc_type == "pi") {
        new_bc.type = BCList::pinf;
      } else if( bc_type == "po") {
        new_bc.type = BCList::pout;
      } else if (bc_type == "nsw") {
        new_bc.type = BCList::nsw;
      } else if (bc_type == "eb") {
        new_bc.type = BCList::eb;
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

      if ( new_bc.type == BCList::nsw ) {

        // Walls need a normal because they might not be on a domain extent.
        amrex::Vector<amrex::Real> normal_in(3);
        ppRegion.getarr("normal", normal_in, 0, 3);

        normal={normal_in[0], normal_in[1], normal_in[2]};

      } else if ( new_bc.type != BCList::eb ) {

        // This covers mass inflows (mi), pressure inflows (pi), and
        // pressure outflows (po). These all must align with a domain
        // extent.

        int sum_same_loc(0);
        for (int dir(0); dir<3; ++dir){

          int same_loc = amrex::Math::abs(new_bc.region->lo(dir) - new_bc.region->hi(dir)) < tolerance ? 1 : 0;

          if (same_loc ){

            if ( amrex::Math::abs(new_bc.region->lo(dir) - plo[dir] ) < tolerance ){

              point[dir] = plo[dir]+1.0e-15;
              normal[dir] =  1.0;

              dir_int = 2*dir;

              if( new_bc.type == BCList::pinf || new_bc.type == BCList::pout) {
                ppe_lobc[dir] = amrex::LinOpBCType::Dirichlet;
                diff_vel_lobc[dir] = amrex::LinOpBCType::Neumann;
                diff_scal_lobc[dir] = amrex::LinOpBCType::Neumann;
              }

              if( new_bc.type == BCList::pout) {
                diff_temperature_lobc[dir] = amrex::LinOpBCType::Neumann;
                diff_species_lobc[dir] = amrex::LinOpBCType::Neumann;
              }

            } else if ( amrex::Math::abs( new_bc.region->hi(dir) - phi[dir] ) < tolerance ){

              point[dir] = phi[dir]-1.0e-15;
              normal[dir] = -1.0;

              dir_int = 2*dir+1;

              if( new_bc.type == BCList::pinf || new_bc.type == BCList::pout) {
                ppe_hibc[dir] = amrex::LinOpBCType::Dirichlet;
                diff_vel_hibc[dir] = amrex::LinOpBCType::Neumann;
                diff_scal_hibc[dir] = amrex::LinOpBCType::Neumann;
              }

              if( new_bc.type == BCList::pout) {
                diff_temperature_hibc[dir] = amrex::LinOpBCType::Neumann;
                diff_species_hibc[dir] = amrex::LinOpBCType::Neumann;
              }

            } else {
              amrex::Print() << "Flow BCs must be located on domain extents!" << std::endl;
              amrex::Print() << "BC Name: " << regions[bcv] <<  std::endl;
              amrex::Print() << "  Invalid direction: " << dir << std::endl;
              amrex::Abort("Fix the inputs file!");
            }

            // Flag that the level set should see these domain extents
            // as walls for particle collisions.
            if( new_bc.type == BCList::minf || new_bc.type == BCList::pinf ||
                (new_bc.type == BCList::pinf && po_noParOut == 1 )) {
              amrex::Print() << "X> FLAGGING FLOW PLANE " << dir_int << std::endl;
              flow_plane.flip(dir_int);
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
      if( new_bc.type == BCList::pout && po_noParOut == 0 ){
        domain_bc[dir_int ] = 0;
      }


      // flow_planes are used to create 'walls' in the level-set so that
      // particles don't fall out inflows (or outflows if desired).
      if (new_bc.type == BCList::minf) {
        BC::flow_planes.emplace_back(point, normal, false);

      } else if (new_bc.type == BCList::pinf) {
        BC::flow_planes.emplace_back(point, normal, false);

      } else if (new_bc.type == BCList::pout && po_noParOut == 1) {
        BC::flow_planes.emplace_back(point, normal, false);

      } else if (new_bc.type == BCList::nsw) {
        BC::wall_planes.emplace_back(point, normal, false);

      }


      // Get EB data.
      if (new_bc.type == BCList::eb) {
        std::string field = "bc."+regions[bcv]+".eb";
        amrex::ParmParse ppEB(field.c_str());

        if (ppEB.contains("temperature")) {
          EB::fix_temperature = 1;
          ppEB.get("temperature", new_bc.eb.temperature);
        }

        if (ppEB.contains("normal")) {
          amrex::Print() << "------------------------------------------------------------ FOUND NORMAL\n";
          new_bc.eb.has_normal = 1;

          ppEB.getarr("normal", new_bc.eb.normal, 0, 3);
          const amrex::Real nx = new_bc.eb.normal[0];
          const amrex::Real ny = new_bc.eb.normal[1];
          const amrex::Real nz = new_bc.eb.normal[2];

          const amrex::Real nmag = std::sqrt(nx*nx+ny*ny+nz*nz);
          AMREX_ALWAYS_ASSERT(nmag > 0.);

          new_bc.eb.normal[0] = nx / nmag;
          new_bc.eb.normal[1] = ny / nmag;
          new_bc.eb.normal[2] = nz / nmag;

          amrex::Real tol_deg(0.);
          ppEB.query("normal_tol", tol_deg);

          new_bc.eb.normal_tol = tol_deg*M_PI/amrex::Real(180.);

        }

      }

      // Get fluid data.
      if(fluid.solve) {

        std::string field = "bc."+regions[bcv]+"."+fluid.name;
        amrex::ParmParse ppFluid(field.c_str());

        // Mass inflows need fluid velocity and volume fraction.
        if( new_bc.type == BCList::minf ) {
          ppFluid.get("volfrac", new_bc.fluid.volfrac);
          volfrac_total += new_bc.fluid.volfrac;

          read_bc_velocity(ppFluid, &new_bc.fluid);

        } else if ( new_bc.type == BCList::eb ) {

          if (ppFluid.contains("velocity")) {
             new_bc.fluid.flow_thru_eb    = 1;
             new_bc.fluid.eb_has_velocity = 1;
             read_bc_velocity(ppFluid, &new_bc.fluid);

          } else if (ppFluid.contains("volflow")) {
             new_bc.fluid.flow_thru_eb   = 1;
             new_bc.fluid.eb_has_volflow = 1;
             read_bc_volflow(ppFluid, &new_bc.fluid);
          }

          if ( new_bc.fluid.flow_thru_eb ) {
            ppFluid.get("volfrac", new_bc.fluid.volfrac);
            volfrac_total += new_bc.fluid.volfrac;

            EB::has_flow = 1;
            EB::compute_area = EB::compute_area || new_bc.fluid.eb_has_volflow;
          }

        }

        new_bc.fluid.pressure = -1.0;

        // Read in fluid pressure
        new_bc.fluid.pressure_defined =
          ppFluid.query("pressure", new_bc.fluid.pressure);
        if(( new_bc.type == BCList::pinf || new_bc.type == BCList::pout) && (!new_bc.fluid.pressure_defined)) {
          amrex::Print() << "Pressure BCs must have pressure defined!" << std::endl;
          amrex::Print() << "BC region: " << regions[bcv] << std::endl;
          amrex::Abort("Fix the inputs file!");
        }

        // Read in fluid density
        read_bc_density(ppFluid, &new_bc.fluid);

        if( new_bc.type == BCList::minf || new_bc.type == BCList::pinf ||
              (new_bc.type == BCList::eb && new_bc.fluid.flow_thru_eb) ) {
          if (fluid.solve_enthalpy) {
            read_bc_temperature(ppFluid, &new_bc.fluid);
          }

          // Get species data.
          if (fluid.solve_species) {
            read_bc_species(ppFluid, fluid.species, &new_bc.fluid);
          }
        }

        if (fluid.constraint_type == ConstraintType::IncompressibleFluid) {

          if ((new_bc.type == BCList::minf || new_bc.type == BCList::pinf) && (!new_bc.fluid.density_defined)) {
            amrex::Print() << "Density BCs must have density defined!" << std::endl;
            amrex::Print() << "BC region: " << regions[bcv] << std::endl;
            amrex::Abort("Fix the inputs file!");
          }

        } else if (fluid.constraint_type == ConstraintType::IdealGasOpenSystem ||
                   fluid.constraint_type == ConstraintType::IdealGasClosedSystem) {

          // Read in fluid pressure
          read_bc_thermodynamic_pressure(ppFluid, &new_bc.fluid);

          if((new_bc.type == BCList::pinf || new_bc.type == BCList::pout) &&
             (!new_bc.fluid.thermodynamic_pressure_defined) &&
             fluid.constraint_type == ConstraintType::IdealGasOpenSystem) {
            amrex::Print() << "Pressure BCs must have pressure defined!" << std::endl;
            amrex::Print() << "BC region: " << regions[bcv] << std::endl;
            amrex::Abort("Fix the inputs file!");
          }

          if(new_bc.type == BCList::minf || new_bc.type == BCList::pinf ||
             (new_bc.type == BCList::eb && new_bc.fluid.flow_thru_eb)) {

            const int sum_defined = int(new_bc.fluid.density_defined) +
                                    int(new_bc.fluid.temperature_defined) +
                                    int(new_bc.fluid.thermodynamic_pressure_defined);

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sum_defined == 2,
                "Boundary conditions inputs must provide exactly two quantities"
                " among fluid temperature, density, and thermodynamic pressure");
          }

        }
      }

      // Get solid data.
      if(DEM::solve || PIC::solve) {

        // Get the list of solids used in defining the BC region
        std::vector<std::string> solids_types;
        {
          ppRegion.queryarr("solids", solids_types);
        }

        for(size_t lcs(0); lcs < solids_types.size(); ++ lcs){

          SolidsPhase::SOLIDS_t new_solid;

          std::string field = "bc."+regions[bcv]+"."+solids_types[lcs];
          amrex::ParmParse ppSolid(field.c_str());

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

          if (solids.solve_species) {

            const int nspecies_s = solids.nspecies;
            new_solid.species.resize(nspecies_s);

            std::string species_field = field+".species";
            amrex::ParmParse ppSpecies(species_field.c_str());

            amrex::Real total_mass_fraction(0);

            for (int n(0); n < solids.nspecies; n++) {
              // Get the name of the solid species we want to get the BC
              std::string dem_specie = solids.species[n];
              // Get the BC mass fraction for the current species
              ppSpecies.query(dem_specie.c_str(), new_solid.species[n].mass_fraction);

              total_mass_fraction += new_solid.species[n].mass_fraction;
            }

            // Sanity check that the input species mass fractions sum up to 1
            if (!(amrex::Math::abs(total_mass_fraction-1) < 1.e-15)) {
              std::string message = "Error: SOLID type " + solids_types[lcs]
                + " species BCs mass fractions in region " + regions[bcv]
                + " sum up to " + std::to_string(total_mass_fraction) + "\n";

              amrex::Abort(message);
            }
          }

          new_bc.solids.push_back(new_solid);
        }
      }

      BC::bc.push_back(new_bc);

    }

//    if (fluid.constraint_type == ConstraintType::IdealGasOpenSystem) {
//
//      for (int bcv(1); bcv < BC::bc.size(); ++bcv) {
//
//        auto bc_type = BC::bc[bcv].type;
//
//        if (bc_type == minf_ || bc_type == pinf_ || bc_type == pout_) {
//          const amrex::Real diff = std::abs(BC::bc[0].fluid.pressure -
//                                            BC::bc[bcv].fluid.pressure);
//
//          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(diff < 1.e-15,
//              "BCs for thermodynamic pressure are not uniform in space");
//        }
//      }
//    }

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

      if(fluid.solve){

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


      if(DEM::solve || PIC::solve){

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


void
read_bc_velocity (amrex::ParmParse pp, FluidPhase::FLUID_t *fluid)
{
  amrex::Vector<amrex::Real> vel_in;
  pp.getarr("velocity", vel_in);


  // Set flag indicating this is an EB boundary condition for fluid
  const int is_eb = fluid->flow_thru_eb;

  const int rcomps = vel_in.size();

  fluid->constant_velocity = ((rcomps == 3) || (is_eb && rcomps == 1) );
  if (fluid->constant_velocity) {
    fluid->velocity= vel_in;

  } else {

    int found;
    int k=0;
    do {
      amrex::Vector<amrex::Real> kth_input;
      found = pp.queryktharr("velocity", k, kth_input);
      if (found){
        // Assert that it has the correct length.
        if ( ( kth_input.size() != 4 ) || ( is_eb && kth_input.size() != 2 )) {
          amrex::Print() << "Bad velocity inputs specification:\n";
          for( int lc=0; lc<kth_input.size(); lc++)
            amrex::Print()  << "  " << kth_input[lc];
          amrex::Print() << std::endl;
          amrex::Abort("Fix input deck.");
        }
        const amrex::Real new_time = kth_input[0];
        const int len = fluid->vel_table.size();
        if (len == 0){
          fluid->vel_table.push_back(kth_input);
        } else {
          for (int lc=0; lc<len; lc++) {
            if ( new_time <= fluid->vel_table[lc][0]){
              fluid->vel_table.insert(fluid->vel_table.begin()+lc, kth_input);
            } else if (lc == len-1){
              fluid->vel_table.push_back(kth_input);
            }
          }
        }
      }
      k++;
    } while(found);
  }

  const int tsize = fluid->vel_table.size();
  if ( fluid->velocity.size() == 0 && tsize == 0){
    amrex::Print() << "No velocity inputs found for fluid\n";
    amrex::Abort("No velocity inputs found for fluid");
  }

  // Verify that all entries in the same are the same size
  if( tsize > 0) {
    const int ncomps = fluid->vel_table[0].size();
    for (int lc=1; lc<tsize; lc++) {
      if(fluid->vel_table[lc].size() != ncomps) {
        amrex::Print() << "Velocity inputs must be the same dimensions\n";
        amrex::Abort("Velocity inputs must be the same dimensions\n");
      }
    }
    fluid->eb_vel_is_mag = (is_eb && ncomps == 2);
  } else {
    fluid->eb_vel_is_mag = (is_eb && (fluid->velocity.size() == 1));
  }

}// end read_velocity


void
read_bc_volflow (amrex::ParmParse pp, FluidPhase::FLUID_t *fluid)
{

  amrex::Vector<amrex::Real> volflow_in;
  pp.getarr("volflow", volflow_in);

  fluid->constant_volflow = (volflow_in.size() == 1);
  if (fluid->constant_volflow) {

    fluid->volflow = volflow_in[0];

  } else {

    int found;
    int k=0;
    do {

      amrex::Vector<amrex::Real> kth_input;
      found = pp.queryktharr("volflow", k, kth_input);

      if (found){
        // Assert that it has the correct length.
        if ( kth_input.size() != 2 ) {
          amrex::Print() << "Bad volflow inputs specification:\n";
          for( int lc=0; lc<kth_input.size(); lc++)
            amrex::Print()  << "  " << kth_input[lc];
          amrex::Print() << std::endl;
          amrex::Abort("Fix input deck.");
        }

        const amrex::Real new_time = kth_input[0];
        const int len = fluid->volflow_table.size();
        if (len == 0){
          fluid->volflow_table.push_back(kth_input);
        } else {
          for (int lc=0; lc<len; lc++) {
            if ( new_time <= fluid->volflow_table[lc][0]){
              fluid->volflow_table.insert(fluid->volflow_table.begin()+lc, kth_input);
            } else if (lc == len-1){
              fluid->volflow_table.push_back(kth_input);
            }
          }
        }
      }
      k++;
    } while(found);

  }

}// end read_bc_volflow


void
read_bc_density (amrex::ParmParse pp, FluidPhase::FLUID_t *fluid)
{
  amrex::Vector<amrex::Real> rog_in;
  pp.queryarr("density", rog_in);

  if (rog_in.size() == 0) {
    fluid->density_defined = 0;
    return;
  }

  fluid->density_defined = 1;

  fluid->constant_density = (rog_in.size() == 1);
  if (fluid->constant_density) {

    fluid->density = rog_in[0];

  } else {

    int found;
    int k=0;
    do {

      amrex::Vector<amrex::Real> kth_input;
      found = pp.queryktharr("density", k, kth_input);

      if (found){
        // Assert that it has the correct length.
        if ( kth_input.size() != 2 ) {
          amrex::Print() << "Bad density inputs specification:\n";
          for( int lc=0; lc<kth_input.size(); lc++)
            amrex::Print()  << "  " << kth_input[lc];
          amrex::Print() << std::endl;
          amrex::Abort("Fix input deck.");
        }

        const amrex::Real new_time = kth_input[0];
        const int len = fluid->density_table.size();
        if (len == 0){
          fluid->density_table.push_back(kth_input);
        } else {
          for (int lc=0; lc<len; lc++) {
            if ( new_time <= fluid->density_table[lc][0]){
              fluid->density_table.insert(fluid->density_table.begin()+lc, kth_input);
            } else if (lc == len-1){
              fluid->density_table.push_back(kth_input);
            }
          }
        }
      }
      k++;
    } while(found);

  }

}// end read_bc_density


void
read_bc_thermodynamic_pressure (amrex::ParmParse pp, FluidPhase::FLUID_t *fluid)
{
  amrex::Vector<amrex::Real> pressure_in;
  pp.queryarr("thermodynamic_pressure", pressure_in);

  if (pressure_in.size() == 0) {
    fluid->thermodynamic_pressure_defined = 0;
    return;
  }

  fluid->thermodynamic_pressure_defined = 1;

  fluid->constant_thermodynamic_pressure = (pressure_in.size() == 1);
  if (fluid->constant_thermodynamic_pressure) {

    fluid->thermodynamic_pressure = pressure_in[0];

  } else {

    int found;
    int k=0;
    do {

      amrex::Vector<amrex::Real> kth_input;
      found = pp.queryktharr("thermodynamic_pressure", k, kth_input);

      if (found){
        // Assert that it has the correct length.
        if ( kth_input.size() != 2 ) {
          amrex::Print() << "Bad thermodynamic pressure inputs specification:\n";
          for( int lc=0; lc<kth_input.size(); lc++)
            amrex::Print()  << "  " << kth_input[lc];
          amrex::Print() << std::endl;
          amrex::Abort("Fix input deck.");
        }

        const amrex::Real new_time = kth_input[0];
        const int len = fluid->thermodynamic_pressure_table.size();
        if (len == 0){
          fluid->thermodynamic_pressure_table.push_back(kth_input);
        } else {
          for (int lc=0; lc<len; lc++) {
            if ( new_time <= fluid->thermodynamic_pressure_table[lc][0]){
              fluid->thermodynamic_pressure_table.insert(fluid->thermodynamic_pressure_table.begin()+lc, kth_input);
            } else if (lc == len-1){
              fluid->thermodynamic_pressure_table.push_back(kth_input);
            }
          }
        }
      }

      k++;
    } while(found);
  }
} // end read_bc_thermodynamic_pressure




void
read_bc_temperature (amrex::ParmParse pp, FluidPhase::FLUID_t *fluid)
{

  amrex::Vector<amrex::Real> tg_in;
  pp.queryarr("temperature", tg_in);

  if (tg_in.size() == 0) {
    fluid->temperature_defined = 0;
    return;
  }

  fluid->temperature_defined = 1;

  fluid->constant_temperature = (tg_in.size() == 1);
  if (fluid->constant_temperature) {

    fluid->temperature = tg_in[0];

  } else {

    int found;
    int k=0;
    do {

      amrex::Vector<amrex::Real> kth_input;
      found = pp.queryktharr("temperature", k, kth_input);

      if (found){
        // Assert that it has the correct length.
        if ( kth_input.size() != 2 ) {
          amrex::Print() << "Bad temperature inputs specification:\n";
          for( int lc=0; lc<kth_input.size(); lc++)
            amrex::Print()  << "  " << kth_input[lc];
          amrex::Print() << std::endl;
          amrex::Abort("Fix input deck.");
        }

        const amrex::Real new_time = kth_input[0];
        const int len = fluid->tg_table.size();
        if (len == 0){
          fluid->tg_table.push_back(kth_input);
        } else {
          for (int lc=0; lc<len; lc++) {
            if ( new_time <= fluid->tg_table[lc][0]){
              fluid->tg_table.insert(fluid->tg_table.begin()+lc, kth_input);
            } else if (lc == len-1){
              fluid->tg_table.push_back(kth_input);
            }
          }
        }


      }

      k++;
    } while(found);

  }

}// end read_bc_temperature


void
read_bc_species (amrex::ParmParse pp,
                 const amrex::Vector<std::string> species,
                 FluidPhase::FLUID_t *fluid)
{

  const int nspecies = species.size();

  fluid->species.resize(nspecies);
  fluid->constant_species.resize(nspecies);
  fluid->species_table.resize(nspecies);

  for (int n(0); n < nspecies; n++) {

    amrex::Vector<amrex::Real> species_in;

    // Get the name of the fluid species we want to get the IC
    std::string species_field = "species."+species[n];

    pp.getarr(species_field.c_str(), species_in);

    int constant_species = (species_in.size() == 1 ? 1 : 0);
    fluid->constant_species[n] = constant_species;

    if (constant_species) {
      fluid->species[n].mass_fraction = species_in[0];
    } else {

      int found;
      int k=0;

      do {

        amrex::Vector<amrex::Real> kth_input;
        found = pp.queryktharr(species_field.c_str(), k, kth_input);

        if (found){
          // Assert that it has the correct length.
          if ( kth_input.size() != 2 ) {
            amrex::Print() << "Bad bc species inputs specification:\n";
            for( int lc=0; lc<kth_input.size(); lc++)
              amrex::Print()  << "  " << kth_input[lc];
            amrex::Print() << std::endl;
            amrex::Abort("Fix input deck.");
          }

          const amrex::Real new_time = kth_input[0];
          const int len =  fluid->species_table[n].size();
          if (len == 0){
            fluid->species_table[n].push_back(kth_input);
          } else {
            for (int lc=0; lc<len; lc++) {
              if ( new_time <= fluid->species_table[n][lc][0]) {
                fluid->species_table[n].insert(fluid->species_table[n].begin()+lc, kth_input);
              } else if (lc == len-1){
                fluid->species_table[n].push_back(kth_input);
              }
            }
          }
        }
        k++;
      } while(found);
    }
  }


  int sum_constant_species = 0;
  for (int n(0); n < nspecies; n++) {
    sum_constant_species += fluid->constant_species[n];
  }

  bool err_found = false;

  // Verify that mass fractions sum to one.
  if (sum_constant_species == nspecies) {

    amrex::Real sum_xg = 0.0;
    for (int n(0); n < nspecies; n++) {
      sum_xg += fluid->species[n].mass_fraction;
    }
    err_found = (amrex::Math::abs(1.0 - sum_xg) > 1.e-15);

  }else if (sum_constant_species == 0 ) {

    // First verify that all species are given for all times
    const int entries = fluid->species_table[0].size();

    // Verify that all species have the same number
    // of time entries.
    for (int n(0); n < nspecies; n++) {
      err_found = err_found || (fluid->species_table[n].size() != entries);
      // amrex::Print() << " n: " << n << "  entries: "
      // << fluid->species_table[n].size() << "   err: " << err_found << "\n";
    }

    // Verify that all time entries are the same and that the
    // mass fractions sum to one.
    if(!err_found){
      for (int entry(0); entry<entries; entry++) {

        const amrex::Real time = fluid->species_table[0][entry][0];
        amrex::Real sum_xg = 0.0;
        for (int n(0); n < nspecies; n++) {
          sum_xg += fluid->species_table[n][entry][1];
          err_found = err_found || (amrex::Math::abs(fluid->species_table[n][entry][0] - time) > 1.e-15);
        }
        err_found = (amrex::Math::abs(1.0 - sum_xg) > 1.e-15);
        // amrex::Print() << " time " << time << "   sum_xg: " << sum_xg << "   err: " << err_found << "\n";
      }
    }

  } else{
    err_found = true;
  }

  if(err_found){
    amrex::Print() << "\n\n"
                   << " Bad species bc inputs specification. Check the inputs file:\n"
                   << "  > all species mass fractions sum to one\n"
                   << "  > all flow bcs for species are either constant or time dependent\n"
                   << "  > time dependent bcs must all start/end at the same times\n";
    amrex::Abort("Fix input deck.");
  }

}// end read_bc_species




}
