#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_Geometry.H>

#include <AMReX_ParmParse.H>

#include <mfix_bc.H>
#include <mfix_eb.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_regions.H>
#include <mfix_species.H>

using namespace amrex;


// Constructor
BCList::BCList (const int nlev_in)
  : nlev(nlev_in)
{}


// Destructor
BCList::~BCList()
{
  for (int lev(0); lev < nlev; ++lev) {
    delete bc_ilo[lev];
    delete bc_ihi[lev];
    delete bc_jlo[lev];
    delete bc_jhi[lev];
    delete bc_klo[lev];
    delete bc_khi[lev];
  }
}


void
BCList::MakeBCArrays (int nghost,
                      bool ooo_debug,
                      amrex::Vector<amrex::Geometry>& geom)
{
  for (int lev = 0; lev < bc_ilo.size(); lev++)
  {
    if (bc_ilo[lev] != nullptr) delete bc_ilo[lev];
    if (bc_ihi[lev] != nullptr) delete bc_ihi[lev];
    if (bc_jlo[lev] != nullptr) delete bc_jlo[lev];
    if (bc_jhi[lev] != nullptr) delete bc_jhi[lev];
    if (bc_klo[lev] != nullptr) delete bc_klo[lev];
    if (bc_khi[lev] != nullptr) delete bc_khi[lev];
  }

  if (ooo_debug) amrex::Print() << "MakeBCArrays" << std::endl;

  bc_ilo.clear(); bc_ilo.resize(nlev, nullptr);
  bc_ihi.clear(); bc_ihi.resize(nlev, nullptr);
  bc_jlo.clear(); bc_jlo.resize(nlev, nullptr);
  bc_jhi.clear(); bc_jhi.resize(nlev, nullptr);
  bc_klo.clear(); bc_klo.resize(nlev, nullptr);
  bc_khi.clear(); bc_khi.resize(nlev, nullptr);

  for (int lev = 0; lev < nlev; lev++)
  {
    // Define and allocate the integer MultiFab that is the outside adjacent
    // cells of the problem domain.
    Box domainx(geom[lev].Domain());
    domainx.grow(1,nghost);
    domainx.grow(2,nghost);
    Box box_ilo = amrex::adjCellLo(domainx,0,1);
    Box box_ihi = amrex::adjCellHi(domainx,0,1);

    Box domainy(geom[lev].Domain());
    domainy.grow(0,nghost);
    domainy.grow(2,nghost);
    Box box_jlo = amrex::adjCellLo(domainy,1,1);
    Box box_jhi = amrex::adjCellHi(domainy,1,1);

    Box domainz(geom[lev].Domain());
    domainz.grow(0,nghost);
    domainz.grow(1,nghost);
    Box box_klo = amrex::adjCellLo(domainz,2,1);
    Box box_khi = amrex::adjCellHi(domainz,2,1);

    // Note that each of these is a single IArrayBox so every process has a copy of them
    bc_ilo[lev] = new IArrayBox(box_ilo,2);
    bc_ihi[lev] = new IArrayBox(box_ihi,2);
    bc_jlo[lev] = new IArrayBox(box_jlo,2);
    bc_jhi[lev] = new IArrayBox(box_jhi,2);
    bc_klo[lev] = new IArrayBox(box_klo,2);
    bc_khi[lev] = new IArrayBox(box_khi,2);
  }
}


MFIXBoundaryConditions::MFIXBoundaryConditions (amrex::Vector<amrex::Geometry>& geom,
                                                BCList& bc_list,
                                                MFIXEmbeddedBoundaries& embedded_boundaries)
  : m_geom(geom)
  , m_bc_list(bc_list)
  , m_bc_u_g(50, 0)
  , m_bc_v_g(50, 0)
  , m_bc_w_g(50, 0)
  , m_bc_t_g(50, 0)
  , m_bc_h_g(50, 0)
  , m_bc_ro_g(50, 0)
  , m_bc_tracer(50, 0)
  , m_bc_ep_g(50, 0)
  , m_bc_p_g(50, 0)
  , m_delp_dir(-1)
  , m_flow_plane(std::string("000000"))
  , m_embedded_boundaries(embedded_boundaries)
{
  m_vel_g_bc_types["Dirichlet"] = {BCList::minf};
  m_vel_g_bc_types["Neumann"] = {BCList::pinf, BCList::pout};

  m_ro_g_bc_types["Dirichlet"] = {BCList::minf};
  m_ro_g_bc_types["Neumann"] = {BCList::pinf, BCList::pout};

  m_T_g_bc_types["Dirichlet"] = {BCList::minf, BCList::pinf};
  m_T_g_bc_types["Neumann"] = {BCList::pout};

  m_trac_g_bc_types["Dirichlet"] = {BCList::minf};
  m_trac_g_bc_types["Neumann"] = {BCList::pinf, BCList::pout};

  m_X_gk_bc_types["Dirichlet"] = {BCList::minf, BCList::pinf};
  m_X_gk_bc_types["Neumann"] = {BCList::pout};
}


void
MFIXBoundaryConditions::Initialize (amrex::Geometry& geom,
                                    const MFIXRegions& regions,
                                    MFIXFluidPhase& fluid,
                                    const MFIXSolidsPhase& solids,
                                    MFIXDEM& dem,
                                    MFIXPIC& pic)
{
  // Set flag to keep particles from leaving unless periodic.
  for (int dir(0); dir<3; ++dir) {
    if (geom.isPeriodic(dir)) {
      m_domain_bc[2*dir  ] = 0;
      m_domain_bc[2*dir+1] = 0;
    } else {
      m_domain_bc[2*dir  ] = 1;
      m_domain_bc[2*dir+1] = 1;
    }
  }

  // Default all sides of the domain to Neumann
  for (int dir=0; dir < 3; ++dir) {
    if (geom.isPeriodic(dir)) {

      m_ppe_lobc[dir] = amrex::LinOpBCType::Periodic;
      m_ppe_hibc[dir] = amrex::LinOpBCType::Periodic;

      m_diff_vel_lobc[dir] = amrex::LinOpBCType::Periodic;
      m_diff_vel_hibc[dir] = amrex::LinOpBCType::Periodic;

      m_diff_scal_lobc[dir] = amrex::LinOpBCType::Periodic;
      m_diff_scal_hibc[dir] = amrex::LinOpBCType::Periodic;

      m_diff_temperature_lobc[dir] = amrex::LinOpBCType::Periodic;
      m_diff_temperature_hibc[dir] = amrex::LinOpBCType::Periodic;

      m_diff_species_lobc[dir] = amrex::LinOpBCType::Periodic;
      m_diff_species_hibc[dir] = amrex::LinOpBCType::Periodic;

    } else {

      m_ppe_lobc[dir] = amrex::LinOpBCType::Neumann;
      m_ppe_hibc[dir] = amrex::LinOpBCType::Neumann;

      m_diff_vel_lobc[dir] = amrex::LinOpBCType::Dirichlet;
      m_diff_vel_hibc[dir] = amrex::LinOpBCType::Dirichlet;

      m_diff_scal_lobc[dir] = amrex::LinOpBCType::Dirichlet;
      m_diff_scal_hibc[dir] = amrex::LinOpBCType::Dirichlet;

      m_diff_temperature_lobc[dir] = amrex::LinOpBCType::Dirichlet;
      m_diff_temperature_hibc[dir] = amrex::LinOpBCType::Dirichlet;

      m_diff_species_lobc[dir] = amrex::LinOpBCType::Dirichlet;
      m_diff_species_hibc[dir] = amrex::LinOpBCType::Dirichlet;

    }
  }

  // Initialize the periodic pressure drop to zero in all directions
  for (int dir=0; dir < 3; ++dir)
    m_delp[dir] = 0.0;

  amrex::ParmParse pp("bc");

  amrex::Real delp_in(0.0);
  pp.query("delp", delp_in);

  if (delp_in != 0.0) {

    pp.query("delp_dir", m_delp_dir);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0 <= m_delp_dir && m_delp_dir <= 2,
       "Direction of pressure drop (bc.delp_dir) must be specified.");

    m_delp[m_delp_dir] = delp_in;
  }

  constexpr amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

  // Flag to see if particles 'see' a wall at pressure outflows
  int po_noParOut = 0; // default behavior for PO's -- letting particles exit the domain
  pp.query("po_no_par_out", po_noParOut);

  // Get the list of region names used to define BCs
  std::vector<std::string> input_regions;
  pp.queryarr("regions", input_regions);

  // Loop over BCs
  for(size_t bcv=0; bcv < input_regions.size(); bcv++) {

    amrex::Real volfrac_total(0.0);

    BC_t new_bc(&fluid);

    // Set the region for the initial condition.
    new_bc.region = regions.getRegion(input_regions[bcv]);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(new_bc.region != nullptr, "Invalid bc region!");

    // Get the BC type (MI/PI/PO/NSW/EB)
    std::string bc_type;
    pp.get(input_regions[bcv].c_str(), bc_type);

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

    std::string field_bcRegion = "bc."+input_regions[bcv];
    amrex::ParmParse ppRegion(field_bcRegion.c_str());

    int face = -1;

    if (new_bc.type == BCList::nsw) {

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

        if (same_loc) {

          if (amrex::Math::abs(new_bc.region->lo(dir) - plo[dir] ) < tolerance) {

            point[dir] = plo[dir]+1.0e-15;
            normal[dir] =  1.0;

            face = 2*dir;

            if (new_bc.type == BCList::pinf || new_bc.type == BCList::pout) {
              m_ppe_lobc[dir] = amrex::LinOpBCType::Dirichlet;
              m_diff_vel_lobc[dir] = amrex::LinOpBCType::Neumann;
              m_diff_scal_lobc[dir] = amrex::LinOpBCType::Neumann;
            }

            if (new_bc.type == BCList::pout) {
              m_diff_temperature_lobc[dir] = amrex::LinOpBCType::Neumann;
              m_diff_species_lobc[dir] = amrex::LinOpBCType::Neumann;
            }

          } else if (amrex::Math::abs(new_bc.region->hi(dir) - phi[dir]) < tolerance) {

            point[dir] = phi[dir]-1.0e-15;
            normal[dir] = -1.0;

            face = 2*dir+1;

            if (new_bc.type == BCList::pinf || new_bc.type == BCList::pout) {
              m_ppe_hibc[dir] = amrex::LinOpBCType::Dirichlet;
              m_diff_vel_hibc[dir] = amrex::LinOpBCType::Neumann;
              m_diff_scal_hibc[dir] = amrex::LinOpBCType::Neumann;
            }

            if (new_bc.type == BCList::pout) {
              m_diff_temperature_hibc[dir] = amrex::LinOpBCType::Neumann;
              m_diff_species_hibc[dir] = amrex::LinOpBCType::Neumann;
            }

          } else {
            amrex::Print() << "Flow BCs must be located on domain extents!" << std::endl;
            amrex::Print() << "BC Name: " << input_regions[bcv] <<  std::endl;
            amrex::Print() << "  Invalid direction: " << dir << std::endl;
            amrex::Abort("Fix the inputs file!");
          }

          // Flag that the level set should see these domain extents
          // as walls for particle collisions.
          if (new_bc.type == BCList::minf || new_bc.type == BCList::pinf ||
              (new_bc.type == BCList::pinf && po_noParOut == 1)) {
            m_flow_plane.flip(face);
          }
        }

        sum_same_loc += same_loc;
      }

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE( sum_same_loc == 1, "Invalid bc region!");
    }

    // Store the BC ID for quick look-up when setting BC types
    if (face == 0) {
      m_bc_lo[0].push_back(bcv);
    } else if (face == 1) {
      m_bc_hi[0].push_back(bcv);
    } else if (face == 2) {
      m_bc_lo[1].push_back(bcv);
    } else if (face == 3) {
      m_bc_hi[1].push_back(bcv);
    } else if (face == 4) {
      m_bc_lo[2].push_back(bcv);
    } else if (face == 5) {
      m_bc_hi[2].push_back(bcv);
    }

    m_dir.push_back(face);

    // Enforce the boundary for pressure outflows if specified.
    if (new_bc.type == BCList::pout && po_noParOut == 0) {
      m_domain_bc[face] = 0;
    }

    // m_flow_planes are used to create 'walls' in the level-set so that
    // particles don't fall out inflows (or outflows if desired).
    if (new_bc.type == BCList::minf) {
      m_flow_planes.emplace_back(point, normal, false);

    } else if (new_bc.type == BCList::pinf) {
      m_flow_planes.emplace_back(point, normal, false);

    } else if (new_bc.type == BCList::pout && po_noParOut == 1) {
      m_flow_planes.emplace_back(point, normal, false);

    } else if (new_bc.type == BCList::nsw) {
      m_wall_planes.emplace_back(point, normal, false);

    }

    // Get EB data.
    if (new_bc.type == BCList::eb) {
      std::string field = "bc."+input_regions[bcv]+".eb";
      amrex::ParmParse ppEB(field.c_str());

      if (ppEB.contains("temperature")) {
        m_embedded_boundaries.set_fix_temperature(1);
        ppEB.get("temperature", new_bc.eb.temperature);
      }

      if (ppEB.contains("normal")) {
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
    if(fluid.solve()) {

      std::string field = "bc."+input_regions[bcv]+"."+fluid.names(0);
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

          m_embedded_boundaries.set_has_flow(1);
        }
      }

      // Read in density, temperature and species BC inputs only if BC type is
      // minf, pinf or flow_through_eb
      if (new_bc.type == BCList::minf || new_bc.type == BCList::pinf ||
          (new_bc.type == BCList::eb && new_bc.fluid.flow_thru_eb)) {

        // Read in fluid density
        read_bc_density(ppFluid, &new_bc.fluid);

        read_bc_temperature(ppFluid, &new_bc.fluid);

        // Get species data.
        if (fluid.solve_species())
          read_bc_species(ppFluid, fluid.species_names(), &new_bc.fluid);
      }

      // Read in fluid pressure BC inputs only if BC type is pinf or pout
      if (new_bc.type == BCList::pinf || new_bc.type == BCList::pout) {

        // Read in fluid pressure
        new_bc.fluid.pressure_defined = ppFluid.query("pressure", new_bc.fluid.pressure);

        if (!new_bc.fluid.pressure_defined) {
          amrex::Print() << "Pressure BCs must have pressure defined!" << std::endl;
          amrex::Print() << "BC region: " << input_regions[bcv] << std::endl;
          amrex::Abort("Fix the inputs file!");
        }
      }

      // ********************************************************************
      // Check BCs consistency depending on fluid constraint type
      // ********************************************************************
      if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IncompressibleFluid) {

        if (new_bc.type == BCList::minf || new_bc.type == BCList::pinf ||
            (new_bc.type == BCList::eb && new_bc.fluid.flow_thru_eb)) {

          if (!new_bc.fluid.density_defined) {
            amrex::Print() << "BCs must have density defined!" << std::endl;
            amrex::Print() << "BC region: " << input_regions[bcv] << std::endl;
            amrex::Abort("Fix the inputs file!");
          }
        }

      } else if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasOpenSystem) {

        if (new_bc.type == BCList::minf || new_bc.type == BCList::pinf ||
            (new_bc.type == BCList::eb && new_bc.fluid.flow_thru_eb)) {

          if (new_bc.fluid.density_defined) {
            amrex::Print() << "BCs must NOT have density defined!" << std::endl;
            amrex::Print() << "BC region: " << input_regions[bcv] << std::endl;
            amrex::Abort("Fix the inputs file!");
          }

          if (!new_bc.fluid.temperature_defined) {
            amrex::Print() << "BCs must have temperature defined!" << std::endl;
            amrex::Print() << "BC region: " << input_regions[bcv] << std::endl;
            amrex::Abort("Fix the inputs file!");
          }
        }

      } else if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem) {

        if (new_bc.type == BCList::minf || new_bc.type == BCList::pinf || new_bc.type == BCList::pout) {

          amrex::Print() << "minf, pinf or pout BCs cannot be defined in a closed system!" << std::endl;
          amrex::Print() << "BC region: " << input_regions[bcv] << std::endl;
          amrex::Abort("Fix the inputs file!");
        }

        if (new_bc.type == BCList::eb && new_bc.fluid.flow_thru_eb) {

          amrex::Print() << "Flux through EB cannot be defined in a closed system!" << std::endl;
          amrex::Print() << "BC region: " << input_regions[bcv] << std::endl;
          amrex::Abort("Fix the inputs file!");
        }
      }

    }

    // Get solid data.
    if(dem.solve() || pic.solve()) {

      // Get the list of solids used in defining the BC region
      std::vector<std::string> solids_names;
      {
        ppRegion.queryarr("solids", solids_names);
      }

      for(size_t lcs(0); lcs < solids_names.size(); ++ lcs){

        SOLIDS_t new_solid;

        new_solid.name = solids_names[lcs];
        new_solid.phase = solids.name_to_phase(solids_names[lcs]);

        std::string field = "bc."+input_regions[bcv]+"."+solids_names[lcs];
        amrex::ParmParse ppSolid(field.c_str());

        ppSolid.get("volfrac", new_solid.volfrac);
        volfrac_total += new_bc.fluid.volfrac;

        bool const is_mi(new_bc.type == BCList::minf);
        bool const is_eb(new_bc.type == BCList::eb);

        if (new_solid.volfrac > tolerance) {

          if (is_mi && dem.solve()) {

            amrex::Print() << "Mass inflows for DEM solids has not bee implemented yet!\n";
            amrex::Abort("Mass inflows for DEM solids has not bee implemented yet!");

          } else if (is_eb || (is_mi && pic.solve()) ) {

            int const dir_lohi(face%2); // lo=0, hi=1
            int const dir((face-dir_lohi)/((int)2)); // dir=0,1,2 (x,y,z)

            new_solid.statwt = 1.0;
            ppSolid.query("statwt", new_solid.statwt);

            if (ppSolid.contains("velocity")) {

              amrex::Vector<amrex::Real> vel_in;
              ppSolid.queryarr("velocity", vel_in);

              const int rcomps = vel_in.size();

              if (rcomps == 3) {
                new_solid.velocity = vel_in;
                new_solid.velmag = std::sqrt(vel_in[0]*vel_in[0]
                                           + vel_in[1]*vel_in[1]
                                           + vel_in[2]*vel_in[2]);

              } else if (rcomps == 1) {
                new_solid.velmag = amrex::Math::abs(vel_in[0]);

                new_solid.velocity.push_back( (dir == 0 ? vel_in[0] : 0.) );
                new_solid.velocity.push_back( (dir == 1 ? vel_in[0] : 0.) );
                new_solid.velocity.push_back( (dir == 2 ? vel_in[0] : 0.) );

              } else {
                std::string message = " Error: BC region " + input_regions[bcv]
                  + " has an invalid velocity for solids " +  solids_names[lcs] + "\n"
                  + " Solids velocity can be a scalar (velocity magnitude) or all"
                  + " three velocity components can be supplied.\n Solids inflow"
                  + " cannot not be a function of time.\n";
                amrex::Print() << message;
                amrex::Abort(message);
              }
            } else if (ppSolid.contains("volflow")) {
              ppSolid.get("volflow", new_solid.volflow);

              // +1 low side; -1 high side
              Real const r_normal_dir = (dir_lohi == 0) ? 1. : -1.;

              new_solid.velocity.push_back( ( dir == 0 ? r_normal_dir : 0.) );
              new_solid.velocity.push_back( ( dir == 1 ? r_normal_dir : 0.) );
              new_solid.velocity.push_back( ( dir == 2 ? r_normal_dir : 0.) );

            } else {
                std::string message = " Error: BC region " + input_regions[bcv]
                  + " does not have velocity or volflow specified for solids "
                  +  solids_names[lcs] + "\n";
                amrex::Print() << message;
                amrex::Abort(message);
            }

            // This probably needs a check if we are solving energy equations
            // and if so, require temperature be defined.
            ppSolid.query("temperature", new_solid.temperature);

            { // Get information about diameter distribution.
              int err = new_solid.diameter.set(field, "diameter");
              if(err) {
                std::string message = " Error: BC region " + input_regions[bcv]
                  + " has invalid diameter parameters for solids " + solids_names[lcs] + "\n";
                amrex::Print() << message;
                amrex::Abort(message);
              }
            }

            { // Get information about density distribution
              int err = new_solid.density.set(field, "density");
              if(err) {
                std::string message = " Error: BC region " + input_regions[bcv]
                  + " has invalid density parameters for solids " + solids_names[lcs] + "\n";
                amrex::Print() << message;
                amrex::Abort(message);
              }
            }

            if (solids.solve_species()) {

              const int nspecies_s = solids.nspecies();
              new_solid.species.resize(nspecies_s);

              std::string species_field = field+".species";
              amrex::ParmParse ppSpecies(species_field.c_str());

              amrex::Real total_mass_fraction(0);

              for (int n(0); n < solids.nspecies(); n++) {
                // Get the name of the solid species we want to get the BC
                std::string dem_specie = solids.species_names(n);
                // Get the BC mass fraction for the current species
                ppSpecies.query(dem_specie.c_str(), new_solid.species[n].mass_fraction);

                total_mass_fraction += new_solid.species[n].mass_fraction;
              }

              // Sanity check that the input species mass fractions sum up to 1
              if (!(amrex::Math::abs(total_mass_fraction-1) < 1.e-15)) {
                std::string message = "Error: SOLID phase " + solids_names[lcs]
                  + " species BCs mass fractions in region " + input_regions[bcv]
                  + " sum up to " + std::to_string(total_mass_fraction) + "\n";

                amrex::Abort(message);
              }
            } // end solve species
          } // End EB inflow
        } // volfrac > tolerance

        new_bc.solids.push_back(new_solid);
      }
    }

    m_bc.push_back(new_bc);
  }

  if (fluid.solve() &&
      (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IncompressibleFluid ||
       fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasOpenSystem)) {

    Vector<Real> bc_values(0);

    for (int bcv(0); bcv < m_bc.size(); ++bcv) {

      if (m_bc[bcv].type == BCList::pout)
        bc_values.push_back(m_bc[bcv].fluid.pressure);
    }

    if (bc_values.size() > 0) {
      const Real max_val = *std::max_element(bc_values.begin(), bc_values.end());
      const Real min_val = *std::min_element(bc_values.begin(), bc_values.end());

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(max_val-min_val) < 1.e-15,
        "Error: pressure values at pout boundaries differ");

      if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasOpenSystem) {

        if (fluid.p_therm_defined()) {
          const Real p_therm = fluid.thermodynamic_pressure();

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(p_therm-min_val) < 1.e-15 &&
              std::abs(p_therm-max_val) < 1.e-15 && p_therm > 1.e-15,
              "Error: either pout input pressures differ from fluid thermodynamic "
              "pressure, or given thermodynamic pressure is zero or negative");

        } else {

          fluid.set_thermodynamic_pressure(max_val);
          fluid.set_p_therm_defined(1);

          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.thermodynamic_pressure() > 1.e-15,
              "Error: fluid thermodynamic pressure is zero or negative");

        }
      }
    }
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

    if(fluid.solve()){

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


    if(m_dem.solve() || m_pic.solve()){

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

        amrex::Print() << " density distribution: " << bc[bcv].solids[lcs].density.distribution << std::endl;
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
} // END Initialize


void
MFIXBoundaryConditions::read_bc_velocity (amrex::ParmParse pp,
                                          FLUID_t *fluid)
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
        if ( ( kth_input.size() != 4 ) || ( !is_eb && kth_input.size() == 2 )) {
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

} // end read_velocity


void
MFIXBoundaryConditions::read_bc_volflow (amrex::ParmParse pp,
                                         FLUID_t *fluid)
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

} // end read_bc_volflow


void
MFIXBoundaryConditions::read_bc_density (amrex::ParmParse pp,
                                         FLUID_t *fluid)
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

} // end read_bc_density


void
MFIXBoundaryConditions::read_bc_temperature (amrex::ParmParse pp,
                                             FLUID_t *fluid)
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

} // end read_bc_temperature


void
MFIXBoundaryConditions::read_bc_species (amrex::ParmParse pp,
                                         const amrex::Vector<std::string> species,
                                         FLUID_t *fluid)
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

  } else if (sum_constant_species == 0 ) {

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

} // end read_bc_species



amrex::Real& MFIXBoundaryConditions::
get_bc_area (const int a_bcv)
{
  AMREX_ASSERT( a_bcv < m_bc.size() );
  AMREX_ASSERT( m_bc[a_bcv].type == BCList::pinf ||
                m_bc[a_bcv].type == BCList::pout ||
                m_bc[a_bcv].type == BCList::minf ||
                m_bc[a_bcv].type == BCList::eb );

  return m_area[a_bcv];
}

amrex::Real& MFIXBoundaryConditions::
get_bc_fab_area (const int a_bcv, std::pair<int,int> a_index)
{
  AMREX_ASSERT( a_bcv < m_bc.size() );
  AMREX_ASSERT( m_bc[a_bcv].type == BCList::pinf ||
                m_bc[a_bcv].type == BCList::pout ||
                m_bc[a_bcv].type == BCList::minf ||
                m_bc[a_bcv].type == BCList::eb );

  return m_fab_area[a_bcv][a_index];
}
