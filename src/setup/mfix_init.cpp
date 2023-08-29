#include <AMReX_ParmParse.H>
#include <AMReX_EBAmrUtil.H>

#include <mfix.H>
#include <mfix_init_fluid.H>
#include <mfix_regions.H>
#include <mfix_bc.H>
#include <mfix_ic.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_species.H>
#include <mfix_mlmg_options.H>
#include <mfix_utils.H>
#include <mfix_monitors.H>

using MFIXParIter = MFIXParticleContainer::MFIXParIter;
using PairIndex = MFIXParticleContainer::PairIndex;

void
mfix::InitParams ()
{
  if (ooo_debug)
    amrex::Print() << "InitParams" << std::endl;

  m_timer.Initialize();
  regions.Initialize();

  // Read and process species, fluid and DEM particle model options.
  species.Initialize();
  reactions.Initialize(species);
  fluid.Initialize(species, reactions);
  solids.Initialize(species, reactions);

  m_dem.Initialize();
  m_pic.Initialize();

  // Read in regions, initial and boundary conditions. Note that
  // regions need to be processed first as they define the
  // physical extents of ICs and BCs.
  m_boundary_conditions.Initialize(geom[0], regions, fluid, solids, m_dem, m_pic);
  m_initial_conditions.Initialize(regions, fluid, solids, m_dem, m_pic);

  m_rw->Initialize();

  // set n_error_buf (used in AmrMesh) to default (can overwrite later)
  for (int i = 0; i < n_error_buf.size(); i++)
    n_error_buf[i] = {8,8,8};

  {
    ParmParse pp("fluid.newton_solver");

    pp.query("absolute_tol", newton_abstol);
    pp.query("relative_tol", newton_reltol);
    pp.query("max_iterations", newton_maxiter);
  }

  {
    ParmParse pp("mfix");

    // Verbosity and MLMG parameters are now ParmParse with "nodal_proj" in the
    // inputs file
    // Examples: nodal_proj.verbose = 1
    //           nodal_proj.bottom_verbose = 1
    //           nodal_proj.maxiter
    //           nodal_proj.bottom_maxiter
    //           nodal_proj.bottom_rtol
    //           nodal_proj.bottom_atol
    //           nodal_proj.bottom_solver
    // More info at "AMReX-Hydro/Projections/hydro_NodalProjector.cpp"
    nodalproj_options = std::make_unique<MfixUtil::MLMGOptions>("nodal_proj");

    // Is this a steady-state calculation
    m_steady_state = 0;
    pp.query("steady_state", m_steady_state);

    // Tolerance to check for steady state
    steady_state_tol = -1.;
    pp.query("steady_state_tol", steady_state_tol);

    // Maximum number of iterations allowed to reach steady state
    pp.query("steady_state_maxiter", steady_state_maxiter);

    if (m_steady_state > 0) {
      if (steady_state_tol < 0)
        amrex::Abort("Must set steady_state_tol if running to steady state!");

      amrex::Print() << "Running to steady state with maxiter = "
                     << steady_state_maxiter << " and tolerance "
                     << steady_state_tol << std::endl;
    }
    else if (steady_state_tol > 0)
      amrex::Abort("steady_state_tol set but not steady_state!");

    // Flag to set verbosity
    m_verbose = 0;
    pp.query("verbose", m_verbose);

    pp.query("ooo_debug", ooo_debug);

    // Flag to invoke UDFs
    bool call_usr_bool = false;
    pp.query("call_usr", call_usr_bool);
    call_udf = call_usr_bool ? 1 : 0; // Set global flag

    Array<Real,3> gravity_in{0.0, 0.0, 0.0};
    pp.get("gravity", gravity_in);
    for (int dir = 0; dir < 3; dir++)
      gravity[dir] = gravity_in[dir];

    // frequency and bin size for sorting particles
    pp.query("particle_sorting_bin", particle_sorting_bin);
    sort_particle_int = -1;
    pp.query("sort_particle_int", sort_particle_int);

    // options for load balance
    pp.query("overload_tolerance",  overload_toler);
    pp.query("underload_tolerance", underload_toler);

    // Options to control initial projections (mostly we use these for
    // debugging)
    pp.query("initial_iterations", initial_iterations);
    pp.query("do_initial_proj", do_initial_proj);

    pp.query("test_tracer_conservation", test_tracer_conservation);

    pp.query("ntrac", ntrac);

    if (ntrac < 1)
      amrex::Abort("We currently require at least one tracer");

    // Scalar diffusion coefficients
    mu_s.resize(ntrac);
    for (int i = 0; i < ntrac; i++) mu_s[i] = 0.;
    pp.queryarr("mu_s", mu_s, 0, ntrac);

    amrex::Print() << "Scalar diffusion coefficients " << std::endl;
    for (int i = 0; i < ntrac; i++)
      amrex::Print() << "Tracer" << i << ":" << mu_s[i] << std::endl;

    if (test_tracer_conservation && !fluid.solve_tracer())
      amrex::Abort("No point in testing tracer conservation with fluid.solve_tracer"
          " = false");

    // At the moment, there is no relation between density and species
    //if (solve_species && !fluid.solve_density)
    //  amrex::Abort("Can't advect species mass fraction without advecting density");

    // At the moment, there is no relation between density and temperature
    //if (fluid.solve_enthalpy() )&&)&& !fluid.solve_density)
    //  amrex::Abort("Can't advect enthalpy without advecting density");

    // At the moment, there is no relation between density and tracer
    //if (fluid.solve_tracer && !fluid.solve_density)
    //  amrex::Abort("Can't advect tracer without advecting density");

    // control load balance
    // The default type is "KnapSack"; alternative is "SFC"
    pp.query("load_balance_type",      load_balance_type);
    pp.query("knapsack_weight_type",   knapsack_weight_type);
    pp.query("load_balance_fluid",     load_balance_fluid);
    pp.query("grid_pruning",           m_grid_pruning);


    // Include drag multiplier in projection. (False by default)
    pp.query("use_drag_coeff_in_proj_gp"        , m_use_drag_in_projection);

    // Redistribute before the nodal projection
    pp.query("redistribute_before_nodal_proj"   , m_redistribute_before_nodal_proj);

    // Redistribute after the nodal projection
    pp.query("redistribute_nodal_proj"          , m_redistribute_nodal_proj);

    // Threshold volfrac for correcting small cell velocity in the predictor and corrector
    pp.query("correction_small_volfrac"         , m_correction_small_volfrac);

    // Are we using MOL or Godunov?
    std::string l_advection_type = "Godunov";
    pp.query("advection_type"                   , l_advection_type);
    pp.query("use_ppm"                          , m_godunov_ppm);
    pp.query("godunov_use_forces_in_trans"      , m_godunov_use_forces_in_trans);
    pp.query("godunov_include_diff_in_forcing"  , m_godunov_include_diff_in_forcing);
    pp.query("use_mac_phi_in_godunov"           , m_use_mac_phi_in_godunov);
    pp.query("use_drag_in_godunov"              , m_use_drag_in_godunov);

    // agglomeration for GMG coarse levels
    pp.query("agg_grid_size", agg_grid_size);

    pp.query("redistribution_type"              , m_redistribution_type);
    if (m_redistribution_type != "NoRedist" &&
        m_redistribution_type != "FluxRedist" &&
        m_redistribution_type != "StateRedist")
      amrex::Abort("redistribution type must be FluxRedist, NoRedist or StateRedist");

    // Default to Godunov
    if(amrex::toLower(l_advection_type).compare("mol") == 0) {
      m_advection_type = AdvectionType::MOL;
    } else if(amrex::toLower(l_advection_type).compare("godunov") == 0) {
      m_advection_type = AdvectionType::Godunov;
    } else {
      amrex::Print() << "Advection type must be MOL or Godunov" << std::endl;
      amrex::Abort("Advection type must be MOL or Godunov");
    }

    if (advection_type() == AdvectionType::MOL) m_godunov_include_diff_in_forcing = false;

    // MOL: Explicit predictor / Crank_Nicolson corrector
    // Godunov: Implicit predictor / No corrector
    if (advection_type() == AdvectionType::MOL) {
      m_predictor_diff_type = DiffusionType::Explicit;
      m_corrector_diff_type = DiffusionType::Crank_Nicolson;
    } else {
      m_predictor_diff_type = DiffusionType::Implicit;
      m_corrector_diff_type = DiffusionType::Invalid;
    }

    // Default is true; should we use tensor solve instead of separate solves for each component?
    pp.query("use_tensor_solve",use_tensor_solve);
    pp.query("use_tensor_correction",use_tensor_correction);

    if (use_tensor_solve && use_tensor_correction) {
      amrex::Abort("We cannot have both use_tensor_solve and use_tensor_correction be true");
    }

    if (m_predictor_diff_type != DiffusionType::Implicit && use_tensor_correction) {
      amrex::Abort("We cannot have use_tensor_correction be true and diffusion type not Implicit");
    }

    // Options to control time stepping
    // MOL: default CFl = 0.5
    // Godunov: default CFL = 0.9
    if (advection_type() == AdvectionType::MOL) {
      m_cfl = 0.5;
    } else {
      m_cfl = 0.9;
    }
    pp.query("cfl", m_cfl);

    if (advection_type() == AdvectionType::MOL && m_cfl > 0.5) {
      amrex::Abort("We currently require cfl <= 0.5 when using the MOL advection scheme");
    }
    if (advection_type() != AdvectionType::MOL && m_cfl > 1.0) {
      amrex::Abort("We currently require cfl <= 1.0 when using the Godunov advection scheme");
    }

    // Verbosity and MLMG parameters are now ParmParse with "mac_proj" in the
    // inputs file
    // Examples: mac_proj.verbose = 1
    //           mac_proj.bottom_verbose = 1
    //           mac_proj.maxiter
    //           mac_proj.bottom_maxiter
    //           mac_proj.bottom_rtol
    //           mac_proj.bottom_atol
    //           mac_proj.bottom_solver
    // More info at "AMReX-Hydro/Projections/hydro_MacProjector.cpp"
    macproj_options = std::make_unique<MfixUtil::MLMGOptions>("mac_proj");

    // Checks for hypre namespace
    if (nodalproj_options->bottom_solver_type == "hypre" &&
          macproj_options->bottom_solver_type == "hypre") {
       std::string nodal_ns = nodalproj_options->hypre_namespace;
       std::string mac_ns = macproj_options->hypre_namespace;

       if (nodal_ns == "hypre" && mac_ns != "hypre")
          amrex::Abort("hypre namespace required for nodal projection");

       if (nodal_ns != "hypre" && mac_ns == "hypre")
          amrex::Abort("hypre namespace required for MAC projection");

       if ((nodal_ns == mac_ns) && nodal_ns != "hypre")
          amrex::Abort("same hypre namespace other than \"hypre\" not allowed for nodal and MAC projections");

    }

    AMREX_ALWAYS_ASSERT(load_balance_type.compare("KnapSack") == 0  ||
                        load_balance_type.compare("SFC") == 0 ||
                        load_balance_type.compare("Greedy") == 0);

    AMREX_ALWAYS_ASSERT(knapsack_weight_type.compare("RunTimeCosts") == 0 ||
                        knapsack_weight_type.compare("NumParticles") == 0);

    pp.query("dual_grid", dual_grid);

    if (load_balance_type.compare("KnapSack") == 0)
      pp.query("knapsack_nmax", knapsack_nmax);

    if (load_balance_type.compare("Greedy") == 0) {
      pp.query("greedy_dir",           greedy_dir);
      pp.query("greedy_3d",            greedy_3d);
      pp.query("greedy_min_grid_size", greedy_min_grid_size);
    }

    // fluid grids' distribution map
    pp.queryarr("pmap", pmap);

    // Parameters used be the level-set algorithm. Refer to LSFactory (or
    // mfix.H) for more details:
    //   -> refinement: how well resolved (fine) the (level-set/EB-facet)
    //                  grid needs to be (note: a fine level-set grid means
    //                  that distances and normals are computed accurately)
    //   -> pad:        how many (refined) grid points _outside_ the
    //                  problem domain the grid extends (avoids edge cases
    //                  in physical domain)
    pp.query("levelset__refinement", levelset_refinement);

    // Not needed here... the role of refining EB is filled with AMR level-set
    levelset_eb_refinement = 1;
    // Make sure that a coarsened level-set has a level-set pad of _at least_ 2;
    levelset_pad = 2*levelset_refinement;
    // Ensure that velocity_reconstruction has enough level-set to work off:
    // (2 => EB lives on the same grid resolution as fluid)
    levelset_eb_pad = amrex::max(2, levelset_pad);

    amrex::Print() << "Auto-generating level-set parameters:" << std::endl
                   << "eb_refinement = " << levelset_eb_refinement << std::endl
                   << "levelset_pad  = " << levelset_pad << std::endl
                   << "eb_pad        = " << levelset_eb_pad << std::endl;
  }

  if (m_dem.solve() || m_pic.solve())
  {
    ParmParse pp("particles");

    const int contains_size_x = pp.query("max_grid_size_x", particle_max_grid_size_x);
    const int contains_size_y = pp.query("max_grid_size_y", particle_max_grid_size_y);
    const int contains_size_z = pp.query("max_grid_size_z", particle_max_grid_size_z);

    if (contains_size_x && contains_size_y && contains_size_z) {

      particle_max_grid_size = IntVect(particle_max_grid_size_x,
                                       particle_max_grid_size_y,
                                       particle_max_grid_size_z);
    } else {

      AMREX_ALWAYS_ASSERT(!contains_size_x && !contains_size_y && !contains_size_z);

      // set the particles grid equal to the fluid grid
      particle_max_grid_size = max_grid_size[0];
    }

    // Keep particles that are initially touching the wall. Used by DEM tests.
    pp.query("removeOutOfRange", removeOutOfRange);
    pp.query("reduceGhostParticles", reduceGhostParticles);

    // distribution map for particle grids
    pp.queryarr("pmap", particle_pmap);
  }

  if ((m_dem.solve() || m_pic.solve()) && (!fluid.solve()))
  {
    if (m_timer.timestep_type() != MFIXTimer::TimestepType::Fixed)
      amrex::Abort("If running particle-only must specify a positive fixed_dt"
          " in the inputs file");
  }

  if ((m_dem.solve() || m_pic.solve()) && fluid.solve())
  {
    ParmParse pp("mfix");

    std::string drag_type = "None";
    pp.query("drag_type", drag_type);

    if (drag_type.compare("WenYu") == 0) {
      m_drag_type = DragType::WenYu;
    }
    else if (drag_type.compare("Gidaspow") == 0) {
      m_drag_type = DragType::Gidaspow;
    }
    else if (drag_type.compare("BVK2") == 0) {
      m_drag_type = DragType::BVK2;
    }
    else if (drag_type.compare("SyamOBrien") == 0) {
      m_drag_type = DragType::SyamOBrien;

      ParmParse ppSyamOBrien("mfix.SyamOBrien");
      int found_c1 = ppSyamOBrien.query("c1", m_SyamOBrien_coeff_c1);
      int found_d1 = ppSyamOBrien.query("d1", m_SyamOBrien_coeff_d1);

      if (!found_c1 || !found_d1 ) {
        std::string message;
        message  = " Error: Drag model SyamOBrien requires coefficients c1 and d1.\n";
        message += " Specify the following entries in the inputs file:\n";
        message += " mfix.SyamOBrien.c1 = amrex::Real\n";
        message += " mfix.SyamOBrien.d1 = amrex::Real\n";

        amrex::Print() << message;
        amrex::Abort(message);
      }

    }
    else if (drag_type.compare("UserDrag") == 0) {
      m_drag_type = DragType::UserDrag;
    }
    else {
      amrex::Abort("Don't know this drag_type!");
    }

    // Convection model type
    if (fluid.solve_enthalpy())
    {
      std::string convection_type = "None";

      if(!pp.contains("convection_type"))
      {
        if ( ParallelDescriptor::IOProcessor() )
          amrex::Warning("Convection type not specified in input file. "
              "Assuming RanzMarshall convection type");

        convection_type = "RanzMarshall";
      }
      else {
        pp.get("convection_type", convection_type);
      }

      if (convection_type.compare("RanzMarshall") == 0) {
        m_convection_type = ConvectionType::RanzMarshall;
      }
      else if (convection_type.compare("Gunn") == 0) {
        m_convection_type = ConvectionType::Gunn;
      }
      else if (convection_type.compare("NullConvection") == 0) {
        m_convection_type = ConvectionType::NullConvection;
      }
      else {
        amrex::Abort("Don't know this convection_type!");
      }
    }

    std::string deposition_scheme = "trilinear";
    pp.query("deposition_scheme", deposition_scheme);

    if (deposition_scheme.compare("trilinear") == 0) {
      m_deposition_scheme = DepositionScheme::trilinear;
    }
    else if (deposition_scheme.compare("trilinear-dpvm-square") == 0) {
      m_deposition_scheme = DepositionScheme::square_dpvm;
    }
    else if (deposition_scheme.compare("true-dpvm") == 0) {
      m_deposition_scheme = DepositionScheme::true_dpvm;
    }
    else if (deposition_scheme.compare("centroid") == 0) {
      m_deposition_scheme = DepositionScheme::centroid;
    }
    else {
      amrex::Abort("Don't know this deposition_scheme!");
    }

    m_max_solids_volume_fraction = 0.6;
    pp.query("max_solids_volume_fraction", m_max_solids_volume_fraction);
    //pp.query("close_pack", m_max_solids_volume_fraction);

    m_deposition_scale_factor = 1.;
    pp.query("deposition_scale_factor", m_deposition_scale_factor);

    m_deposition_diffusion_coeff = -1.;
    pp.query("deposition_diffusion_coeff", m_deposition_diffusion_coeff);
  }

  {
    ParmParse reports_pp("mfix.reports");

    reports_pp.query("mass_balance_int", m_rw->mass_balance_report_int);
    reports_pp.query("mass_balance_per_approx", m_rw->mass_balance_report_per_approx);

    if ((m_rw->mass_balance_report_int > 0 && m_rw->mass_balance_report_per_approx > 0) )
      amrex::Abort("Must choose only one of mass_balance_int or mass_balance_report_per_approx");

    // OnAdd check to turn off report if not solving species
    if (fluid.solve_species() && fluid.nspecies() >= 1) {
      m_rw->report_mass_balance = (m_rw->mass_balance_report_int > 0 ||
                                     m_rw->mass_balance_report_per_approx > 0);
    } else {
      m_rw->report_mass_balance = 0;
    }
  }

}


//! Tag using each EB level's volfrac. This requires that the `eb_levels` have
//! already been build.
void
mfix::ErrorEst (int lev,
                TagBoxArray & tags,
                Real /*time*/,
                int /*ngrow*/)
{
    if (ooo_debug) amrex::Print() << "ErrorEst" << std::endl;
    //___________________________________________________________________________
    // Tag all cells with volfrac \in (0, 1)
    MultiFab volfrac(grids[lev], dmap[lev], 1, 1);
    eb_levels[lev]->fillVolFrac(volfrac, geom[lev]);

    amrex::TagVolfrac(tags, volfrac);
}


void
mfix::Init (Real time,
            const bool init_fluid_grids)
{
  if (ooo_debug) amrex::Print() << "Init" << std::endl;
  m_rw->InitIOChkData();
  m_rw->InitIOPltData();

  // Note that finest_level = last level
  finest_level = nlev-1;

  // Fluid phase grids
  if (init_fluid_grids) {

    /****************************************************************************
     *                                                                          *
     * Generate levels using ErrorEst tagging.                                  *
     *                                                                          *
     ***************************************************************************/

    // This tells the AmrMesh class not to iterate when creating the initial
    // grid hierarchy
    SetIterateToFalse();

    // This tells the Cluster routine to use the new chopping routine which
    // rejects cuts if they don't improve the efficiency
    SetUseNewChop();

    /****************************************************************************
     *                                                                          *
     * MFIX-specific grid creation                                              *
     *                                                                          *
     ***************************************************************************/

    // Define coarse level BoxArray and DistributionMap
    const BoxArray& ba = MakeBaseGrids(geom[0].Domain(), max_grid_size[0]);

    DistributionMapping dm;

    if (pmap.empty())
      dm.define(ba, ParallelDescriptor::NProcs());
    else
      dm.define(pmap);

    MakeNewLevelFromScratch(0, time, ba, dm);

    for (int lev = 1; lev <= finest_level; lev++)
    {
       if (m_verbose > 0)
            std::cout << "Setting refined region at level " << lev
                      << " to " << grids[lev] << std::endl;

       MakeNewLevelFromScratch(lev, time, grids[lev], dmap[lev]);
    }
  }

  // Solids phase grids
  if (m_dem.solve() || m_pic.solve()) {

    BoxList             pbl;
    BoxArray            pba;
    DistributionMapping pdm;
    Vector<int>         pboxmap;  // map each particle box to its parent fluid box

    if (!dual_grid) {

      pbl = boxArray(0).boxList();
      pba = boxArray(0);
      pdm = DistributionMap(0);

    } else { // dual_grid is enabled

      // Define coarse level BoxArray and DistributionMap
      const BoxArray& ba = MakeBaseGrids(geom[0].Domain(), particle_max_grid_size);

      DistributionMapping dm;

      if (pmap.empty())
        dm.define(ba, ParallelDescriptor::NProcs());
      else
        dm.define(pmap);

      pbl = ba.boxList();
      pba = ba;
      pdm = dm;

      // chop grids to use all the gpus for particle generation
      if (ba.size() < ParallelDescriptor::NProcs()) {

        IntVect reduced_size = particle_max_grid_size;

        while (pbl.size() < ParallelDescriptor::NProcs()) {

          pbl.clear();
          pboxmap.clear();

          int maxdir = reduced_size.maxDir(false);
          reduced_size[maxdir] /= 2;

          for (auto i=0; i<pba.size(); i++) {

            BoxArray tmpba{pba[i]};
            tmpba.maxSize(reduced_size);
            pbl.join(tmpba.boxList());
            pboxmap.insert(pboxmap.end(), tmpba.size(), i);
          }
        }

        pba.define(pbl);
        pdm.define(pba, ParallelDescriptor::NProcs());
      }
    }

    pc = new MFIXParticleContainer(geom[0], pdm, pba, this->maxLevel()+1,
                                   m_initial_conditions, m_boundary_conditions,
                                   solids, m_dem, m_pic, fluid, reactions);

    pc->setSortingBinSizes(IntVect(particle_sorting_bin));

    if (load_balance_type.compare("Greedy") == 0)
      pc->setGreedyRegrid(greedy_dir, greedy_3d, greedy_min_grid_size);

    if (!pboxmap.empty())
      pc->setParticleFluidGridMap(pboxmap);

    // Updating m_rw pc pointer is needed since mfix pc has changed
    m_rw->set_pc(pc);
  }

  /****************************************************************************
   *                                                                          *
   * MFIX-Specific Initialization                                             *
   *                                                                          *
   ***************************************************************************/

  // ******************************************************
  // We only do these at level 0
  // ******************************************************

  for (int lev = 0; lev < nlev; lev++)
      mfix_set_bc_type(lev,nghost_state());
}


void
mfix::PruneBaseGrids (BoxArray &ba) const
{
    // Use 1 ghost layer
    EBDataCollection ebdc(*eb_levels[0], geom[0],
          ba, DistributionMapping{ba}, {1}, EBSupport::basic);

    const auto &cflag = ebdc.getMultiEBCellFlagFab();
    Vector<Box> uncovered;

    for (MFIter mfi(cflag); mfi.isValid(); ++mfi)
    {
        FabType t = cflag[mfi].getType();
        const Box& vbx = mfi.validbox();
        if (t != FabType::covered) {
            uncovered.push_back(vbx);
        }
    }

    amrex::AllGatherBoxes(uncovered);
    ba = BoxArray(BoxList(std::move(uncovered)));
}


BoxArray
mfix::MakeBaseGrids (const Box& domain,
                     const IntVect& grid_sizes) const
{
    if (ooo_debug)
      amrex::Print() << "MakeBaseGrids" << std::endl;

    BoxArray ba(domain);

    ba.maxSize(grid_sizes);

    if (m_grid_pruning) {
       PruneBaseGrids(ba);
    }

    // We only call ChopGrids if dividing up the grid using max_grid_size didn't
    //    create enough grids to have at least one grid per processor.
    // This option is controlled by "refine_grid_layout" which defaults to true.
    if (refine_grid_layout && ba.size() < ParallelDescriptor::NProcs())
        ChopGrids(geom[0].Domain(), ba, ParallelDescriptor::NProcs());

    if (ba == grids[0]) {
        ba = grids[0];  // to avoid duplicates
    }

    amrex::Print() << "In MakeBaseGrids: BA HAS " << ba.size() << " GRIDS " << std::endl;
    return ba;
}


void
mfix::ChopGrids (const Box& domain,
                 BoxArray& ba,
                 int target_size) const
{
    if (ooo_debug) amrex::Print() << "ChopGrids" << std::endl;
    if ( ParallelDescriptor::IOProcessor() )
       amrex::Warning("Using max_grid_size only did not make enough grids for the number of processors");

    // Here we hard-wire the maximum number of times we divide the boxes.
    int n = 10;

    // Here we hard-wire the minimum size in any one direction the boxes can be
    int min_grid_size = 4;

    IntVect chunk(domain.length(0),domain.length(1),domain.length(2));

    int j = -1;
    for (int cnt = 1; cnt <= n; ++cnt)
    {
        if (chunk[0] >= chunk[1] && chunk[0] >= chunk[2])
        {
            j = 0;
        }
        else if (chunk[1] >= chunk[0] && chunk[1] >= chunk[2])
        {
            j = 1;
        }
        else if (chunk[2] >= chunk[0] && chunk[2] >= chunk[1])
        {
            j = 2;
        }
        chunk[j] /= 2;

        if (chunk[j] >= min_grid_size)
        {
            ba.maxSize(chunk);
        }
        else
        {
            // chunk[j] was the biggest chunk -- if this is too small then we're done
            if ( ParallelDescriptor::IOProcessor() )
               amrex::Warning("ChopGrids was unable to make enough grids for the number of processors");
            return;
        }

        // Test if we now have enough grids
        if (ba.size() >= target_size)
          return;
    }
}


void
mfix::MakeNewLevelFromScratch (int lev, Real /*time*/,
                               const BoxArray& new_grids,
                               const DistributionMapping& new_dmap)
{
    if (ooo_debug)
      amrex::Print() << "MakeNewLevelFromScratch" << std::endl;

    if (m_verbose > 0)
    {
        std::cout << "MAKING NEW LEVEL " << lev << std::endl;
        std::cout << "WITH BOX ARRAY   " << new_grids << std::endl;
    }

    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

    if (m_verbose > 0) {
      amrex::Print() << "SETTING NEW GRIDS IN MAKE NEW LEVEL " << new_grids << std::endl;
      amrex::Print() << "SETTING NEW DMAP IN MAKE NEW LEVEL " << new_dmap << std::endl;
    }

    macproj = std::make_unique<Hydro::MacProjector>(Geom(0,finest_level),
                                   MLMG::Location::FaceCentroid,  // Location of mac_vec
                                   MLMG::Location::FaceCentroid,  // Location of beta
                                   MLMG::Location::CellCenter,    // Location of solution variable phi
                                   MLMG::Location::CellCentroid);// Location of MAC RHS

    // This is being done by mfix::make_eb_geometry,
    // otherwise it would be done here
    if (lev == 0)
      bc_list.MakeBCArrays(nghost_state(), ooo_debug, geom);
}


void
mfix::ReMakeNewLevelFromScratch (int lev,
                                 const BoxArray & new_grids,
                                 const DistributionMapping & new_dmap)
{
    if (ooo_debug) amrex::Print() << "ReMakeNewLevelFromScratch" << std::endl;
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

    if (lev == 0)
      bc_list.MakeBCArrays(nghost_state(), ooo_debug, geom);

    // We need to re-fill these arrays for the larger domain (after replication).
    mfix_set_bc_type(lev,nghost_state());
}



void
mfix::InitLevelData (Real /*time*/)
{
    if (ooo_debug) amrex::Print() << "InitLevelData" << std::endl;
    // Allocate the fluid data, NOTE: this depends on the ebfactories.
    if (fluid.solve())
       for (int lev = 0; lev < nlev; lev++)
          AllocateArrays(lev);

    if (m_rw->only_print_grid_report) {
       return;
    }

    // Allocate the particle data
    if (m_dem.solve() || m_pic.solve())
    {
      Real strt_init_part = ParallelDescriptor::second();

      pc->AllocData();

      if (m_initial_conditions.AutoParticleInit()) {
         amrex::Print() << "Auto generating particles ..." << std::endl;
         // Particles are generated on level 0
         pc->InitParticlesAuto(particle_ebfactory[0].get());

      } else {
        amrex::Print() << "Reading particles from particle_input.dat ..." << std::endl;
        pc->InitParticlesAscii("particle_input.dat");
      }

      pc->Redistribute();

      Real end_init_part = ParallelDescriptor::second() - strt_init_part;
      ParallelDescriptor::ReduceRealMax(end_init_part, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "Time spent in initializing particles " << end_init_part << std::endl;
    }

    // Used in load balancing
    if (m_dem.solve() || m_pic.solve())
    {
      for (int lev(0); lev < particle_cost.size(); lev++) {
        if (particle_cost[lev] != nullptr)  delete particle_cost[lev];
        if (particle_proc[lev] != nullptr)  delete particle_proc[lev];
      }

      particle_cost.clear();
      particle_cost.resize(nlev, nullptr);
      particle_proc.clear();
      particle_proc.resize(nlev, nullptr);

      const Real proc = static_cast<Real>(ParallelDescriptor::MyProc());

      for (int lev = 0; lev < nlev; lev++)
      {
        particle_cost[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                          pc->ParticleDistributionMap(lev), 1, 0);
        particle_cost[lev]->setVal(0.0);
        particle_proc[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                          pc->ParticleDistributionMap(lev), 1, 0);
        particle_proc[lev]->setVal(proc);
      }
    }

    // Used in load balancing
    if (fluid.solve())
    {
      for (int lev(0); lev < fluid_cost.size(); lev++) {
        if (fluid_cost[lev] != nullptr)  delete fluid_cost[lev];
        if (fluid_proc[lev] != nullptr)  delete fluid_proc[lev];
      }

      fluid_cost.clear();
      fluid_cost.resize(nlev, nullptr);
      fluid_proc.clear();
      fluid_proc.resize(nlev, nullptr);

      const Real proc = static_cast<Real>(ParallelDescriptor::MyProc());

      for (int lev = 0; lev < nlev; lev++)
      {
        fluid_cost[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0);
        fluid_cost[lev]->setVal(0.0);
        fluid_proc[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0);
        fluid_proc[lev]->setVal(proc);
      }
    }
}


void
mfix::PostInit (Real& dt,
                const Real time,
                int is_restarting,
                Real stop_time)
{
    if (ooo_debug) amrex::Print() << "PostInit" << std::endl;

    if (m_dem.solve() || m_pic.solve())
    {
        // Auto generated particles may be out of the domain. This call will
        // remove them. Note that this has to occur after the EB geometry is
        // created. if (particle_init_type == "Auto" && !is_restarting &&
        // particle_ebfactory[finest_level])

      if (removeOutOfRange)
        {

          if ((nlev == 1) &&
              ((!is_restarting || m_run_type == RunType::PIC2DEM) && particle_ebfactory[finest_level]))
          {
            //___________________________________________________________________
            // Only 1 refined level-set

            Print() << "Clean up auto-generated particles.\n" << std::endl;

            const MultiFab * ls_data = level_sets[1].get();
            iMultiFab ls_valid(ls_data->boxArray(), ls_data->DistributionMap(),
                               ls_data->nComp(), ls_data->nGrow());

            pc->RemoveOutOfRange(finest_level, particle_ebfactory[finest_level].get(),
                                 ls_data, levelset_refinement);
          }
          else if (!is_restarting && particle_ebfactory[finest_level])
          {
            //___________________________________________________________________
            // Multi-level everything

            Print() << "Clean up auto-generated particles.\n" << std::endl;

            for (int ilev = 0; ilev < nlev; ilev ++)
            {
              const MultiFab * ls_data = level_sets[ilev].get();
              iMultiFab ls_valid(ls_data->boxArray(), ls_data->DistributionMap(),
                                 ls_data->nComp(), ls_data->nGrow());

              pc->RemoveOutOfRange(ilev, particle_ebfactory[ilev].get(),
                                   ls_data, 1);
            }
          }

        }
    }

    if (m_run_type != RunType::PIC2DEM) {

      if (m_dem.solve() || m_pic.solve()) {

        pc->setSortInt(sort_particle_int);
        if (m_dem.solve()) {
            pc->MFIX_PC_InitCollisionParams();
            pc->setReduceGhostParticles(reduceGhostParticles);
        }

        if (dual_grid) { Regrid(); }

        if (!is_restarting) { pc->InitParticlesRuntimeVariables(fluid.solve_enthalpy()); }

        if (!fluid.solve()) { dt = m_timer.dt(); }
      }
    }

    if (fluid.solve()) {
        mfix_init_fluid(is_restarting, time, dt, stop_time);
    }

    for (int lev = 0; lev < nlev; lev++) {

      if(m_dem.solve() || m_pic.solve()) {

        m_boundary_conditions.calc_bc_areas(lev,
           pc->ParticleBoxArray(lev),
           pc->ParticleDistributionMap(lev),
           particle_ebfactory[lev].get());

      } else {

        m_boundary_conditions.calc_bc_areas(lev,
            grids[lev], dmap[lev], ebfactory[lev].get());
      }
    }

    // Call user-defined subroutine to set constants, check data, etc.
    if (call_udf) mfix_usr0();
}


void
mfix::mfix_init_fluid (int is_restarting,
                       const Real time,
                       Real dt,
                       Real stop_time)
{
    if (ooo_debug) amrex::Print() << "mfix_init_fluid" << std::endl;

    Real xlen = geom[0].ProbHi(0) - geom[0].ProbLo(0);
    Real ylen = geom[0].ProbHi(1) - geom[0].ProbLo(1);
    Real zlen = geom[0].ProbHi(2) - geom[0].ProbLo(2);

    // Here we set bc values for p and u,v,w before the IC's are set
    mfix_set_bc0();

    for (int lev = 0; lev < nlev; lev++)
    {
       Box domain(geom[lev].Domain());

       MultiFab& ep_g = *(m_leveldata[lev]->ep_g);

       Real dx = geom[lev].CellSize(0);
       Real dy = geom[lev].CellSize(1);
       Real dz = geom[lev].CellSize(2);

       const GpuArray<Real, 3> plo = geom[lev].ProbLoArray();

       LevelData& ld = *m_leveldata[lev];

       // We deliberately don't tile this loop since we will be looping
       //    over bc's on faces and it makes more sense to do this one grid at a time
       for (MFIter mfi(ep_g, false); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          const Box& sbx = ep_g[mfi].box();

          if (is_restarting) {
            init_fluid_parameters(bx, mfi, ld, fluid);
          } else {
            init_fluid(sbx, bx, domain, mfi, ld, dx, dy, dz, xlen, ylen, zlen, plo,
                       test_tracer_conservation, m_initial_conditions, fluid);
          }
       }
    }

    mfix_set_p0();

    // Here we re-set the bc values for p and u,v,w just in case init_fluid
    //      over-wrote some of the bc values with ic values
    mfix_set_bc0();

    // Make sure to fill the "old state" before we start.
    for (int lev = 0; lev < nlev; lev++)
    {
       LevelData& ld = *m_leveldata[lev];

       // TODO check if commenting the following is correct
       //MultiFab::Copy(*m_leveldata[lev]->vel_go,  *m_leveldata[lev]->vel_g, 0, 0, 3, 0);

       MultiFab::Copy(*ld.ro_go,  *ld.ro_g, 0, 0, 1, ld.ro_g->nGrow());
       MultiFab::Copy(*ld.trac_o, *ld.trac, 0, 0, 1, ld.trac->nGrow());

       if (fluid.solve_enthalpy()) {
         MultiFab::Copy(*ld.T_go, *ld.T_g, 0, 0, 1, ld.T_g->nGrow());
         MultiFab::Copy(*ld.h_go, *ld.h_g, 0, 0, 1, ld.h_g->nGrow());
       }

       if (fluid.solve_species()) {
         MultiFab::Copy(*ld.X_gko, *ld.X_gk, 0, 0, fluid.nspecies(), ld.X_gk->nGrow());
       }

       if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem &&
           fluid.solve_enthalpy()) {
         *(ld.thermodynamic_p_go) = *(ld.thermodynamic_p_g);
       }
    }

    if (is_restarting == 0)
    {
      // Just for reference, we compute the volume inside the EB walls (as if
      // there were no particles)
      m_leveldata[0]->ep_g->setVal(1.0);

      const Real* dx = geom[0].CellSize();
      const Real cell_volume = dx[0] * dx[1] * dx[2];

      sum_vol_orig = Utils::volWgtSum(0, *(m_leveldata[0]->ep_g), 0, ebfactory);

      Print() << "Enclosed domain volume is   " << cell_volume * sum_vol_orig << std::endl;

      Real domain_vol = sum_vol_orig;

      // Now initialize the volume fraction ep_g before the first projection
      mfix_calc_volume_fraction(time, sum_vol_orig);
      Print() << "Setting original sum_vol to " << cell_volume * sum_vol_orig << std::endl;

      Print() << "Difference is   " << cell_volume * (domain_vol - sum_vol_orig) << std::endl;

      m_boundary_conditions.set_density_bcs(time, get_ro_g());
      m_boundary_conditions.set_density_bcs(time, get_ro_g_old());

      m_boundary_conditions.set_tracer_bcs(time, fluid, get_trac());
      m_boundary_conditions.set_tracer_bcs(time, fluid, get_trac_old());

      if (fluid.solve_enthalpy()) {
        m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g());
        m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g_old());
      }

      if (fluid.solve_enthalpy()) {
        m_boundary_conditions.set_enthalpy_bcs(time, fluid,get_h_g());
        m_boundary_conditions.set_enthalpy_bcs(time, fluid,get_h_g_old());
      }

      if (fluid.solve_species()) {
        m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk());
        m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk_old());
      }

      fillpatch_all(get_vel_g(), get_ro_g(), get_h_g(), get_trac(), get_X_gk(), time);
      InitialRedistribution(time);

      // Project the initial velocity field
      if (do_initial_proj)
        mfix_project_velocity();

      // Iterate to compute the initial pressure
      if (initial_iterations > 0)
        mfix_initial_iterations(dt, stop_time);

      m_rw->InitMassBalance(fluid.nspecies());

    }
    else
    {
      const int dir_bc = 1;
      m_boundary_conditions.set_epg_bcs(time, get_ep_g(), dir_bc);

      const Real* dx = geom[0].CellSize();
      const Real cell_volume = dx[0] * dx[1] * dx[2];

      if (m_run_type != RunType::PIC2DEM) {

        //Calculation of sum_vol_orig for a restarting point
        sum_vol_orig = Utils::volWgtSum(0, *(m_leveldata[0]->ep_g), 0, ebfactory);

        Print() << "Setting original sum_vol to " << cell_volume * sum_vol_orig << std::endl;
      }
    }
}


void
mfix::mfix_set_bc0 ()
{
    if (ooo_debug) amrex::Print() << "mfix_set_bc0" << std::endl;

    Real time = 0.0;

    if (fluid.solve_enthalpy()) {
      m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g());
      m_boundary_conditions.set_enthalpy_bcs(time, fluid,get_h_g());
    }

    if (fluid.solve_species())
      m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk());

   // Put velocity Dirichlet bc's on faces
   int extrap_dir_bcs = 0;

   m_boundary_conditions.set_velocity_bcs(time, get_vel_g(), extrap_dir_bcs);
}


void
mfix::mfix_set_p0 ()
{
  if (ooo_debug) amrex::Print() << "mfix_set_p0" << std::endl;

  IntVect press_per = IntVect(geom[0].isPeriodic(0),geom[0].isPeriodic(1),geom[0].isPeriodic(2));

  // Here we set a separate periodicity flag for p0_g because when we use
  // pressure drop (delp) boundary conditions we fill all variables *except* p0
  // periodically
  if (m_boundary_conditions.delp_dir() > -1)
    press_per[m_boundary_conditions.delp_dir()] = 0;

  p0_periodicity = Periodicity(press_per);

  // Initialize gp0 to 0
  gp0[0] = 0.0;
  gp0[1] = 0.0;
  gp0[2] = 0.0;

  for (int lev = 0; lev < nlev; lev++)
  {
     Box domain(geom[lev].Domain());

     // We put this outside the MFIter loop because we need gp0 even on ranks with no boxes
     // because we will use it in computing dt separately on every rank
     set_gp0(lev, domain);

     // We deliberately don't tile this loop since we will be looping
     //    over bc's on faces and it makes more sense to do this one grid at a time
     for (MFIter mfi(*(m_leveldata[lev]->ep_g), false); mfi.isValid(); ++mfi)
     {
       const Box& bx = mfi.validbox();

       set_p0 (bx, &mfi, lev, domain );
     }
   }
}


void
mfix::mfix_set_ls_near_inflow ()
{
    if (ooo_debug) amrex::Print() << "mfix_ls_near_inflow" << std::endl;
    // This function is a bit Wonky... TODO: figure out why we need + nghost
    // (it's late at the moment, so I can't figure it it out right now...)
    const int levelset_nghost = levelset_eb_pad + nghost_state();

    if (nlev > 1)
    {
        // The level set is always at the same level of refinement as the
        // boundary condition arrays
        int n = 1;

        for (int lev = 0; lev < nlev; lev++)
        {
            Box domain(geom[lev].Domain());

            MultiFab* ls_phi = level_sets[lev].get();
            const Real* dx   = geom[lev].CellSize();

            // Don't tile this
            for (MFIter mfi(*ls_phi, false); mfi.isValid(); ++mfi)
            {
                FArrayBox & ls_fab = (* ls_phi)[mfi];

                set_ls_inflow(lev, ls_fab, domain, &levelset_nghost, n, dx);
            }
        }
    }
    else
    {
        // Here we assume that the particles and fluid only exist on one level
        int lev = 0;
        int lev_ref = 1;

        // ... but that the level set may be at a finer resolution.
        int n = levelset_refinement;
        {
            Box domain(geom[lev].Domain());
            const Real* dx   = geom[lev].CellSize();
            MultiFab* ls_phi = level_sets[lev_ref].get();

            // Don't tile this
            for (MFIter mfi(*ls_phi, false); mfi.isValid(); ++mfi)
            {
                FArrayBox& ls_fab = (*ls_phi)[mfi];

                set_ls_inflow( lev, ls_fab, domain, &levelset_nghost, n, dx);
            }
        }

        // ... now also "fix" the level-zero level-set (for consistency with
        // multi-level levelset)
        n = 1;
        {
            Box domain(geom[lev].Domain());
            const Real* dx   = geom[lev].CellSize();
            MultiFab* ls_phi = level_sets[lev].get();

            // Don't tile this
            for (MFIter mfi(*ls_phi, false); mfi.isValid(); ++mfi)
            {
                FArrayBox& ls_fab = (* ls_phi)[mfi];

                set_ls_inflow(lev, ls_fab, domain, &levelset_nghost, n, dx);
            }
        }
    }
}
