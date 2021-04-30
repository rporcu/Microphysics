#include <AMReX_ParmParse.H>

#include <mfix.H>

#include <mfix_init_fluid.H>

#include <AMReX_EBAmrUtil.H>

#include <mfix_regions_parms.H>
#include <mfix_bc_parms.H>
#include <mfix_ic_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_solids_parms.H>
#include <mfix_species_parms.H>

using MFIXParIter = MFIXParticleContainer::MFIXParIter;
using PairIndex = MFIXParticleContainer::PairIndex;

void
mfix::InitParams ()
{
  if (ooo_debug) amrex::Print() << "InitParams" << std::endl;

  // Read and process species, fluid and DEM particle model options.
  SPECIES::Initialize();
  fluid.Initialize();
  solids.Initialize();

  enthalpy_source = solids.enthalpy_source;

  BL_ASSERT(fluid.nspecies <= SPECIES::NMAX);
  BL_ASSERT(solids.nspecies <= SPECIES::NMAX);

  // Read and process chemical reactions inputs.
  REACTIONS::Initialize();
  initialize_chem_reactions();

  BL_ASSERT(REACTIONS::nreactions <= REACTIONS::NMAX);

  DEM::Initialize();
  PIC::Initialize();

  // Need to do this -- We might want to move this to a BC class once we have
  // one
  bcs_X.resize(fluid.nspecies);

  // Read in regions, initial and boundary conditions. Note that
  // regions need to be processed first as they define the
  // physical extents of ICs and BCs.
  REGIONS::Initialize();
  IC::Initialize(fluid, solids);
  BC::Initialize(geom[0], fluid, solids);

  // set n_error_buf (used in AmrMesh) to default (can overwrite later)
  for (int i = 0; i < n_error_buf.size(); i++)
    n_error_buf[i] = {8,8,8};

  {
    ParmParse pp("mfix");

    // Options to control time stepping
    pp.query("cfl", m_cfl);

    fixed_dt = -1.;
    pp.query("fixed_dt", fixed_dt);
    pp.query("dt_min", dt_min);
    pp.query("dt_max", dt_max);

    // Verbosity and MLMG parameters are now ParmParse with "nodal_proj" in the
    // inputs file
    // Examples: nodal_proj.verbose = 1
    //           nodal_proj.bottom_verbose = 1
    //           nodal_proj.maxiter
    //           nodal_proj.bottom_maxiter
    //           nodal_proj.bottom_rtol
    //           nodal_proj.bottom_atol
    //           nodal_proj.bottom_solver
    // More info at "amrex/Src/LinearSolvers/Projections/AMReX_NodalProjector.cpp"
    ParmParse pp_nodal("nodal_proj");
    // Options to control MLMG behavior
    pp_nodal.query("mg_rtol", nodal_mg_rtol);
    pp_nodal.query("mg_atol", nodal_mg_atol);
    pp_nodal.query("mg_max_coarsening_level", nodal_mg_max_coarsening_level);

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

    // Flag to envoke UDFs
    bool call_usr_bool = false;
    pp.query("call_usr", call_usr_bool);
    call_udf = call_usr_bool ? 1 : 0; // Set global flag

    Array<Real,3> gravity_in{0.0, 0.0, 0.0};
    pp.get("gravity", gravity_in);
    for (int dir = 0; dir < 3; dir++)
      gravity[dir] = gravity_in[dir];

    // The default type is "AsciiFile" but we can over-write that in the inputs
    // file with "Random"
    pp.query("particle_init_type", particle_init_type);

    Array<int,3> sorting_bin{0, 0, 0};
    pp.query("particle_sorting_bin", sorting_bin);
    particle_sorting_bin = IntVect(sorting_bin);

    // Options to control initial projections (mostly we use these for
    // debugging)
    pp.query("initial_iterations", initial_iterations);
    pp.query("do_initial_proj", do_initial_proj);

    pp.query("advect_density", advect_density);
    pp.query("advect_tracer" , advect_tracer);
    pp.query("advect_enthalpy", advect_enthalpy);
    pp.query("test_tracer_conservation", test_tracer_conservation);

    // Set the mfix class flag equal to the FLUID parameter
    advect_fluid_species = fluid.solve_species;

    // We can still turn it off explicitly even if we passed species inputs
    pp.query("advect_fluid_species", advect_fluid_species);

    // Set the mfix class flag equal to the REACTIONS parameter
    solve_reactions = REACTIONS::solve && SPECIES::solve;

    // We can still turn it off explicitly even if we passed stoichiometry inputs
    pp.query("solve_reactions", solve_reactions);

    if (advect_fluid_species)
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.solve_species,
          "Advect fluid species flag is on but no fluid species were provided");

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

    if (test_tracer_conservation && !advect_tracer)
      amrex::Abort("No point in testing tracer conservation with advect_tracer"
          " = false");

    // At the moment, there is no relation between density and species
    //if (advect_fluid_species && !advect_density)
    //  amrex::Abort("Can't advect species mass fraction without advecting density");

    // At the moment, there is no relation between density and temperature
    //if (advect_enthalpy && !advect_density)
    //  amrex::Abort("Can't advect enthalpy without advecting density");

    if (advect_tracer && !advect_density)
      amrex::Abort("Can't advect tracer without advecting density");

    // The default type is "KnapSack"; alternative is "SFC"
    pp.query("load_balance_type", load_balance_type);
    pp.query("knapsack_weight_type", knapsack_weight_type);
    pp.query("load_balance_fluid", load_balance_fluid);


    // Include drag multiplier in projection. (False by default)
    pp.query("use_drag_coeff_in_proj_gp"        , m_use_drag_in_projection);

    // Are we using MOL or Godunov?
    std::string l_advection_type = "MOL";
    pp.query("advection_type"                   , l_advection_type);
    pp.query("use_ppm"                          , m_godunov_ppm);
    pp.query("godunov_use_forces_in_trans"      , m_godunov_use_forces_in_trans);
    pp.query("godunov_include_diff_in_forcing"  , m_godunov_include_diff_in_forcing);
    pp.query("use_mac_phi_in_godunov"           , m_use_mac_phi_in_godunov);
    pp.query("use_drag_in_godunov"              , m_use_drag_in_godunov);


    pp.query("redistribution_type"              , m_redistribution_type);
    if (m_redistribution_type != "NoRedist" &&
        m_redistribution_type != "FluxRedist" &&
        m_redistribution_type != "StateRedist")
      amrex::Abort("redistribution type must be FluxRedist, NoRedist or StateRedist");


    // Default to MOL
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
    // Godunov: Crank_Nicolson
    if (advection_type() == AdvectionType::MOL) {
      m_predictor_diff_type = DiffusionType::Explicit;
      m_corrector_diff_type = DiffusionType::Crank_Nicolson;
    } else {
      m_predictor_diff_type = DiffusionType::Crank_Nicolson;
      m_corrector_diff_type = DiffusionType::Crank_Nicolson;
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
    // More info at "amrex/Src/LinearSolvers/Projections/AMReX_MacProjector.cpp"

    ParmParse pp_mac("mac_proj");
    pp_mac.query("mg_rtol", mac_mg_rtol);
    pp_mac.query("mg_atol", mac_mg_atol);
    pp_mac.query("mg_max_coarsening_level", mac_mg_max_coarsening_level);

    AMREX_ALWAYS_ASSERT(load_balance_type.compare("KnapSack") == 0  ||
                        load_balance_type.compare("SFC") == 0);

    AMREX_ALWAYS_ASSERT(knapsack_weight_type.compare("RunTimeCosts") == 0 ||
                        knapsack_weight_type.compare("NumParticles") == 0);

    ParmParse amr_pp("amr");
    amr_pp.query("dual_grid", dual_grid);

    if (load_balance_type.compare("KnapSack") == 0)
      pp.query("knapsack_nmax", knapsack_nmax);

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

  if (DEM::solve || PIC::solve)
  {
    ParmParse pp("particles");

    pp.query("max_grid_size_x", particle_max_grid_size_x);
    pp.query("max_grid_size_y", particle_max_grid_size_y);
    pp.query("max_grid_size_z", particle_max_grid_size_z);

    // Keep particles that are initially touching the wall. Used by DEM tests.
    pp.query("removeOutOfRange", removeOutOfRange);
    
    // distribution map for particle grids
    pp.queryarr("pmap", particle_pmap);
  }

  if ((DEM::solve || PIC::solve) && (!fluid.solve))
  {
    if (fixed_dt <= 0.0)
      amrex::Abort("If running particle-only must specify a positive fixed_dt"
          " in the inputs file");
  }

  if ((DEM::solve || PIC::solve) && fluid.solve)
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
    else if (drag_type.compare("UserDrag") == 0) {
      m_drag_type = DragType::UserDrag;
    }
    else {
      amrex::Abort("Don't know this drag_type!");
    }

    // Convection model type
    if (advect_enthalpy)
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
      else {
        amrex::Abort("Don't know this convection_type!");
      }
    }

    // Constraint type
    {
      std::string idealgas_constraint_type = "none";
      pp.query("idealgas_constraint", idealgas_constraint_type);
      idealgas_constraint_type = amrex::toLower(idealgas_constraint_type);

      if (idealgas_constraint_type.compare("none") == 0) {
        m_idealgas_constraint = IdealGasConstraint::None;
      }
      else if (idealgas_constraint_type.compare("opensystem") == 0 ||
               idealgas_constraint_type.compare("open_system") == 0) {
        m_idealgas_constraint = IdealGasConstraint::OpenSystem;
      }
      else if (idealgas_constraint_type.compare("closedsystem") == 0 ||
               idealgas_constraint_type.compare("closed_system") == 0) {
        m_idealgas_constraint = IdealGasConstraint::ClosedSystem;
      }
      else {
        amrex::Abort("Don't know this idealgas_constraint_type!");
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
    ParmParse amr_pp("amr");

    amr_pp.query("restart_from_cold_flow", restart_from_cold_flow);

    amr_pp.query("plot_int", plot_int);
    amr_pp.query("plot_per_exact", plot_per_exact);
    amr_pp.query("plot_per_approx", plot_per_approx);

    if ((plot_int       > 0 && plot_per_exact  > 0) ||
        (plot_int       > 0 && plot_per_approx > 0) ||
        (plot_per_exact > 0 && plot_per_approx > 0) )
      amrex::Abort("Must choose only one of plot_int or plot_per_exact or plot_per_approx");

    amr_pp.queryarr("avg_p_g", avg_p_g);
    amr_pp.queryarr("avg_ep_g", avg_ep_g);
    amr_pp.queryarr("avg_vel_g", avg_vel_g);
    amr_pp.queryarr("avg_T_g", avg_T_g);

    amr_pp.queryarr("avg_vel_p", avg_vel_p);

    amr_pp.queryarr("avg_T_p", avg_T_p);

    // Regions geometry
    amr_pp.queryarr("avg_region_x_e", avg_region_x_e);
    amr_pp.queryarr("avg_region_x_w", avg_region_x_w);
    amr_pp.queryarr("avg_region_y_n", avg_region_y_n);
    amr_pp.queryarr("avg_region_y_s", avg_region_y_s);
    amr_pp.queryarr("avg_region_z_t", avg_region_z_t);
    amr_pp.queryarr("avg_region_z_b", avg_region_z_b);
  }
}


//! Tag using each EB level's volfrac. This requires that the `eb_levels` have
//! already been build.
void mfix::ErrorEst (int lev, TagBoxArray & tags, Real time, int ngrow)
{
    if (ooo_debug) amrex::Print() << "ErrorEst" << std::endl;
    //___________________________________________________________________________
    // Tag all cells with volfrac \in (0, 1)
    MultiFab volfrac(grids[lev], dmap[lev], 1, 1);
    eb_levels[lev]->fillVolFrac(volfrac, geom[lev]);

    amrex::TagVolfrac(tags, volfrac);
}


void mfix::Init (Real time)
{
    if (ooo_debug) amrex::Print() << "Init" << std::endl;
    InitIOChkData();
    InitIOPltData();

    // Note that finest_level = last level
    finest_level = nlev-1;

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
    const BoxArray& ba = MakeBaseGrids();
    DistributionMapping dm;
    if (pmap.empty())
      dm.define(ba, ParallelDescriptor::NProcs());
    else
      dm.define(pmap);
    // DistributionMapping dm(ba, ParallelDescriptor::NProcs());
    MakeNewLevelFromScratch(0, time, ba, dm);


    for (int lev = 1; lev <= finest_level; lev++)
    {
       if (m_verbose > 0)
            std::cout << "Setting refined region at level " << lev
                      << " to " << grids[lev] << std::endl;

       MakeNewLevelFromScratch(lev, time, grids[lev], dmap[lev]);
    }

    /****************************************************************************
     *                                                                          *
     * Create particle container using mfix::ParGDB                             *
     *                                                                          *
     ***************************************************************************/

    if (DEM::solve || PIC::solve) {
      pc = new MFIXParticleContainer(this, solids);
      pc->setSortingBinSizes(IntVect(particle_sorting_bin));
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


BoxArray mfix::MakeBaseGrids () const
{
    if (ooo_debug) amrex::Print() << "MakeBaseGrids" << std::endl;
    BoxArray ba(geom[0].Domain());

    ba.maxSize(max_grid_size[0]);

    // We only call ChopGrids if dividing up the grid using max_grid_size didn't
    //    create enough grids to have at least one grid per processor.
    // This option is controlled by "refine_grid_layout" which defaults to true.
    if ( refine_grid_layout &&
         ba.size() < ParallelDescriptor::NProcs() )
           ChopGrids(geom[0].Domain(), ba, ParallelDescriptor::NProcs());

    if (ba == grids[0]) {
        ba = grids[0];  // to avoid duplicates
    }
    amrex::Print() << "In MakeBaseGrids: BA HAS " << ba.size() << " GRIDS " << std::endl;
    return ba;
}


void mfix::ChopGrids (const Box& domain, BoxArray& ba, int target_size) const
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
        if (ba.size() >= target_size) return;
    }
}


void mfix::MakeNewLevelFromScratch (int lev, Real time,
                                    const BoxArray& new_grids,
                                    const DistributionMapping& new_dmap)
{
    if (ooo_debug) amrex::Print() << "MakeNewLevelFromScratch" << std::endl;
    if (m_verbose > 0)
    {
        std::cout << "MAKING NEW LEVEL " << lev << std::endl;
        std::cout << "WITH BOX ARRAY   " << new_grids << std::endl;
    }

    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    amrex::Print() << "SETTING NEW GRIDS IN MAKE NEW LEVEL " << new_grids << std::endl;
    amrex::Print() << "SETTING NEW DMAP IN MAKE NEW LEVEL " << new_dmap << std::endl;

    macproj = std::make_unique<MacProjector>(Geom(0,finest_level),
                                   MLMG::Location::FaceCentroid,  // Location of mac_vec
                                   MLMG::Location::FaceCentroid,  // Location of beta
                                   MLMG::Location::CellCenter,    // Location of solution variable phi
                                   MLMG::Location::CellCentroid);// Location of MAC RHS

    // This is being done by mfix::make_eb_geometry,
    // otherwise it would be done here
    if (lev == 0) MakeBCArrays(nghost_state());
}


void mfix::ReMakeNewLevelFromScratch (int lev,
                                      const BoxArray & new_grids,
                                      const DistributionMapping & new_dmap)
{
    if (ooo_debug) amrex::Print() << "ReMakeNewLevelFromScratch" << std::endl;
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

    if (lev == 0) MakeBCArrays(nghost_state());

    // We need to re-fill these arrays for the larger domain (after replication).
    mfix_set_bc_type(lev,nghost_state());
}



void mfix::InitLevelData (Real time)
{
    if (ooo_debug) amrex::Print() << "InitLevelData" << std::endl;
    // Allocate the fluid data, NOTE: this depends on the ebfactories.
    if (fluid.solve)
       for (int lev = 0; lev < nlev; lev++)
          AllocateArrays(lev);

    // Allocate the particle data
    if (DEM::solve || PIC::solve)
    {
      Real strt_init_part = ParallelDescriptor::second();

      pc->AllocData();

      if (particle_init_type == "AsciiFile")
      {
        amrex::Print() << "Reading particles from particle_input.dat ..." << std::endl;
        pc->InitParticlesAscii("particle_input.dat");

      } else if (particle_init_type == "Random")
      {
        int n_per_cell = 1;

        amrex::Print() << "Randomly initializing " << n_per_cell
                       << " particles per cell ..."
                       << std::endl;

        Real      radius = 1.0;
        Real      volume = 1.0;
        Real        mass = 1.0;
        Real     density = 1.0;
        Real        omoi = 1.0;
        Real        velx = 0.0;
        Real        vely = 0.0;
        Real        velz = 0.0;
        Real   dragcoeff = 0.0;
        Real       dragx = 0.0;
        Real       dragy = 0.0;
        Real       dragz = 0.0;
        Real      omegax = 0.0;
        Real      omegay = 0.0;
        Real      omegaz = 0.0;
        Real      statwt = 1.0;
        Real        cp_s = 1.0;
        Real temperature = 0.0;
        Real  convection = 0.0;
        int        phase = 1;
        int        state = 0;

        MFIXParticleContainer::ParticleInitData pdata =
          {{}, {}, {radius,volume,mass,density,omoi,velx,vely,velz,
            omegax,omegay,omegaz,dragx,dragy,dragz,dragcoeff,statwt,cp_s,
            temperature,convection}, {phase,state}};

        pc->InitNRandomPerCell(n_per_cell, pdata);
        pc->WriteAsciiFileForInit("random_particles");
        exit(0);

      } else if (particle_init_type == "Auto") {

         amrex::Print() << "Auto generating particles ..." << std::endl;

         pc->InitParticlesAuto();

      } else {

         amrex::Abort("Bad particle_init_type");
      }

      pc->Redistribute();

      Real end_init_part = ParallelDescriptor::second() - strt_init_part;
      ParallelDescriptor::ReduceRealMax(end_init_part, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "Time spent in initializing particles " << end_init_part << std::endl;
    }

    // Used in load balancing
    if (DEM::solve || PIC::solve)
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
    if (fluid.solve)
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
mfix::PostInit (Real& dt, Real time, int restart_flag, Real stop_time)
{
    if (ooo_debug) amrex::Print() << "PostInit" << std::endl;

    if (DEM::solve || PIC::solve)
    {
        // Auto generated particles may be out of the domain. This call will
        // remove them. Note that this has to occur after the EB geometry is
        // created. if (particle_init_type == "Auto" && !restart_flag &&
        // particle_ebfactory[finest_level])

      if (removeOutOfRange)
        {

          if ((nlev == 1) && (!restart_flag && particle_ebfactory[finest_level]))
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
          else if (!restart_flag && particle_ebfactory[finest_level])
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

        // We need to do this *after* restart (hence putting this here not
        // in Init) because we may want to change the particle_max_grid_size on restart.
        if ( dual_grid && particle_max_grid_size_x > 0
                       && particle_max_grid_size_y > 0
                       && particle_max_grid_size_z > 0)
        {
          IntVect particle_max_grid_size(particle_max_grid_size_x,
                                         particle_max_grid_size_y,
                                         particle_max_grid_size_z);

          for (int lev = 0; lev < nlev; lev++)
          {
            BoxArray particle_ba(geom[lev].Domain());
            particle_ba.maxSize(particle_max_grid_size);

            DistributionMapping particle_dm; 
            if (particle_pmap.empty())
              particle_dm.define(particle_ba, ParallelDescriptor::NProcs());
            else
              particle_dm.define(particle_pmap);

            pc->Regrid(particle_dm, particle_ba);

            if (particle_cost[lev] != nullptr)
              delete particle_cost[lev];

            particle_cost[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                              pc->ParticleDistributionMap(lev), 1, 0);
            particle_cost[lev]->setVal(0.0);

            // initialize the ranks of particle grids
            if (particle_proc[lev] != nullptr)
              delete particle_proc[lev];
            //
            const Real proc = Real(ParallelDescriptor::MyProc());
            particle_proc[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                              pc->ParticleDistributionMap(lev), 1, 0);
            particle_proc[lev]->setVal(proc);


            // This calls re-creates a proper particle_ebfactories
            //  and regrids all the multifabs that depend on it
            if (DEM::solve || PIC::solve)
                RegridLevelSetArray(lev);

          }
        }

        if (DEM::solve) {
            pc->MFIX_PC_InitCollisionParams();
        }

        pc->InitParticlesRuntimeVariables(advect_enthalpy, solids.solve_species);

        if (!fluid.solve){
            dt = fixed_dt;
        }
    }

    if (fluid.solve)
        mfix_init_fluid(restart_flag, dt, stop_time);

    // Call user-defined subroutine to set constants, check data, etc.
    if (call_udf) mfix_usr0();
}

void
mfix::MakeBCArrays (int nghost)
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

void
mfix::mfix_init_fluid (int is_restarting, Real dt, Real stop_time)
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

          if ( is_restarting ) {
            init_fluid_parameters(bx, mfi, ld, advect_enthalpy, advect_fluid_species, fluid);
          } else {
            init_fluid(sbx, bx, domain, mfi, ld, dx, dy, dz, xlen, ylen, zlen, plo,
                test_tracer_conservation, advect_enthalpy, advect_fluid_species,
                m_idealgas_constraint, fluid);
          }
       }

    }

    mfix_set_p0();

    // Here we re-set the bc values for p and u,v,w just in case init_fluid
    //      over-wrote some of the bc values with ic values
    mfix_set_bc0();

    for (int lev = 0; lev < nlev; lev++)
    {
      m_leveldata[lev]->ep_g->FillBoundary(geom[lev].periodicity());
      m_leveldata[lev]->ro_g->FillBoundary(geom[lev].periodicity());

      if (advect_enthalpy)
      {
        m_leveldata[lev]->h_g->FillBoundary(geom[lev].periodicity());
        m_leveldata[lev]->T_g->FillBoundary(geom[lev].periodicity());
      }

      if (advect_tracer)
        m_leveldata[lev]->trac->FillBoundary(geom[lev].periodicity());

      if (advect_fluid_species)
      {
        m_leveldata[lev]->X_gk->FillBoundary(geom[lev].periodicity());
      }

      m_leveldata[lev]->vel_g->FillBoundary(geom[lev].periodicity());
    }


    // Make sure to fill the "old state" before we start.
    for (int lev = 0; lev < nlev; lev++)
    {
       LevelData& ld = *m_leveldata[lev];

       // TODO check if commenting the following is correct
       //MultiFab::Copy(*m_leveldata[lev]->vel_go,  *m_leveldata[lev]->vel_g, 0, 0, 3, 0);

       MultiFab::Copy(*ld.ro_go,  *ld.ro_g, 0, 0, 1, ld.ro_g->nGrow());
       MultiFab::Copy(*ld.trac_o, *ld.trac, 0, 0, 1, ld.trac->nGrow());

       if (advect_enthalpy) {
         MultiFab::Copy(*ld.T_go, *ld.T_g, 0, 0, 1, ld.T_g->nGrow());
         MultiFab::Copy(*ld.h_go, *ld.h_g, 0, 0, 1, ld.h_g->nGrow());
       }

       if (advect_fluid_species) {
         MultiFab::Copy(*ld.X_gko, *ld.X_gk, 0, 0, fluid.nspecies, ld.X_gk->nGrow());
       }

       if (m_idealgas_constraint == IdealGasConstraint::ClosedSystem) {
         MultiFab::Copy(*ld.pressure_go, *ld.pressure_g, 0, 0, 1, ld.pressure_g->nGrow());
       }
    }


    if (is_restarting == 0)
    {
      // Just for reference, we compute the volume inside the EB walls (as if
      // there were no particles)
      m_leveldata[0]->ep_g->setVal(1.0);

      const Real* dx = geom[0].CellSize();
      const Real cell_volume = dx[0] * dx[1] * dx[2];

      sum_vol_orig = volWgtSum(0,*(m_leveldata[0]->ep_g),0);

      Print() << "Enclosed domain volume is   " << cell_volume * sum_vol_orig << std::endl;

      Real domain_vol = sum_vol_orig;

      // Now initialize the volume fraction ep_g before the first projection
      mfix_calc_volume_fraction(sum_vol_orig);
      Print() << "Setting original sum_vol to " << cell_volume * sum_vol_orig << std::endl;

      Print() << "Difference is   " << cell_volume * (domain_vol - sum_vol_orig) << std::endl;

      // This sets bcs for ep_g
      Real time = 0.0;

      mfix_set_density_bcs(time, get_ro_g());
      mfix_set_density_bcs(time, get_ro_g_old());

      mfix_set_tracer_bcs(time, get_trac());
      mfix_set_tracer_bcs(time, get_trac_old());

      if (advect_enthalpy) {
        mfix_set_temperature_bcs(time, get_T_g());
        mfix_set_temperature_bcs(time, get_T_g_old());
      }

      if (advect_enthalpy) {
        mfix_set_enthalpy_bcs(time, get_h_g());
        mfix_set_enthalpy_bcs(time, get_h_g_old());
      }

      if (advect_fluid_species) {
        mfix_set_species_bcs(time, get_X_gk());
        mfix_set_species_bcs(time, get_X_gk_old());
      }

      // Project the initial velocity field
      if (do_initial_proj)
        mfix_project_velocity();

      // Iterate to compute the initial pressure
      if (initial_iterations > 0)
        mfix_initial_iterations(dt,stop_time);
    }
    else
    {
      const int dir_bc = 1;
      mfix_set_epg_bcs(get_ep_g(), dir_bc);

      const Real* dx = geom[0].CellSize();
      const Real cell_volume = dx[0] * dx[1] * dx[2];

      //Calculation of sum_vol_orig for a restarting point
      sum_vol_orig = volWgtSum(0,*(m_leveldata[0]->ep_g),0);

      Print() << "Setting original sum_vol to " << cell_volume * sum_vol_orig << std::endl;
    }
}

void
mfix::mfix_set_bc0 ()
{
    if (ooo_debug) amrex::Print() << "mfix_set_bc0" << std::endl;
    for (int lev = 0; lev < nlev; lev++)
    {
     Box domain(geom[lev].Domain());

     MultiFab& ep_g = *(m_leveldata[lev]->ep_g);

     // Don't tile this -- at least for now
     for (MFIter mfi(ep_g, false); mfi.isValid(); ++mfi)
     {
       const Box& sbx = ep_g[mfi].box();

       if (advect_enthalpy)
         set_temperature_bc0(sbx, &mfi, lev, domain);

       if (advect_fluid_species)
         set_species_bc0(sbx, &mfi, lev, domain);
     }

     m_leveldata[lev]->ep_g->FillBoundary(geom[lev].periodicity());
     m_leveldata[lev]->ro_g->FillBoundary(geom[lev].periodicity());

     if (advect_enthalpy) {
       m_leveldata[lev]->h_g->FillBoundary(geom[lev].periodicity());
       m_leveldata[lev]->T_g->FillBoundary(geom[lev].periodicity());
     }

     if (advect_tracer)
       m_leveldata[lev]->trac->FillBoundary(geom[lev].periodicity());

     if (advect_fluid_species)
       m_leveldata[lev]->X_gk->FillBoundary(geom[lev].periodicity());
   }

   // Put velocity Dirichlet bc's on faces
   Real time = 0.0;
   int extrap_dir_bcs = 0;

   mfix_set_velocity_bcs(time, get_vel_g(), extrap_dir_bcs);

   for (int lev = 0; lev < nlev; lev++)
     m_leveldata[lev]->vel_g->FillBoundary(geom[lev].periodicity());
}

void
mfix::mfix_set_p0 ()
{
  if (ooo_debug) amrex::Print() << "mfix_set_p0" << std::endl;

  IntVect press_per = IntVect(geom[0].isPeriodic(0),geom[0].isPeriodic(1),geom[0].isPeriodic(2));

  // Here we set a separate periodicity flag for p0_g because when we use
  // pressure drop (delp) boundary conditions we fill all variables *except* p0
  // periodically
  if (BC::delp_dir > -1) press_per[BC::delp_dir] = 0;
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

     m_leveldata[lev]->p0_g->FillBoundary(p0_periodicity);
   }
}

void mfix::mfix_set_ls_near_inflow ()
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
