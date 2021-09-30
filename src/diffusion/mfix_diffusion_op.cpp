#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB_utils.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <mfix_diffusion_op.H>
#include <mfix_eb_parms.H>
#include <mfix_bc_parms.H>
#include <mfix_species_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_mf_helpers.H>

using namespace amrex;

//
// Constructor:
// We set up everything which doesn't change between timesteps here
//
DiffusionOp::DiffusionOp (AmrCore* _amrcore,
                          Vector< const EBFArrayBoxFactory* >const& _ebfactory,
                          FluidPhase& _fluid,
                          std::array<amrex::LinOpBCType,3> a_velbc_lo,
                          std::array<amrex::LinOpBCType,3> a_velbc_hi,
                          std::array<amrex::LinOpBCType,3> a_scalbc_lo,
                          std::array<amrex::LinOpBCType,3> a_scalbc_hi,
                          std::array<amrex::LinOpBCType,3> a_temperaturebc_lo,
                          std::array<amrex::LinOpBCType,3> a_temperaturebc_hi,
                          std::array<amrex::LinOpBCType,3> a_speciesbc_lo,
                          std::array<amrex::LinOpBCType,3> a_speciesbc_hi,
                          int _nghost)
  : nghost(_nghost)
  , m_velbc_lo(a_velbc_lo)
  , m_velbc_hi(a_velbc_hi)
  , m_scalbc_lo(a_scalbc_lo)
  , m_scalbc_hi(a_scalbc_hi)
  , m_temperaturebc_lo(a_temperaturebc_lo)
  , m_temperaturebc_hi(a_temperaturebc_hi)
  , m_speciesbc_lo(a_speciesbc_lo)
  , m_speciesbc_hi(a_speciesbc_hi)
  , eb_is_dirichlet(false) // We default to Neumann bc's for scalarson EB walls
  , fluid(_fluid)
{
    if(verbose > 0)
        amrex::Print() << "Constructing DiffusionOp class" << std::endl;


    // Get inputs from ParmParse
    readParameters();

    // Actually do the setup work here
    setup(_amrcore, _ebfactory);
}

void DiffusionOp::setup (AmrCore* _amrcore,
                         Vector< const EBFArrayBoxFactory* >const& _ebfactory)
{
    // The amrcore boxArray and DistributionMap change when we regrid so we must
    // pass the new object in here.
    amrcore = _amrcore;

    // The ebfactory changes when we regrid so we must pass it in here.
    ebfactory = _ebfactory;

    geom  = amrcore->Geom();
    grids = amrcore->boxArray();
    dmap  = amrcore->DistributionMap();

    int max_level = amrcore->maxLevel();

    const int nspecies_g = fluid.nspecies;

    // Resize and reset data
    b.resize(max_level + 1);

    phi.resize(max_level + 1);
    rhs.resize(max_level + 1);
    vel_eb.resize(max_level + 1);

    if (fluid.solve_species)
    {
      species_phi.resize(max_level + 1);
      species_rhs.resize(max_level + 1);
    }

    for(int lev = 0; lev <= max_level; lev++)
    {
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            BoxArray edge_ba = grids[lev];
            edge_ba.surroundingNodes(dir);
            b[lev][dir] = std::make_unique<MultiFab>(edge_ba, dmap[lev], 1, nghost,
                                           MFInfo(), *ebfactory[lev]);
            b[lev][dir]->setVal(0.);
        }

        phi[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], 3, 1,
                                    MFInfo(), *ebfactory[lev]);
        phi[lev]->setVal(0.);

        // No ghost cells needed for rhs
        rhs[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], 3, 0,
                                    MFInfo(), *ebfactory[lev]);
        rhs[lev]->setVal(0.);

        vel_eb[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], 3, nghost,
                                       MFInfo(), *ebfactory[lev]);
        vel_eb[lev]->setVal(0.0);

        if (fluid.solve_species)
        {
          species_phi[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], nspecies_g, 1,
                                              MFInfo(), *ebfactory[lev]);
          species_phi[lev]->setVal(0.);

          // No ghost cells needed for rhs
          species_rhs[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], nspecies_g, 0,
                                              MFInfo(), *ebfactory[lev]);
          species_rhs[lev]->setVal(0.);
        }
    }

    //
    // Define the matrix for the viscous tensor solve.
    //
    LPInfo info;
    info.setMaxCoarseningLevel(mg_max_coarsening_level);
    info.setAgglomerationGridSize(mg_agg_grid_size);
    vel_matrix = std::make_unique<MLEBTensorOp>(geom, grids, dmap, info, ebfactory);

    // It is essential that we set MaxOrder to 2 if we want to use the standard
    // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the gradient.
    vel_matrix->setMaxOrder(2);

    // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
    vel_matrix->setDomainBC(m_velbc_lo, m_velbc_hi);

    //
    // Define the matrix for the scalar diffusion solve.
    //
    scal_matrix = std::make_unique<MLEBABecLap>(geom, grids, dmap, info, ebfactory);
    temperature_matrix = std::make_unique<MLEBABecLap>(geom, grids, dmap, info, ebfactory);

    if (fluid.solve_species) {
      species_matrix = std::make_unique<MLEBABecLap>(geom, grids, dmap, info, ebfactory, nspecies_g);
    }

    // It is essential that we set MaxOrder to 2 if we want to use the standard
    // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the gradient.
    scal_matrix->setMaxOrder(2);
    temperature_matrix->setMaxOrder(2);

    if (fluid.solve_species) {
      species_matrix->setMaxOrder(2);
    }

    // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
    scal_matrix->setDomainBC(m_scalbc_lo, m_scalbc_hi);
    temperature_matrix->setDomainBC(m_temperaturebc_lo, m_temperaturebc_hi);

    if (fluid.solve_species)
    {
      species_matrix->setDomainBC(m_speciesbc_lo, m_speciesbc_hi);

      species_b.resize(max_level + 1);

      for(int lev = 0; lev <= max_level; lev++)
      {
          for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
          {
              BoxArray edge_ba = grids[lev];
              edge_ba.surroundingNodes(dir);
              species_b[lev][dir] = std::make_unique<MultiFab>(edge_ba, dmap[lev], nspecies_g, nghost,
                                                     MFInfo(), *ebfactory[lev]);
              species_b[lev][dir]->setVal(0.);
          }
      }
    }

}

DiffusionOp::~DiffusionOp ()
{}

void DiffusionOp::readParameters ()
{
    ParmParse pp("diffusion");

    pp.query("verbose_solver", verbose);
    pp.query("verbose", mg_verbose);
    pp.query("bottom_verbose", mg_bottom_verbose);
    pp.query("maxiter", mg_maxiter);
    pp.query("bottom_maxiter", mg_bottom_maxiter);
    pp.query("mg_max_fmg_iter", mg_max_fmg_iter);
    pp.query("mg_max_coarsening_level", mg_max_coarsening_level);
    pp.query("rtol", mg_rtol);
    pp.query("atol", mg_atol);
    pp.query("bottom_solver", bottom_solver_type);
    pp.query("agg_grid_size", mg_agg_grid_size);
}


//
// Set the user-supplied settings for the MLMG solver
// (this must be done every time step, since MLMG is created after updating matrix
//
void DiffusionOp::setSolverSettings (MLMG& solver)
{
    // The default bottom solver is BiCG
    if(bottom_solver_type == "smoother")
    {
        solver.setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if(bottom_solver_type == "hypre")
    {
        solver.setBottomSolver(MLMG::BottomSolver::hypre);
    }
        // Maximum iterations for MultiGrid / ConjugateGradients
        solver.setMaxIter(mg_maxiter);
        solver.setMaxFmgIter(mg_max_fmg_iter);
        solver.setBottomMaxIter(mg_bottom_maxiter);

        // Verbosity for MultiGrid / ConjugateGradients
        solver.setVerbose(mg_verbose);
        solver.setBottomVerbose(mg_bottom_verbose);

        // This ensures that ghost cells of phi are correctly filled when
        // returned from the solver
        solver.setFinalFillBC(true);
}

void DiffusionOp::ComputeDivTau (const Vector< MultiFab* >& divtau_out,
                                 const Vector< MultiFab* >& vel_in,
                                 const Vector< MultiFab* >& ep_in,
                                 const Vector< MultiFab* >& T_g_in,
                                 const int advect_enthalpy)
{
    BL_PROFILE("DiffusionOp::ComputeDivTau");

    int finest_level = amrcore->finestLevel();

    Vector< MultiFab* > divtau_aux(finest_level+1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
       divtau_aux[lev] = new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                      MFInfo(), *ebfactory[lev]);
       divtau_aux[lev]->setVal(0.0);
    }

    // We want to return div (mu grad)) phi
    vel_matrix->setScalars(0.0, -1.0);

    Vector<BCRec> bcs_s; // This is just to satisfy the call to EB_interp...
    bcs_s.resize(3);
    // Compute the coefficients
    for (int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab mu_g(ep_in[lev]->boxArray(), ep_in[lev]->DistributionMap(),
                      ep_in[lev]->nComp(), ep_in[lev]->nGrow(), MFInfo(),
                      ep_in[lev]->Factory());

        mu_g.setVal(0);

        const Real mu_g0 = fluid.mu_g0;

        auto& fluid_parms = *fluid.parameters;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.growntilebox(vel_in[lev]->nGrowVect());

          if (bx.ok())
          {
            Array4<Real      > const& mu_g_array = mu_g.array(mfi);
            Array4<Real const> const& T_g_array  = advect_enthalpy ?
              T_g_in[lev]->const_array(mfi) : Array4<const Real>();

            ParallelFor(bx, [mu_g_array,T_g_array,advect_enthalpy,mu_g0,fluid_parms]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              if (advect_enthalpy)
                mu_g_array(i,j,k) = fluid_parms.calc_mu_g(T_g_array(i,j,k));
              else
                mu_g_array(i,j,k) = mu_g0;
            });
          }
        }

        mu_g.FillBoundary(geom[lev].periodicity());

        //EB_set_covered(mu_g, 0, mu_g.nComp(), mu_g.nGrow(), covered_val);
        EB_set_covered(mu_g, 0, mu_g.nComp(), mu_g.nGrow(), 1.e40);

        // average_cellcenter_to_face( GetArrOfPtrs(b[lev]), *eta_in[lev], geom[lev] );
        EB_interp_CellCentroid_to_FaceCentroid (mu_g, GetArrOfPtrs(b[lev]), 0, 0, 1, geom[lev], bcs_s);

        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
             b[lev][dir]->FillBoundary(geom[lev].periodicity());

        vel_matrix->setShearViscosity(lev, GetArrOfConstPtrs(b[lev]), MLMG::Location::FaceCentroid);
        vel_matrix->setEBShearViscosity(lev, mu_g);
        vel_matrix->setLevelBC(lev, GetVecOfConstPtrs(vel_in)[lev]);
    }

    MLMG solver(*vel_matrix);

    solver.apply(divtau_aux, vel_in);

    for(int lev = 0; lev <= finest_level; lev++)
    {
       amrex::single_level_weighted_redistribute(*divtau_aux[lev],
                                                 *divtau_out[lev],
                                                 *ep_in[lev],
                                                 0,
                                                 AMREX_SPACEDIM,
                                                 geom[lev]);

       EB_set_covered(*divtau_out[lev], 0, divtau_out[lev]->nComp(), divtau_out[lev]->nGrow(), 0.);
    }

    for(int lev = 0; lev <= finest_level; lev++)
      delete divtau_aux[lev];
}

void DiffusionOp::ComputeLapT (const Vector< MultiFab*      >& lapT_out,
                               const Vector< MultiFab*      >& T_g,
                               const Vector< MultiFab const*>& ep_g,
                               const Vector< MultiFab const*>& T_g_on_eb)
{
  BL_PROFILE("DiffusionOp::ComputeLapT");

  int finest_level = amrcore->finestLevel();

  Vector< MultiFab* > lapT_aux(finest_level+1);

  for(int lev = 0; lev <= finest_level; lev++)
  {
    lapT_aux[lev] = new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(),
        *ebfactory[lev]);

    lapT_aux[lev]->setVal(0.0);
  }

  // We want to return div (ep_g k_g grad)) T_g
  temperature_matrix->setScalars(0.0, -1.0);

  Vector<BCRec> bcs_s; // This is just to satisfy the call to EBterp...
  bcs_s.resize(3);

  // Compute the coefficients
  for (int lev = 0; lev <= finest_level; lev++) {

    MultiFab ep_k_g(ep_g[lev]->boxArray(), ep_g[lev]->DistributionMap(),
                    ep_g[lev]->nComp(), ep_g[lev]->nGrow(), MFInfo(),
                    ep_g[lev]->Factory());

    // Initialize to 0
    ep_k_g.setVal(0.);

    auto& fluid_parms = *fluid.parameters;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.growntilebox(IntVect(1,1,1));

      if (bx.ok())
      {
        Array4<Real      > const& ep_k_g_array = ep_k_g.array(mfi);
        Array4<Real const> const& ep_g_array  = ep_g[lev]->const_array(mfi);
        Array4<Real const> const& T_g_array   = T_g[lev]->const_array(mfi);

        amrex::ParallelFor(bx, [ep_g_array,T_g_array,ep_k_g_array,fluid_parms]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          ep_k_g_array(i,j,k) = ep_g_array(i,j,k)*fluid_parms.calc_k_g(T_g_array(i,j,k));
        });
      }
    }

//    ep_k_g.FillBoundary(geom[lev].periodicity());

    EB_interp_CellCentroid_to_FaceCentroid (ep_k_g, GetArrOfPtrs(b[lev]), 0, 0, 1, geom[lev], bcs_s);

//    b[lev][0]->FillBoundary(geom[lev].periodicity());
//    b[lev][1]->FillBoundary(geom[lev].periodicity());
//    b[lev][2]->FillBoundary(geom[lev].periodicity());

    if (EB::fix_temperature) {
      // The following is a WIP in AMReX
      //temperature_matrix->setPhiOnCentroid();

      MultiFab k_g_on_eb(T_g_on_eb[lev]->boxArray(), T_g_on_eb[lev]->DistributionMap(),
                         T_g_on_eb[lev]->nComp(), T_g_on_eb[lev]->nGrow(), MFInfo(),
                         T_g_on_eb[lev]->Factory());

      k_g_on_eb.setVal(0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*T_g_on_eb[lev]); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.growntilebox(T_g_on_eb[lev]->nGrowVect());

        if (bx.ok()) {
          Array4<Real      > const& k_g_on_eb_array = k_g_on_eb.array(mfi);
          Array4<Real const> const& T_g_on_eb_array = T_g_on_eb[lev]->const_array(mfi);

          amrex::ParallelFor(bx, [k_g_on_eb_array,T_g_on_eb_array,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            if (T_g_on_eb_array(i,j,k) > 0)
              k_g_on_eb_array(i,j,k) = fluid_parms.calc_k_g(T_g_on_eb_array(i,j,k));
          });
        }
      }

      k_g_on_eb.FillBoundary(geom[lev].periodicity());

      temperature_matrix->setEBDirichlet(lev, *T_g_on_eb[lev], k_g_on_eb);
    }

    temperature_matrix->setBCoeffs(lev, GetArrOfConstPtrs(b[lev]),
        MLMG::Location::FaceCentroid);

    temperature_matrix->setLevelBC(lev, GetVecOfConstPtrs(T_g)[lev]);
  }

  MLMG solver(*temperature_matrix);

  solver.apply(lapT_aux, T_g);

  for(int lev = 0; lev <= finest_level; lev++) {
    amrex::single_level_redistribute(*lapT_aux[lev], *lapT_out[lev], 0, 1, geom[lev]);
    EB_set_covered(*lapT_out[lev], 0, lapT_out[lev]->nComp(), lapT_out[lev]->nGrow(), 0.);
  }

  for(int lev = 0; lev <= finest_level; lev++)
  {
    delete lapT_aux[lev];
  }
}

void DiffusionOp::ComputeLapS (const Vector< MultiFab* >& laps_out,
                               const Vector< MultiFab* >& scal_in,
                               const Vector< MultiFab const*>& /*ro_in*/,
                               const Vector< MultiFab const*>& /*ep_in*/,
                               const Vector< Real >& mu_s)
{
    BL_PROFILE("DiffusionOp::ComputeLapS");

    int finest_level = amrcore->finestLevel();

    int ntrac = scal_in[0]->nComp();

    Vector< MultiFab* >  laps_aux(finest_level+1);
    Vector< MultiFab* >    phi_eb(finest_level+1);
    for(int lev = 0; lev <= finest_level; lev++) {

      laps_aux[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, nghost,
                                    MFInfo(), *ebfactory[lev]);
       laps_aux[lev]->setVal(0.0);

       phi_eb[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, 0,
                                  MFInfo(), *ebfactory[lev]);
       phi_eb[lev]->setVal(0.);

       // This value was just for testing
       // if (eb_is_dirichlet)
       //    phi_eb[lev]->setVal(1.0);
    }

    // We want to return div (mu grad)) phi
    scal_matrix->setScalars(0.0, -1.0);

    // Compute the coefficients
    for (int lev = 0; lev <= finest_level; lev++) {

      for(int dir = 0; dir < 3; dir++)
           for(int n = 0; n < ntrac; n++)
             b[lev][dir]->setVal(mu_s[n],n,1);

        if (eb_is_dirichlet)
            scal_matrix->setEBDirichlet(lev, *phi_eb[lev], mu_s);

        scal_matrix->setBCoeffs(lev, GetArrOfConstPtrs(b[lev]),MLMG::Location::FaceCentroid);
        scal_matrix->setLevelBC(lev, GetVecOfConstPtrs(scal_in)[lev]);
    }

    MLMG solver(*scal_matrix);

    solver.apply(laps_aux, scal_in);

    for(int lev = 0; lev <= finest_level; lev++) {
      amrex::single_level_redistribute(*laps_aux[lev], *laps_out[lev], 0, ntrac, geom[lev]);
      EB_set_covered(*laps_out[lev], 0, laps_out[lev]->nComp(), laps_out[lev]->nGrow(), 0.);
    }

    for(int lev = 0; lev <= finest_level; lev++)
    {
       delete laps_aux[lev];
       delete   phi_eb[lev];
    }
}

void DiffusionOp::ComputeLapX (const Vector< MultiFab*      >& lapX_out,
                               const Vector< MultiFab*      >& X_gk_in,
                               const Vector< MultiFab const*>& ro_g_in,
                               const Vector< MultiFab const*>& ep_g_in,
                               const Vector< MultiFab const*>& T_g_in)
{
  BL_PROFILE("DiffusionOp::ComputeLapX");

  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

//  // TODO: check on this
//  const bool already_on_centroids = true;

  int finest_level = amrcore->finestLevel();

  // Number of fluid species
  const int nspecies_g = fluid.nspecies;

  Vector<BCRec> bcs_dummy; // This is just to satisfy the call to EB_interp...
  bcs_dummy.resize(3*nspecies_g);

  // Auxiliary data where we store Div{ep_g ro_g D_gk Grad{X_gk}}
  Vector< MultiFab* > lapX_aux(finest_level+1);

  // Allocate space for lapX_aux and set it to 0
  for(int lev = 0; lev <= finest_level; lev++)
  {
    lapX_aux[lev] = new MultiFab(grids[lev], dmap[lev], nspecies_g, nghost, MFInfo(),
        *ebfactory[lev]);

    lapX_aux[lev]->setVal(0.0);
  }

  // We want to return div (ep_g ro_g D_gk grad)) phi
  species_matrix->setScalars(0.0, -1.0);

  // Compute the coefficients
  for (int lev = 0; lev <= finest_level; lev++)
  {
    MultiFab b_coeffs(ep_g_in[lev]->boxArray(), ep_g_in[lev]->DistributionMap(),
        nspecies_g, 1, MFInfo(), ep_g_in[lev]->Factory());

    b_coeffs.setVal(0.);

    auto& fluid_parms = *fluid.parameters;

    // b_coeffs  = ep_g ro_g D_gk
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.growntilebox(IntVect(1,1,1));

      Array4<Real      > const& b_coeffs_arr = b_coeffs.array(mfi);
      Array4<Real const> const& ep_g_arr     = ep_g_in[lev]->const_array(mfi);
      Array4<Real const> const& ro_g_arr     = ro_g_in[lev]->const_array(mfi);
      Array4<Real const> const& T_g_arr      = T_g_in[lev]->const_array(mfi);

      amrex::ParallelFor(bx, [ep_g_arr,ro_g_arr,T_g_arr,b_coeffs_arr,nspecies_g,
          fluid_parms,run_on_device]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const Real ep_g = ep_g_arr(i,j,k);
        const Real ro_g = ro_g_arr(i,j,k);
        const Real T_g  = T_g_arr(i,j,k);

        for (int n(0); n < nspecies_g; ++n) {
          const Real D_gk = run_on_device ?
            fluid_parms.calc_D_gk<RunOn::Device>(T_g,n) :
            fluid_parms.calc_D_gk<RunOn::Host>(T_g,n);

          b_coeffs_arr(i,j,k,n) = ep_g*ro_g*D_gk;
        }
      });
    }

    // species_b = interp(b_coeffs)
    EB_interp_CellCentroid_to_FaceCentroid (b_coeffs, GetArrOfPtrs(species_b[lev]), 0, 0,
                                            nspecies_g, geom[lev], bcs_dummy);

    // Set BCoeffs
    species_matrix->setBCoeffs(lev, GetArrOfConstPtrs(species_b[lev]), MLMG::Location::FaceCentroid);

    // Set LevelBC
    species_matrix->setLevelBC(lev, GetVecOfConstPtrs(X_gk_in)[lev]);
  }

  MLMG solver(*species_matrix);
//  setSolverSettings(solver);
//
//  // This ensures that ghost cells of sol are correctly filled when returned
//  // from the solver
//  solver.setFinalFillBC(true);

  // Compute div (ep_g ro_g D_gk grad)) phi
  solver.apply(lapX_aux, X_gk_in);

#ifdef AMREX_DEBUG
  {
    for (int lev(0); lev <= finest_level; ++lev) {
      MultiFab temp(grids[lev], dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
      temp.setVal(0.);

      for(int n(0); n < nspecies_g; ++n)
        MultiFab::Add(temp, *lapX_aux[lev], n, 0, 1, 0);

      Print() << "lev = " << lev << std::endl;
      Print() << "summed div fluxes max = " << temp.max(0, 0) << std::endl;
      Print() << "summed div fluxes min = " << temp.min(0, 0) << std::endl;
    }
  }
#endif

  // CORRECT FLUXES
  // Compute Fluxes for correcting the result
  // div(ep_g ro_g D_gk grad(phi)) - div(phi sum(ep_g ro_g D_gk grad(phi)))
  for (int species_k(0); species_k < nspecies_g; ++species_k) {

    Vector<BCRec> loc_bcs_dummy; // This is just to satisfy the call to EB_interp...
    loc_bcs_dummy.resize(3*nspecies_g);

    // Auxiliary data where we store Div{ep_g ro_g D_gk Grad{X_gk}}
    Vector< MultiFab* > correction_aux(finest_level+1);

    // Allocate space for lapX_aux and set it to 0
    for(int lev = 0; lev <= finest_level; lev++)
    {
      correction_aux[lev] = new MultiFab(grids[lev], dmap[lev], nspecies_g, nghost, MFInfo(),
                                         *ebfactory[lev]);

      correction_aux[lev]->setVal(0.0);
    }

    // We want to return div (ep_g ro_g D_gk grad)) phi
    species_matrix->setScalars(0.0, -1.0);

    // Compute the coefficients
    for (int lev = 0; lev <= finest_level; lev++)
    {
      MultiFab Xb_coeffs(ep_g_in[lev]->boxArray(), ep_g_in[lev]->DistributionMap(),
                         nspecies_g, 1, MFInfo(), ep_g_in[lev]->Factory());

      Xb_coeffs.setVal(0.);

      auto& fluid_parms = *fluid.parameters;

      // Xb_coeffs  = ep_g ro_g X_gk D_gm
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*ep_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.growntilebox(IntVect(1,1,1));

        Array4<Real      > const& Xb_coeffs_arr = Xb_coeffs.array(mfi);
        Array4<Real const> const& ep_g_arr      = ep_g_in[lev]->const_array(mfi);
        Array4<Real const> const& ro_g_arr      = ro_g_in[lev]->const_array(mfi);
        Array4<Real const> const& T_g_arr       = T_g_in[lev]->const_array(mfi);
        Array4<Real const> const& X_gk_arr      = X_gk_in[lev]->const_array(mfi);

        amrex::ParallelFor(bx, [ep_g_arr,ro_g_arr,T_g_arr,Xb_coeffs_arr,
            X_gk_arr,nspecies_g,fluid_parms,species_k]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const Real ep_g = ep_g_arr(i,j,k);
          const Real ro_g = ro_g_arr(i,j,k);
          const Real T_g  = T_g_arr(i,j,k);
          const Real X_gk = X_gk_arr(i,j,k,species_k);

          const Real val = ep_g*ro_g*X_gk;

          for (int m(0); m < nspecies_g; ++m) {
            Xb_coeffs_arr(i,j,k,m) = val*fluid_parms.calc_D_gk<RunOn::Gpu>(T_g,m);
          }
        });
      }

      for(int dir = 0; dir < 3; dir++) {
        species_b[lev][dir]->setVal(0);
      }

      // species_b = interp(b_coeffs)
      EB_interp_CellCentroid_to_FaceCentroid (Xb_coeffs, GetArrOfPtrs(species_b[lev]), 0, 0,
                                              nspecies_g, geom[lev], loc_bcs_dummy);

      // Set BCoeffs
      species_matrix->setBCoeffs(lev, GetArrOfConstPtrs(species_b[lev]), MLMG::Location::FaceCentroid);

      // Set LevelBC
      species_matrix->setLevelBC(lev, GetVecOfConstPtrs(X_gk_in)[lev]);
    }

    MLMG solver(*species_matrix);

    // Compute div (ep_g ro_g D_gk grad)) phi
    solver.apply(correction_aux, X_gk_in);

    for (int lev(0); lev <= finest_level; ++lev) {

      for (int m(0); m < nspecies_g; ++m) {
        MultiFab::Subtract(*lapX_aux[lev], *correction_aux[lev], m, species_k, 1, lapX_aux[lev]->nGrow());
      }
    }

    for (int lev(0); lev <= finest_level; ++lev) {
      delete correction_aux[lev];
    }

  } // correct_fluxes

#ifdef AMREX_DEBUG
  {
    for (int lev(0); lev <= finest_level; ++lev) {
      MultiFab temp(grids[lev], dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
      temp.setVal(0.);

      for(int n(0); n < nspecies_g; ++n)
        MultiFab::Add(temp, *lapX_aux[lev], n, 0, 1, 0);

      Print() << "lev = " << lev << std::endl;
      Print() << "summed div fluxes max = " << temp.max(0, 0) << std::endl;
      Print() << "summed div fluxes min = " << temp.min(0, 0) << std::endl;
    }
  }
#endif

  // Redistribute lapX_aux into lapX_out
  for(int lev = 0; lev <= finest_level; lev++)
  {
    amrex::single_level_redistribute(*lapX_aux[lev], *lapX_out[lev], 0, nspecies_g, geom[lev]);
    EB_set_covered(*lapX_out[lev], 0, lapX_out[lev]->nComp(), lapX_out[lev]->nGrow(), 0.);
  }

  // Free lapX_aux memory
  for(int lev = 0; lev <= finest_level; lev++)
  {
    delete lapX_aux[lev];
  }
}


void DiffusionOp::SubtractDiv_XGradX (const Vector< MultiFab*      >& X_gk_in,
                                      const Vector< MultiFab const*>& ro_g_in,
                                      const Vector< MultiFab const*>& ep_g_in,
                                      const Vector< MultiFab const*>& T_g_in,
                                      const Real& dt)
{
  BL_PROFILE("DiffusionOp::SubtractDiv_XGradX");

  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

  // TODO: check on this
  const bool already_on_centroids = true;

  int finest_level = amrcore->finestLevel();

  // Number of fluid species
  const int nspecies_g = fluid.nspecies;

  Vector<BCRec> bcs_dummy; // This is just to satisfy the call to EB_interp...
  bcs_dummy.resize(3*nspecies_g);

  // Weaset it up for Div{rho_g D_gk Grad{X_gk}}
  species_matrix->setScalars(0.0, -1.0);

  // Compute the coefficients
  for (int lev = 0; lev <= finest_level; lev++)
  {
    MultiFab b_coeffs(ep_g_in[lev]->boxArray(), ep_g_in[lev]->DistributionMap(),
        nspecies_g, 1, MFInfo(), ep_g_in[lev]->Factory());

    // FROM HERE
    b_coeffs.setVal(0.);

    auto& fluid_parms = *fluid.parameters;

    // b_coeffs  = ep_g ro_g D_gk
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.growntilebox(IntVect(1,1,1));

      Array4<Real      > const& b_coeffs_arr = b_coeffs.array(mfi);
      Array4<Real const> const& ep_g_arr     = ep_g_in[lev]->const_array(mfi);
      Array4<Real const> const& ro_g_arr     = ro_g_in[lev]->const_array(mfi);
      Array4<Real const> const& T_g_arr      = T_g_in[lev]->const_array(mfi);

      amrex::ParallelFor(bx, [ep_g_arr,ro_g_arr,T_g_arr,b_coeffs_arr,nspecies_g,
          fluid_parms,run_on_device]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const Real ep_g = ep_g_arr(i,j,k);
        const Real ro_g = ro_g_arr(i,j,k);
        const Real T_g  = T_g_arr(i,j,k);

        for (int n(0); n < nspecies_g; ++n) {
          const Real D_gk = run_on_device ?
            fluid_parms.calc_D_gk<RunOn::Device>(T_g,n) :
            fluid_parms.calc_D_gk<RunOn::Host>(T_g,n);

          b_coeffs_arr(i,j,k,n) = ep_g*ro_g*D_gk;
        }
      });
    }

    // species_b = interp(b_coeffs)
    EB_interp_CellCentroid_to_FaceCentroid (b_coeffs, GetArrOfPtrs(species_b[lev]), 0, 0,
                                            nspecies_g, geom[lev], bcs_dummy);

    // Set BCoeffs
    species_matrix->setBCoeffs(lev, GetArrOfConstPtrs(species_b[lev]), MLMG::Location::FaceCentroid);

    // Set LevelBC
    species_matrix->setLevelBC(lev, GetVecOfConstPtrs(X_gk_in)[lev]);
  }

  MLMG solver(*species_matrix);
//  setSolverSettings(solver);
//
//  // This ensures that ghost cells of sol are correctly filled when returned
//  // from the solver
//  solver.setFinalFillBC(true);

  // Allocate fluxes
  Vector<Array<MultiFab*, 3>> fluxes(finest_level+1);

  for (int lev(0); lev <= finest_level; ++lev) {
    for(int dir = 0; dir < 3; dir++) {
      BoxArray edge_ba = amrex::convert(grids[lev], IntVect::TheDimensionVector(dir));
      fluxes[lev][dir] = new MultiFab(edge_ba, dmap[lev], nspecies_g, 1, MFInfo(), *ebfactory[lev]);
      fluxes[lev][dir]->setVal(0.);
    }
  }

//  setSolverSettings(solver);
//
//  // This ensures that ghost cells of sol are correctly filled when returned
//  // from the solver
//  solver.setFinalFillBC(true);

  // Copy X_gk MultiFabs into temporary variables
  Vector<MultiFab*> X_gk_copy(finest_level+1);
  for (int lev(0); lev <= finest_level; ++lev) {
    X_gk_copy[lev] = (MFHelpers::createFrom(*X_gk_in[lev])).release();
  }

  // Compute fluxes
  solver.getFluxes(fluxes, X_gk_copy, MLLinOp::Location::FaceCentroid);

  for (int lev(0); lev <= finest_level; ++lev) {
    delete X_gk_copy[lev];
  }

  // Correct the first term computed
  for (int lev(0); lev <= finest_level; ++lev) {
    // Auxiliary data
    Array<MultiFab*, 3> X_gk_faces;

    for(int dir = 0; dir < 3; dir++) {
      BoxArray edge_ba = amrex::convert(grids[lev], IntVect::TheDimensionVector(dir));
      X_gk_faces[dir] = new MultiFab(edge_ba, dmap[lev], nspecies_g, 1, MFInfo(), *ebfactory[lev]);
      X_gk_faces[dir]->setVal(0.);
    }

    // Interpolate
    EB_interp_CellCentroid_to_FaceCentroid(*X_gk_in[lev], X_gk_faces, 0, 0,
                                           nspecies_g, geom[lev], bcs_dummy);

    // Compute fluxes_gk = X_gk_faces sum{fluxes_gk}
    for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*fluxes[lev][dir],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.growntilebox(IntVect::TheDimensionVector(dir));

        Array4<Real const> const& X_gk_faces_arr = X_gk_faces[dir]->const_array(mfi);
        Array4<Real      > const& fluxes_arr     = fluxes[lev][dir]->array(mfi);

        amrex::ParallelFor(bx, [X_gk_faces_arr,fluxes_arr,nspecies_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real sum(0);

          for (int n(0); n < nspecies_g; ++n) {
            sum += fluxes_arr(i,j,k,n);
          }

          for (int n(0); n < nspecies_g; ++n) {
            fluxes_arr(i,j,k,n) = X_gk_faces_arr(i,j,k,n)*sum;
          }
        });
      } // MFIter
    } // direction loop

    // Data for storing the divergence of the auxiliary data
    MultiFab DivXGX_aux(grids[lev], dmap[lev], nspecies_g, nghost, MFInfo(), *ebfactory[lev]);

    // Set the DivXGX_aux data to 0
    DivXGX_aux.setVal(0.);

    // Compute the divergence
    EB_computeDivergence(DivXGX_aux, GetArrOfConstPtrs(fluxes[lev]), geom[lev],
        already_on_centroids);

    // Data for storing the divergence
    MultiFab DivXGX(grids[lev], dmap[lev], nspecies_g, nghost, MFInfo(), *ebfactory[lev]);

    // Set the DivXGX data to 0
    DivXGX.setVal(0.);

    // Redistribute the divergence
    amrex::single_level_redistribute(DivXGX_aux, DivXGX, 0, nspecies_g, geom[lev]);

    // Set divergence covered values
    EB_set_covered(DivXGX, 0, nspecies_g, DivXGX.nGrow(), 0.);

    // Correction: X_gk_in - div(X_gk_in sum{j_gk})
    MultiFab::Saxpy(*X_gk_in[lev], dt, DivXGX, 0, 0, nspecies_g, 1);

    // Free up space
    for(int dir = 0; dir < 3; dir++) {
      delete fluxes[lev][dir];
      delete X_gk_faces[dir];
    }
  } // lev
}


void DiffusionOp::ComputeLaphX (const Vector< MultiFab*       >& laphX_out,
                                const Vector< MultiFab*       >& X_gk_in,
                                const Vector< MultiFab const* >& ro_g_in,
                                const Vector< MultiFab const* >& ep_g_in,
                                const Vector< MultiFab const* >& T_g_in)
{
  BL_PROFILE("DiffusionOp::ComputeLaphX");

  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

//  // TODO: check on this
//  const bool already_on_centroids = true;

  int finest_level = amrcore->finestLevel();

  // Number of fluid species
  const int nspecies_g = fluid.nspecies;

  // Auxiliary data where we store Div{ep_g ro_g h_gk D_gk Grad{X_gk}}
  Vector< MultiFab* > laphX_aux(finest_level+1);

  // Allocate space for laphX_aux and set it to 0
  for(int lev = 0; lev <= finest_level; lev++)
  {
    laphX_aux[lev] = new MultiFab(grids[lev], dmap[lev], nspecies_g, nghost, MFInfo(),
        *ebfactory[lev]);

    laphX_aux[lev]->setVal(0.0);
  }

  // We want to return div (ep_g ro_g h_gk D_gk grad)) phi
  species_matrix->setScalars(0.0, -1.0);

  Vector<BCRec> bcs_dummy; // This is just to satisfy the call to EB_interp...
  bcs_dummy.resize(3*nspecies_g);

  // Compute the coefficients
  for (int lev = 0; lev <= finest_level; lev++)
  {
    // Local temporary data in case h_gk is not null
    MultiFab hb_coeffs(ep_g_in[lev]->boxArray(), ep_g_in[lev]->DistributionMap(),
        nspecies_g, 1, MFInfo(), ep_g_in[lev]->Factory());

    hb_coeffs.setVal(0.);

    auto& fluid_parms = *fluid.parameters;

    // hb_coeffs = ep_g ro_g h_gk D_gk
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.growntilebox(IntVect(1,1,1));

      Array4<Real      > const& hb_coeffs_arr = hb_coeffs.array(mfi);
      Array4<Real const> const& ep_g_arr      = ep_g_in[lev]->const_array(mfi);
      Array4<Real const> const& ro_g_arr      = ro_g_in[lev]->const_array(mfi);
      Array4<Real const> const& T_g_arr       = T_g_in[lev]->const_array(mfi);

      amrex::ParallelFor(bx, [ep_g_arr,ro_g_arr,T_g_arr,hb_coeffs_arr,
          nspecies_g,fluid_parms,run_on_device]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const Real ep_g = ep_g_arr(i,j,k);
        const Real ro_g = ro_g_arr(i,j,k);
        const Real T_g  = T_g_arr(i,j,k);

        for (int n(0); n < nspecies_g; ++n) {
          const Real D_gk = run_on_device ?
            fluid_parms.calc_D_gk<RunOn::Device>(T_g,n) :
            fluid_parms.calc_D_gk<RunOn::Host>(T_g,n);

          const Real val = ep_g*ro_g*D_gk;

          hb_coeffs_arr(i,j,k,n) = run_on_device ?
            fluid_parms.calc_h_gk<RunOn::Device>(T_g,n) * val :
            fluid_parms.calc_h_gk<RunOn::Host>(T_g,n) * val;
        }
      });
    }

    // if h_gk is nullptr  species_b = b_coeffs
    // else                species_b = hb_coeffs
    EB_interp_CellCentroid_to_FaceCentroid (hb_coeffs, GetArrOfPtrs(species_b[lev]), 0, 0,
                                            nspecies_g, geom[lev], bcs_dummy);

    // Set BCoeffs
    species_matrix->setBCoeffs(lev, GetArrOfConstPtrs(species_b[lev]), MLMG::Location::FaceCentroid);

    // Set LevelBC
    species_matrix->setLevelBC(lev, GetVecOfConstPtrs(X_gk_in)[lev]);
  }

  {
    MLMG solver(*species_matrix);

    // Compute div (ep_g ro_g [h_gk] D_gk grad)) phi
    solver.apply(laphX_aux, X_gk_in);
  }

  // CORRECT FLUXES
  // Compute Fluxes for correcting the result
  // div(ep_g h_gk ro_g D_gk grad(phi_gk)) - div(phi_gk sum(ep_g ro_g D_gm grad(phi_gm)))
  for (int species_k(0); species_k < nspecies_g; ++species_k) {

    // Auxiliary data where we store Div{ep_g ro_g h_gk D_gk Grad{X_gk}}
    Vector< MultiFab* > correction_aux(finest_level+1);

    // Allocate space for laphX_aux and set it to 0
    for(int lev = 0; lev <= finest_level; lev++)
    {
      correction_aux[lev] = new MultiFab(grids[lev], dmap[lev], nspecies_g, nghost, MFInfo(),
                                         *ebfactory[lev]);

      correction_aux[lev]->setVal(0.0);
    }

    // We want to return div (ep_g ro_g h_gk X_gk D_gm grad)) phi
    species_matrix->setScalars(0.0, -1.0);

    Vector<BCRec> loc_bcs_dummy; // This is just to satisfy the call to EB_interp...
    loc_bcs_dummy.resize(3*nspecies_g);

    // Compute the coefficients
    for (int lev = 0; lev <= finest_level; lev++)
    {
      // Local temporary data in case h_gk is not null
      MultiFab hXb_coeffs(ep_g_in[lev]->boxArray(), ep_g_in[lev]->DistributionMap(),
                          nspecies_g, 1, MFInfo(), ep_g_in[lev]->Factory());

      hXb_coeffs.setVal(0.);

      auto& fluid_parms = *fluid.parameters;

      // hXb_coeffs = ep_g ro_g h_gk X_gk D_gm
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*ep_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.growntilebox(IntVect(1,1,1));

        Array4<Real      > const& hXb_coeffs_arr = hXb_coeffs.array(mfi);
        Array4<Real const> const& ep_g_arr       = ep_g_in[lev]->const_array(mfi);
        Array4<Real const> const& ro_g_arr       = ro_g_in[lev]->const_array(mfi);
        Array4<Real const> const& T_g_arr        = T_g_in[lev]->const_array(mfi);
        Array4<Real const> const& X_gk_arr       = X_gk_in[lev]->const_array(mfi);

        amrex::ParallelFor(bx, [ep_g_arr,ro_g_arr,T_g_arr,hXb_coeffs_arr,X_gk_arr,
            nspecies_g,fluid_parms,species_k,run_on_device]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const Real ep_g = ep_g_arr(i,j,k);
          const Real ro_g = ro_g_arr(i,j,k);
          const Real T_g  = T_g_arr(i,j,k);

          const Real hgk  = run_on_device ?
            fluid_parms.calc_h_gk<RunOn::Device>(T_g,species_k) :
            fluid_parms.calc_h_gk<RunOn::Host>(T_g,species_k);

          const Real Xgk  = X_gk_arr(i,j,k,species_k);

          for (int m(0); m < nspecies_g; ++m) {
            const Real val = ep_g*ro_g*hgk*Xgk;

            hXb_coeffs_arr(i,j,k,m) = run_on_device ?
              val*fluid_parms.calc_D_gk<RunOn::Device>(T_g,m) :
              val*fluid_parms.calc_D_gk<RunOn::Host>(T_g,m);
          }
        });
      }

      for(int dir = 0; dir < 3; dir++) {
        species_b[lev][dir]->setVal(0);
      }

      // if h_gk is nullptr  species_b = b_coeffs
      // else                species_b = hb_coeffs
      EB_interp_CellCentroid_to_FaceCentroid (hXb_coeffs, GetArrOfPtrs(species_b[lev]), 0, 0,
                                              nspecies_g, geom[lev], loc_bcs_dummy);

      // Set BCoeffs
      species_matrix->setBCoeffs(lev, GetArrOfConstPtrs(species_b[lev]), MLMG::Location::FaceCentroid);

      // Set LevelBC
      species_matrix->setLevelBC(lev, GetVecOfConstPtrs(X_gk_in)[lev]);
    } // lev

    MLMG solver(*species_matrix);

    // Compute div (ep_g ro_g [h_gk] D_gk grad)) phi
    solver.apply(correction_aux, X_gk_in);

    for (int lev(0); lev <= finest_level; ++lev) {

      for (int m(0); m < nspecies_g; ++m) {
        MultiFab::Subtract(*laphX_aux[lev], *correction_aux[lev], m, species_k, 1, laphX_aux[lev]->nGrow());
      }
    }

    for (int lev(0); lev <= finest_level; ++lev) {
      delete correction_aux[lev];
    }

  } // correct_fluxes

  // Redistribute laphX_aux into laphX_out
  for(int lev = 0; lev <= finest_level; lev++)
  {
    amrex::single_level_redistribute(*laphX_aux[lev], *laphX_out[lev], 0, nspecies_g, geom[lev]);
    EB_set_covered(*laphX_aux[lev], 0, laphX_aux[lev]->nComp(), laphX_aux[lev]->nGrow(), 0.);
  }

  // Free laphX_aux memory
  for(int lev = 0; lev <= finest_level; lev++)
  {
    delete laphX_aux[lev];
  }
}
