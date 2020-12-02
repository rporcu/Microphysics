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

using namespace amrex;

//
// Constructor:
// We set up everything which doesn't change between timesteps here
//
DiffusionOp::DiffusionOp (AmrCore* _amrcore,
                          Vector< const EBFArrayBoxFactory* >const& _ebfactory,
                          std::array<amrex::LinOpBCType,3> a_velbc_lo,
                          std::array<amrex::LinOpBCType,3> a_velbc_hi,
                          std::array<amrex::LinOpBCType,3> a_scalbc_lo,
                          std::array<amrex::LinOpBCType,3> a_scalbc_hi,
                          std::array<amrex::LinOpBCType,3> a_temperaturebc_lo,
                          std::array<amrex::LinOpBCType,3> a_temperaturebc_hi,
                          std::array<amrex::LinOpBCType,3> a_speciesbc_lo,
                          std::array<amrex::LinOpBCType,3> a_speciesbc_hi,
                          int _nghost)
{
    if(verbose > 0)
        amrex::Print() << "Constructing DiffusionOp class" << std::endl;

    nghost = _nghost;

    m_velbc_lo = a_velbc_lo;
    m_velbc_hi = a_velbc_hi;
    m_scalbc_lo = a_scalbc_lo;
    m_scalbc_hi = a_scalbc_hi;
    m_temperaturebc_lo = a_temperaturebc_lo;
    m_temperaturebc_hi = a_temperaturebc_hi;
    m_speciesbc_lo = a_speciesbc_lo;
    m_speciesbc_hi = a_speciesbc_hi;

    // Get inputs from ParmParse
    readParameters();

    // Actually do the setup work here
    setup(_amrcore, _ebfactory);

    // We default to Neumann bc's for scalarson EB walls
    eb_is_dirichlet = false;
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

    const int nspecies_g = FLUID::nspecies;

    // Resize and reset data
    b.resize(max_level + 1);

    phi.resize(max_level + 1);
    rhs.resize(max_level + 1);
    vel_eb.resize(max_level + 1);

    if (SPECIES::solve)
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
            b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], 1, nghost,
                                           MFInfo(), *ebfactory[lev]));
        }

        phi[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 1,
                                    MFInfo(), *ebfactory[lev]));

        // No ghost cells needed for rhs
        rhs[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0,
                                    MFInfo(), *ebfactory[lev]));

        vel_eb[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                       MFInfo(), *ebfactory[lev]));
        vel_eb[lev]->setVal(0.0);

        if (SPECIES::solve)
        {
          species_phi[lev].reset(new MultiFab(grids[lev], dmap[lev], nspecies_g, 1,
                                              MFInfo(), *ebfactory[lev]));

          // No ghost cells needed for rhs
          species_rhs[lev].reset(new MultiFab(grids[lev], dmap[lev], nspecies_g, 0,
                                              MFInfo(), *ebfactory[lev]));
        }
    }

    //
    // Define the matrix for the viscous tensor solve.
    //
    LPInfo info;
    info.setMaxCoarseningLevel(mg_max_coarsening_level);
    vel_matrix.reset(new MLEBTensorOp(geom, grids, dmap, info, ebfactory));

    // It is essential that we set MaxOrder to 2 if we want to use the standard
    // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the gradient.
    vel_matrix->setMaxOrder(2);

    // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
    vel_matrix->setDomainBC(m_velbc_lo, m_velbc_hi);

    //
    // Define the matrix for the scalar diffusion solve.
    //
    scal_matrix.reset(new MLEBABecLap(geom, grids, dmap, info, ebfactory));
    temperature_matrix.reset(new MLEBABecLap(geom, grids, dmap, info, ebfactory));

    if (SPECIES::solve) {
      species_matrix.reset(new MLEBABecLap(geom, grids, dmap, info, ebfactory, nspecies_g));
    }

    // It is essential that we set MaxOrder to 2 if we want to use the standard
    // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the gradient.
    scal_matrix->setMaxOrder(2);
    temperature_matrix->setMaxOrder(2);
    
    if (SPECIES::solve) {
      species_matrix->setMaxOrder(2);
    }

    // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
    scal_matrix->setDomainBC(m_scalbc_lo, m_scalbc_hi);
    temperature_matrix->setDomainBC(m_temperaturebc_lo, m_temperaturebc_hi);

    if (SPECIES::solve)
    {
      species_matrix->setDomainBC(m_speciesbc_lo, m_speciesbc_hi);
      
      species_b.resize(max_level + 1);

      for(int lev = 0; lev <= max_level; lev++)
      {
          for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
          {
              BoxArray edge_ba = grids[lev];
              edge_ba.surroundingNodes(dir);
              species_b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], nspecies_g, nghost,
                                                     MFInfo(), *ebfactory[lev]));
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
                                 const Vector< MultiFab* >& ro_in,
                                 const Vector< MultiFab* >& ep_in,
                                 const Vector< MultiFab* >& eta_in)
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
        // average_cellcenter_to_face( GetArrOfPtrs(b[lev]), *eta_in[lev], geom[lev] );
        EB_interp_CellCentroid_to_FaceCentroid (*eta_in[lev], GetArrOfPtrs(b[lev]), 0, 0, 1, geom[lev], bcs_s);

        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
             b[lev][dir]->FillBoundary(geom[lev].periodicity());

        vel_matrix->setShearViscosity  ( lev, GetArrOfConstPtrs(b[lev]), MLMG::Location::FaceCentroid);
        vel_matrix->setEBShearViscosity( lev, (*eta_in[lev]));
        vel_matrix->setLevelBC         ( lev, GetVecOfConstPtrs(vel_in)[lev] );
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

       // Divide by density
       for (int n = 0; n < 3; n++)
       {
           MultiFab::Divide( *divtau_out[lev], *ro_in[lev], 0, n, 1, 0 );
           MultiFab::Divide( *divtau_out[lev], *ep_in[lev], 0, n, 1, 0 );
       }
    }

    for(int lev = 0; lev <= finest_level; lev++)
       delete divtau_aux[lev];
}

void DiffusionOp::ComputeLapT (const Vector< MultiFab* >& lapT_out,
                               const Vector< MultiFab* >& T_g,
                               const Vector< MultiFab* >& ep_g,
                               const Vector< MultiFab* >& k_g,
                               const Vector< MultiFab* >& T_g_on_eb,
                               const Vector< MultiFab* >& k_g_on_eb)
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
  for (int lev = 0; lev <= finest_level; lev++)
  {
    MultiFab ep_g_k_g(ep_g[lev]->boxArray(), ep_g[lev]->DistributionMap(), 1, 1,
        MFInfo(), ep_g[lev]->Factory());

    // Initialize to 0
    ep_g_k_g.setVal(0.);

    MultiFab::Copy(ep_g_k_g, *ep_g[lev], 0, 0, 1, 1);
    MultiFab::Multiply(ep_g_k_g, *k_g[lev], 0, 0, 1, 1);

    EB_interp_CellCentroid_to_FaceCentroid (ep_g_k_g, GetArrOfPtrs(b[lev]), 0,
        0, 1, geom[lev], bcs_s);

    if (EB::fix_temperature) {
      // The following is a WIP in AMReX
      //temperature_matrix->setPhiOnCentroid();
      temperature_matrix->setEBDirichlet(lev, *T_g_on_eb[lev], *k_g_on_eb[lev]);
    }

    temperature_matrix->setBCoeffs(lev, GetArrOfConstPtrs(b[lev]),
        MLMG::Location::FaceCentroid);

    temperature_matrix->setLevelBC(lev, GetVecOfConstPtrs(T_g)[lev]);
  }

  MLMG solver(*temperature_matrix);

  solver.apply(lapT_aux, T_g);

  for(int lev = 0; lev <= finest_level; lev++)
  {
    amrex::single_level_redistribute(*lapT_aux[lev], *lapT_out[lev], 0, 1,
        geom[lev]);
  }

  for(int lev = 0; lev <= finest_level; lev++)
  {
    delete lapT_aux[lev];
  }
}

void DiffusionOp::ComputeLapS (const Vector< MultiFab* >& laps_out,
                               const Vector< MultiFab* >& scal_in,
                               const Vector< MultiFab* >& ro_in,
                               const Vector< MultiFab* >& ep_in,
                               const Vector< Real      >& mu_s)
{
    BL_PROFILE("DiffusionOp::ComputeLapS");

    int finest_level = amrcore->finestLevel();

    int ntrac = scal_in[0]->nComp();

    Vector< MultiFab* >  laps_aux(finest_level+1);
    Vector< MultiFab* >    phi_eb(finest_level+1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
       laps_aux[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, nghost,
                                    MFInfo(), *ebfactory[lev]);
       laps_aux[lev]->setVal(0.0);

       phi_eb[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, 0,
                                  MFInfo(), *ebfactory[lev]);

       // This value was just for testing
       // if (eb_is_dirichlet)
       //    phi_eb[lev]->setVal(1.0);
    }

    // We want to return div (mu grad)) phi
    scal_matrix->setScalars(0.0, -1.0);

    // Compute the coefficients
    for (int lev = 0; lev <= finest_level; lev++)
    {
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

    for(int lev = 0; lev <= finest_level; lev++)
    {
       amrex::single_level_redistribute(*laps_aux[lev], *laps_out[lev], 0, ntrac, geom[lev]);
    }

    for(int lev = 0; lev <= finest_level; lev++)
    {
       delete laps_aux[lev];
       delete   phi_eb[lev];
    }
}

void DiffusionOp::ComputeLapX (const Vector< MultiFab* >& lapX_out,
                               const Vector< MultiFab* >& X_gk_in,
                               const Vector< MultiFab* >& ro_g_in,
                               const Vector< MultiFab* >& ep_g_in,
                               const Vector< MultiFab* >& D_gk_in)
{
  BL_PROFILE("DiffusionOp::ComputeLapX");

  int finest_level = amrcore->finestLevel();

  // Number of fluid species
  const int nspecies_g = FLUID::nspecies;

  Vector< MultiFab* > lapX_aux(finest_level+1);

  for(int lev = 0; lev <= finest_level; lev++)
  {
    lapX_aux[lev] = new MultiFab(grids[lev], dmap[lev], nspecies_g, nghost, MFInfo(),
        *ebfactory[lev]);

    lapX_aux[lev]->setVal(0.0);
  }

  // We want to return div (ep_g ro_g D_gk grad)) phi
  species_matrix->setScalars(0.0, -1.0);

  Vector<BCRec> bcs_X; // This is just to satisfy the call to EB_interp...
  bcs_X.resize(3*nspecies_g);

  // Compute the coefficients
  for (int lev = 0; lev <= finest_level; lev++)
  {
    MultiFab b_coeffs(ep_g_in[lev]->boxArray(), ep_g_in[lev]->DistributionMap(),
        nspecies_g, 1, MFInfo(), ep_g_in[lev]->Factory());

    b_coeffs.setVal(0.);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.growntilebox(IntVect(1,1,1));

      if (bx.ok())
      {
        Array4<Real const> const& ep_g_arr     = ep_g_in[lev]->const_array(mfi);
        Array4<Real const> const& ro_g_arr     = ro_g_in[lev]->const_array(mfi);
        Array4<Real const> const& D_gk_arr     = D_gk_in[lev]->const_array(mfi);
        Array4<Real      > const& b_coeffs_arr = b_coeffs.array(mfi);

        amrex::ParallelFor(bx, [ep_g_arr,ro_g_arr,D_gk_arr,b_coeffs_arr,nspecies_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const Real ep_g = ep_g_arr(i,j,k);
          const Real ro_g = ro_g_arr(i,j,k);

          for (int n(0); n < nspecies_g; ++n)
            b_coeffs_arr(i,j,k,n) = ep_g*ro_g*D_gk_arr(i,j,k,n);
        });
      }
    }

    EB_interp_CellCentroid_to_FaceCentroid (b_coeffs, GetArrOfPtrs(species_b[lev]), 0,
        0, nspecies_g, geom[lev], bcs_X);

    species_matrix->setBCoeffs(lev, GetArrOfConstPtrs(species_b[lev]), MLMG::Location::FaceCentroid);

    species_matrix->setLevelBC(lev, GetVecOfConstPtrs(X_gk_in)[lev]);
  }

  MLMG solver(*species_matrix);

  solver.apply(lapX_aux, X_gk_in);

  for(int lev = 0; lev <= finest_level; lev++)
  {
    amrex::single_level_redistribute(*lapX_aux[lev], *lapX_out[lev], 0, nspecies_g, geom[lev]);
  }

  for(int lev = 0; lev <= finest_level; lev++)
  {
    delete lapX_aux[lev];
  }

}
