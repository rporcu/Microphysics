#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB_utils.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <DiffusionOp.H>

using namespace amrex;

//
// Constructor:
// We set up everything which doesn't change between timesteps here
//
DiffusionOp::DiffusionOp (AmrCore* _amrcore,
                          Vector< const EBFArrayBoxFactory* >* _ebfactory,
                          std::array<amrex::LinOpBCType,3> a_velbc_lo,
                          std::array<amrex::LinOpBCType,3> a_velbc_hi,
                          std::array<amrex::LinOpBCType,3> a_scalbc_lo,
                          std::array<amrex::LinOpBCType,3> a_scalbc_hi,
                          int _nghost)
{
    if(verbose > 0)
        amrex::Print() << "Constructing DiffusionOp class" << std::endl;

    nghost = _nghost;

    m_velbc_lo = a_velbc_lo;
    m_velbc_hi = a_velbc_hi;
    m_scalbc_lo = a_scalbc_lo;
    m_scalbc_hi = a_scalbc_hi;

    // Get inputs from ParmParse
    readParameters();

    // Actually do the setup work here
    setup(_amrcore, _ebfactory);
}

void DiffusionOp::setup (AmrCore* _amrcore, 
                         Vector< const EBFArrayBoxFactory* >* _ebfactory)
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

    // Resize and reset data
    b.resize(max_level + 1);
    phi.resize(max_level + 1);
    rhs.resize(max_level + 1);
    vel_eb.resize(max_level + 1);

    for(int lev = 0; lev <= max_level; lev++)
    {
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            BoxArray edge_ba = grids[lev];
            edge_ba.surroundingNodes(dir);
            b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], 1, nghost,
                                           MFInfo(), *(*ebfactory)[lev]));
        }
        phi[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 1,
                                    MFInfo(), *(*ebfactory)[lev]));

        // No ghost cells needed for rhs
        rhs[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0,
                                    MFInfo(), *(*ebfactory)[lev]));

        vel_eb[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                       MFInfo(), *(*ebfactory)[lev]));
        vel_eb[lev]->setVal(0.0);
    }

    //
    // Define the matrix for the viscous tensor solve.
    //
    LPInfo info;
    info.setMaxCoarseningLevel(mg_max_coarsening_level);
    vel_matrix.reset(new MLEBTensorOp(geom, grids, dmap, info, *ebfactory));

    // It is essential that we set MaxOrder to 2 if we want to use the standard
    // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the gradient.
    vel_matrix->setMaxOrder(2);

    // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
    vel_matrix->setDomainBC(m_velbc_lo, m_velbc_hi);

    //
    // Define the matrix for the scalar diffusion solve.
    //
    scal_matrix.reset(new MLEBABecLap(geom, grids, dmap, info, *ebfactory));

     // It is essential that we set MaxOrder to 2 if we want to use the standard
    // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the gradient.
    scal_matrix->setMaxOrder(2);

    // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
    scal_matrix->setDomainBC(m_scalbc_lo, m_scalbc_hi);
}

DiffusionOp::~DiffusionOp ()
{}

void DiffusionOp::readParameters ()
{
    ParmParse pp("diffusion");

    pp.query("verbose_solver", verbose);
    pp.query("verbose", mg_verbose);
    pp.query("bottom_verbose", mg_cg_verbose);
    pp.query("maxiter", mg_maxiter);
    pp.query("bottom_maxiter", mg_cg_maxiter);
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
        solver.setCGMaxIter(mg_cg_maxiter);

        // Verbosity for MultiGrid / ConjugateGradients
        solver.setVerbose(mg_verbose);
        solver.setCGVerbose(mg_cg_verbose);

        // This ensures that ghost cells of phi are correctly filled when
        // returned from the solver
        solver.setFinalFillBC(true);
}

void DiffusionOp::ComputeDivTau (Vector< MultiFab* >& divtau_out,
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
                                      MFInfo(), *(*ebfactory)[lev]);
       divtau_aux[lev]->setVal(0.0);
    }
 
    // Whole domain
    Box domain(geom[0].Domain());
 
    // We want to return div (mu grad)) phi
    vel_matrix->setScalars(0.0, -1.0);
 
    // Compute the coefficients
    for (int lev = 0; lev <= finest_level; lev++)
    {
        average_cellcenter_to_face( GetArrOfPtrs(b[lev]), *eta_in[lev], geom[lev] );
 
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
             b[lev][dir]->FillBoundary(geom[lev].periodicity());
 
        vel_matrix->setShearViscosity  ( lev, GetArrOfConstPtrs(b[lev]));
        vel_matrix->setEBShearViscosity( lev, (*eta_in[lev]));
        vel_matrix->setLevelBC         ( lev, GetVecOfConstPtrs(vel_in)[lev] );
    }
 
    MLMG solver(*vel_matrix);
 
    solver.apply(divtau_aux, vel_in);
 
    for(int lev = 0; lev <= finest_level; lev++)
    {
       amrex::single_level_weighted_redistribute(lev,
                                                 *divtau_aux[lev],
                                                 *divtau_out[lev],
                                                 *ep_in[lev],
                                                 0,
                                                 AMREX_SPACEDIM,
                                                 geom);
 
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

void DiffusionOp::ComputeLapS (Vector< MultiFab* >& laps_out,
                               const Vector< MultiFab* >& scal_in,
                               const Vector< MultiFab* >& ro_in,
                               const Vector< MultiFab* >& ep_in,
                               const Vector< Real > mu_s) 
{
    BL_PROFILE("DiffusionOp::ComputeLapS");

    int finest_level = amrcore->finestLevel();

    int ntrac = scal_in[0]->nComp();

    Vector< MultiFab* > laps_aux(finest_level+1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
       laps_aux[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, nghost,
                                    MFInfo(), *(*ebfactory)[lev]);
       laps_aux[lev]->setVal(0.0);
    }
 
    // Whole domain
    Box domain(geom[0].Domain());
 
    // We want to return div (mu grad)) phi
    scal_matrix->setScalars(0.0, -1.0);
 
    // Compute the coefficients
    for (int lev = 0; lev <= finest_level; lev++)
    {
        for(int dir = 0; dir < 3; dir++)
           for(int n = 0; n < ntrac; n++)
             b[lev][dir]->setVal(mu_s[n],n,1);
 
        scal_matrix->setBCoeffs(lev, GetArrOfConstPtrs(b[lev]));
        scal_matrix->setLevelBC(lev, GetVecOfConstPtrs(scal_in)[lev]);
    }
 
    MLMG solver(*scal_matrix);
 
    solver.apply(laps_aux, scal_in);
 
    for(int lev = 0; lev <= finest_level; lev++)
    {
       amrex::single_level_redistribute(lev, *laps_aux[lev], *laps_out[lev], 0, ntrac, geom);
    }
    
    for(int lev = 0; lev <= finest_level; lev++)
       delete laps_aux[lev];
}
