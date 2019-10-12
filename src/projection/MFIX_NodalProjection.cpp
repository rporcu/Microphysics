#include <MFIX_NodalProjection.H>
#include <AMReX.H>
#include <mfix.H>
#include <AMReX_EBMultiFabUtil.H>
#include <MFIX_MFHelpers.H>



void
NodalProjection::define (const mfix* a_mfix,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_hi )
{
    m_mfix  = a_mfix;
    m_bc_lo = a_bc_lo;
    m_bc_hi = a_bc_hi;
    m_ok    = true;
}

//
// Perform projection:
//
//     vel = vel - grad(phi)/ro
//
//  where phi is the solution of
//
//   div( sigma * grad(phi) ) = div(ep*vel) + d(ep)/dt
//
//  ep, ro, vel, and d(ep)/dt are cell-centered
//
//  sigma = ep/ro is cell-centered as well
//
//  phi is node-centered
//
// If a_scale_factor is passed in, phi is return as phi/a_scale_factor
//
void
NodalProjection::project (      Vector< std::unique_ptr< amrex::MultiFab > >& a_vel,
                          const Vector< std::unique_ptr< amrex::MultiFab > >& a_ep,
                          const Vector< std::unique_ptr< amrex::MultiFab > >& a_ro,
                          const Vector< std::unique_ptr< amrex::MultiFab > >& a_depdt,
                                Real a_time, Real a_scale_factor )
{
    AMREX_ALWAYS_ASSERT(m_ok);
    BL_PROFILE("NodalProjection::project");

    // Setup object for projection and initialize value sof internals
    setup();

    // Compute RHS
    computeRHS(a_vel, a_ep, a_depdt, a_time);

    // Print value of d(ep)/dt + div(ep*u) before projection:
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        amrex::Print() << "AT LEVEL " << lev << " BEFORE NODAL PROJECTION: \n";
        amrex::Print() << "  max(abs(dep/dt + diveu)) = "
                       << m_rhs[lev]->norm0(0,0,false,true) << "\n";
    }

    // Compute and set matrix coefficients
    for (int lev(0); lev < a_ep.size(); ++lev)
    {
        // Compute the PPE coefficients = (ep / ro)
        MultiFab::Copy(*m_sigma[lev],*a_ep[lev],0,0,1,0);
        MultiFab::Divide(*m_sigma[lev],*a_ro[lev],0,0,1,0);

        // Set matrix coefficients
        m_matrix -> setSigma(lev, *m_sigma[lev]);
    }

    // Solve
    m_solver -> solve( GetVecOfPtrs(m_phi), GetVecOfConstPtrs(m_rhs), m_mg_rtol, m_mg_atol );

    // Get fluxes -- fluxes = - (ep/ro)*grad(phi)
    m_solver -> getFluxes( GetVecOfPtrs(m_fluxes) );

    // Perform projection
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        // The fluxes currently hold MINUS (ep/ro) * grad(phi) so we divide by ep
        MultiFab::Divide( *m_fluxes[lev], *a_ep[lev], 0, 0, 1, 0 );
        MultiFab::Divide( *m_fluxes[lev], *a_ep[lev], 0, 1, 1, 0 );
        MultiFab::Divide( *m_fluxes[lev], *a_ep[lev], 0, 2, 1, 0 );

        // Now we correct the velocity with MINUS (1/rho) * grad(phi),
        MultiFab::Add( *a_vel[lev], *m_fluxes[lev], 0, 0, 3, 0);

        // Account for scale factor
        m_fluxes[lev] -> mult(- 1.0/a_scale_factor, m_fluxes[lev]->nGrow() );

        // Finally we get rid of ro and MINUS so that m_fluxes = grad(phi)
        for (int n(0); n < 3; ++n)
            MultiFab::Multiply(*m_fluxes[lev], *a_ro[lev], 0, n, 1, m_fluxes[lev]->nGrow() );

        // Fill boundaries and apply scale factor to phi
        m_phi[lev] -> FillBoundary( m_mfix -> geom[lev].periodicity() );
        m_phi[lev] -> mult(1.0/a_scale_factor, m_fluxes[lev] -> nGrow());

    }

    // Compute RHS -- this is only needed to print out post projection values
    computeRHS(a_vel, a_ep, a_depdt, a_time);

    // Print value of d(ep)/dt + div(ep*u) before projection:
    for (int lev(0); lev < m_phi.size(); ++lev)
    {
        amrex::Print() << "AT LEVEL " << lev << " AFTER NODAL PROJECTION: \n";
        amrex::Print() << "  max(abs(dep/dt + diveu)) = "
                       << m_rhs[lev]->norm0(0,0,false,true) << "\n";    }

}


//
// Read from input file
//
void
NodalProjection::readParameters ()
{
    ParmParse pp("mfix");
    pp.query( "mg_verbose"             , m_mg_verbose );
    pp.query( "mg_cg_verbose"          , m_mg_cg_verbose );
    pp.query( "mg_maxiter"             , m_mg_maxiter );
    pp.query( "mg_cg_maxiter"          , m_mg_cg_maxiter );
    pp.query( "mg_rtol"                , m_mg_rtol );
    pp.query( "mg_atol"                , m_mg_atol );
    pp.query( "mg_max_coarsening_level", m_mg_max_coarsening_level );
    pp.query( "bottom_solver_type"     , m_bottom_solver_type );
}


//
// Setup object before solve
//
void
NodalProjection::setup ()
{
    BL_PROFILE("NodalProjection::setup");
    AMREX_ALWAYS_ASSERT(m_ok);

    readParameters();

    // Set number of levels
    int nlev( m_mfix -> grids.size() );

    // Resize member data if necessary
    if ( nlev != m_phi.size() )
    {
        m_phi.resize(nlev);
        m_fluxes.resize(nlev);
        m_sigma.resize(nlev);
        m_rhs.resize(nlev);
    }

    // Regrid if necessary
    int nghost(1);      // We use 1 ghost node only -- it should be enough

    bool need_regrid(false);  // if BA and DM changed on any level, we need
                              // to update the matrix and the solver as well

    for (int lev(0); lev < nlev; ++lev )
    {
        const auto& ba = m_mfix -> grids[lev];
        const auto& dm = m_mfix -> dmap[lev];
        const auto& eb = *(m_mfix -> ebfactory[lev]);

        if ( (m_phi[lev] == nullptr)                 ||
             (m_phi[lev] -> boxArray()        != ba) ||
             (m_phi[lev] -> DistributionMap() != dm)  )
        {
            // Cell-centered data
            m_fluxes[lev].reset(new MultiFab(ba, dm, 3, nghost, MFInfo(), eb));
            m_sigma[lev].reset(new MultiFab(ba, dm, 1, nghost, MFInfo(), eb));

            // Node-centered data
            const auto& ba_nd = amrex::convert(ba, IntVect{1,1,1});
            m_phi[lev].reset(new MultiFab(ba_nd, dm, 1, nghost, MFInfo(), eb));
            m_rhs[lev].reset(new MultiFab(ba_nd, dm, 1, nghost, MFInfo(), eb));

            need_regrid = true;
        }

    }

    // Setup matrix and solver
    if ( (m_matrix == nullptr) || need_regrid )
    {

        //
        // Setup Matrix
        //
        LPInfo                       info;
        info.setMaxCoarseningLevel(m_mg_max_coarsening_level);
        m_matrix.reset(new MLNodeLaplacian(m_mfix->geom, m_mfix->grids, m_mfix->dmap, info,
                                           GetVecOfConstPtrs(m_mfix->ebfactory)));

        m_matrix->setGaussSeidel(true);
        m_matrix->setHarmonicAverage(false);
        m_matrix->setDomainBC( m_bc_lo, m_bc_hi );

        //
        // Setup solver
        //
        m_solver.reset(new MLMG(*m_matrix));

        m_solver->setMaxIter(m_mg_maxiter);
        m_solver->setVerbose(m_mg_verbose);
        m_solver->setCGVerbose(m_mg_cg_verbose);
        m_solver->setCGMaxIter(m_mg_cg_maxiter);

        if (m_bottom_solver_type == "smoother")
        {
            m_solver->setBottomSolver(MLMG::BottomSolver::smoother);
        }
        else if (m_bottom_solver_type == "bicg")
        {
            m_solver->setBottomSolver(MLMG::BottomSolver::bicgstab);
        }
        else if (m_bottom_solver_type == "cg")
        {
            m_solver->setBottomSolver(MLMG::BottomSolver::cg);
        }
        else if (m_bottom_solver_type == "bicgcg")
        {
            m_solver->setBottomSolver(MLMG::BottomSolver::bicgcg);
        }
        else if (m_bottom_solver_type == "cgbicg")
        {
            m_solver->setBottomSolver(MLMG::BottomSolver::cgbicg);
        }
        else if (m_bottom_solver_type == "hypre")
        {
#ifdef AMREX_USE_HYPRE
            m_solver->setBottomSolver(MLMG::BottomSolver::hypre);
#else
            amrex::Abort("AMReX was not built with HYPRE support");
#endif
        }

    }

    // Initialize all variables
    for (int lev(0); lev < nlev; ++lev)
    {
        m_phi[lev] -> setVal(0.0);
        m_fluxes[lev] -> setVal(0.0);
        m_sigma[lev] -> setVal(1.0);
        m_rhs[lev] ->  setVal(0.0);
    }

}



//
// Compute RHS: div(ep*u) + d(ep)/dt
//
void
NodalProjection::computeRHS (       Vector< std::unique_ptr< amrex::MultiFab > >& a_vel,
                              const Vector< std::unique_ptr< amrex::MultiFab > >& a_ep,
                              const Vector< std::unique_ptr< amrex::MultiFab > >& a_depdt,
                                    Real a_time )
{
    AMREX_ALWAYS_ASSERT(m_ok);
    BL_PROFILE("NodalProjection::computeRHS");

    // Swap ghost cells and apply BCs to velocity -- we need to do this to make sure
    //      that the divergence operator can see inflow values
    m_mfix -> mfix_set_velocity_bcs(a_time, a_vel, 0);

    // Number of levels
    int nlevs(a_vel.size());

    // Note that the solver imposes the boundary conditions with the right scalings so we don't
    //      fill any ghost cells here.
    Vector< std::unique_ptr<MultiFab> > epu(nlevs);

    for (int lev(0); lev < nlevs; ++lev)
    {
        // We only need one ghost cell here -- so no need to make it bigger
        int nghost(1);
        epu[lev].reset(new MultiFab( a_vel[lev]->boxArray(), a_vel[lev]->DistributionMap(),
                                     a_vel[lev]->nComp()   , nghost , MFInfo(),
                                     a_vel[lev]->Factory() ) );

        epu[lev] -> setVal(1.e200);

        MultiFab::Copy(*epu[lev], *a_vel[lev], 0, 0, 3, epu[lev]->nGrow() );

        for (int n(0); n < 3; n++)
            MultiFab::Multiply( *epu[lev], *a_ep[lev], 0, n, 1, epu[lev]->nGrow() );

        epu[lev] -> FillBoundary( m_mfix -> geom[lev].periodicity() );

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        // Extrapolate Dirichlet values to ghost cells -- but do it differently in that
        //  no-slip walls are treated exactly like slip walls --
        // Note that this routine is essential to impose the correct inflow bc's on
        //  the product ep  * vel
        for (MFIter mfi((*epu[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Why are we using this instead of simply multiplying vel and ep with their BCs in place
            // already?
            m_mfix -> set_vec_bcs(lev, (*epu[lev])[mfi], m_mfix -> geom[lev].Domain() );
        }

        epu[lev] -> FillBoundary( m_mfix -> geom[lev].periodicity() );

        // We set these to zero because if the values in the covered cells are undefined,
        //   even though they are multiplied by zero in the divu computation, we can still get NaNs
        EB_set_covered(*epu[lev], 0, epu[lev]->nComp(), 1, 0.0);

    }

    // Restore velocities to carry Dirichlet values on faces -- Do we still need this?
    int extrap_dir_bcs = 0;
    m_mfix -> mfix_set_velocity_bcs(a_time, a_vel, extrap_dir_bcs);

    // Compute div(epu) + a_depdt
    m_matrix -> compRHS( GetVecOfPtrs(m_rhs),  GetVecOfPtrs(epu), {},
                         GetVecOfPtrs(a_depdt) );

}
