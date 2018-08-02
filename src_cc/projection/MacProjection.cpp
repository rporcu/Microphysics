#include <MacProjection.H>
#include <mfix_mac_F.H>
#include <mfix_proj_F.H>
#include <mfix_F.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_ParmParse.H>

// For multigrid
using namespace amrex;

// Define unit vectors for easily convert indeces
IntVect MacProjection::e_x(1,0,0);
IntVect MacProjection::e_y(0,1,0);
IntVect MacProjection::e_z(0,0,1);

//
//
// 
MacProjection::MacProjection (AmrCore* a_amrcore, int a_nghost )
{
    m_amrcore = a_amrcore;
    m_nghost  = a_nghost;

    read_inputs();
}


//
//
// 
MacProjection::~MacProjection ()
{}


//
//
//
const MultiFab&
MacProjection::get_phi (int lev)
{
    return *m_phi[lev];
}


//
//
//
void
MacProjection::read_inputs ()
{
    ParmParse pp("mac");
    
    // Option to control MGML behavior
    pp.query( "verbose", verbose );
    pp.query( "mg_verbose", m_mg_verbose );
    pp.query( "mg_cg_verbose", m_mg_cg_verbose );
    pp.query( "mg_max_iter",  m_mg_max_iter );
    pp.query( "mg_cg_maxiter",  m_mg_cg_maxiter );
    pp.query( "mg_max_fmg_iter", m_mg_max_fmg_iter );
    pp.query( "mg_rtol",  m_mg_rtol );
    pp.query( "mg_atol",  m_mg_atol );
}

//
// 
// 
void
MacProjection::set_bcs ( IArrayBox* a_bc_ilo, IArrayBox* a_bc_ihi,
			 IArrayBox* a_bc_jlo, IArrayBox* a_bc_jhi,
			 IArrayBox* a_bc_klo, IArrayBox* a_bc_khi )
{
    m_bc_ilo = a_bc_ilo;
    m_bc_ihi = a_bc_ihi;
    m_bc_jlo = a_bc_jlo;
    m_bc_jhi = a_bc_jhi;
    m_bc_klo = a_bc_klo;
    m_bc_khi = a_bc_khi;

    int bc_lo[3], bc_hi[3];
    Box domain( m_amrcore->Geom(0).Domain() );
    
    set_ppe_bc ( bc_lo, bc_hi,
		 domain.loVect(), domain.hiVect(),
		 &m_nghost,
		 m_bc_ilo->dataPtr(), m_bc_ihi->dataPtr(),
		 m_bc_jlo->dataPtr(), m_bc_jhi->dataPtr(),
		 m_bc_klo->dataPtr(), m_bc_khi->dataPtr() );

    m_lobc = {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]};
    m_hibc = {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]};
    
}


//
// redefine working arrays if amrcore has changed
//
void
MacProjection::update_internals ()
{
   
    if ( m_diveu.size() != ( m_amrcore->finestLevel()+1) )
    {
	m_diveu.resize(m_amrcore->finestLevel()+1);
	m_phi.resize(m_amrcore->finestLevel()+1);
	m_b.resize(m_amrcore->finestLevel()+1);
	m_gradphi.resize(m_amrcore->finestLevel()+1);
    }
    
    for (int lev=0; lev <= m_amrcore->finestLevel(); ++lev )
    {
	

	if ( m_diveu[lev] == nullptr ||
	     ! BoxArray::SameRefs(m_diveu[lev]->boxArray(), m_amrcore -> boxArray(lev) ) ||
	     ! DistributionMapping::SameRefs(m_diveu[lev]->DistributionMap(),
					     m_amrcore->DistributionMap(lev)) )
	{
	    
	    m_diveu[lev].reset( new MultiFab( m_amrcore -> boxArray(lev),
					      m_amrcore -> DistributionMap(lev),
					      1,m_nghost) );

	    m_phi[lev].reset( new MultiFab( m_amrcore -> boxArray(lev),
					    m_amrcore -> DistributionMap(lev),
					    1,m_nghost) );

	    //
	    // Staggered quantities
	    // NOTE: no ghost node for grad(phi)
	    //
	    m_b[lev].resize(3);
	   
	    BoxArray x_ba = m_amrcore -> boxArray(lev);
	    x_ba = x_ba.surroundingNodes(0);
	    m_b[lev][0].reset( new  MultiFab( x_ba, m_amrcore -> DistributionMap(lev),
					      1, m_nghost) );
	    m_gradphi[lev][0].reset( new  MultiFab( x_ba, m_amrcore -> DistributionMap(lev),
						    1, 0) );

	    BoxArray y_ba = m_amrcore -> boxArray(lev);
	    y_ba = y_ba.surroundingNodes(1);
	    m_b[lev][1].reset( new  MultiFab( y_ba, m_amrcore -> DistributionMap(lev),
					      1, m_nghost) );
	    m_gradphi[lev][1].reset( new  MultiFab( y_ba, m_amrcore -> DistributionMap(lev),
						    1, 0) );

	    BoxArray z_ba = m_amrcore -> boxArray(lev);
	    z_ba = z_ba.surroundingNodes(2);
	    m_b[lev][2].reset( new  MultiFab( z_ba, m_amrcore -> DistributionMap(lev),
					      1, m_nghost) );
	    m_gradphi[lev][2].reset( new  MultiFab( z_ba, m_amrcore -> DistributionMap(lev),
						    1, 0) );
	   
	   
	};
    }

}


//
// Computes the following decomposition:
// 
//    u + c*grad(phi)/ro = u*  with  dep/dt + div(ep*u) = 0
//
// Inputs:
// 
//   lev      = the AMR level
//   u,v,w    = the MAC velocity field to be projected
//   ep, epo  = the cell-centered volume fractions
//              at the beginning and end of time step
//   ro       = the cell-centered density
//   c        = a real constant (tipically c = dt)
//   dt       = the time step width to compute dep/dt
//
// Outputs:
//
//  phi       = the projection auxiliary function
//  u,v,w     = the PROJECTED MAC velocity field
//
// Notes:
//
//  phi is computed by solving
//
//       div(ep*grad(phi)/ro) = div(ep * u*)/c + (ep-epo)/dt/c
//
//  WARNING: this method returns the MAC velocity with up-to-date BCs in place
// 
void 
MacProjection::apply_projection ( Vector< std::unique_ptr<MultiFab> >& u, 
				  Vector< std::unique_ptr<MultiFab> >& v,
				  Vector< std::unique_ptr<MultiFab> >& w,
				  Vector< std::unique_ptr<MultiFab> >& ep,
				  Vector< std::unique_ptr<MultiFab> >& epo,
				  const Vector< std::unique_ptr<MultiFab> >& ro,
				  const Real c, const Real dt )
{
    BL_PROFILE("MacProjection::apply_projection()");

    if (verbose)
    	Print() << "MAC Projection:\n";

    // Check that everything is consistent with amrcore
    update_internals();
    
    // Compute right hand side, AKA div(ep * u)
    compute_diveu(u, v, w, ep, c);

    // Add dep/dt to diveu and store in diveu
    Real ocdt = 1.0 / (c*dt);
    
    for (int lev = 0; lev <= m_amrcore->finestLevel(); ++lev )
    {
	set_ccmf_bcs( lev, *epo[lev] ); // BCs for ep are already in place (check compute_diveu)
	MultiFab::Saxpy( *m_diveu[lev],  ocdt,  *ep[lev], 0, 0, 1, m_diveu[lev]->nGrow() ); 
	MultiFab::Saxpy( *m_diveu[lev], -ocdt, *epo[lev], 0, 0, 1, m_diveu[lev]->nGrow() ); 
    }    

    // Print out infos
    if (verbose)
    {
        Print() << " >> Before projection\n" ; 
	for (int lev=0; lev <= m_amrcore -> finestLevel(); ++lev )
	{
	    // Multiply by c because diveu contains (div(eu) + dep/dt)/c at this time
	    Print() << "  * On level "<< lev
		    << " max(abs(dep/dt + diveu)) = " << (m_diveu[lev] -> norm0())*c << "\n";
	}
    }
    
    // Compute beta coefficients ( div(beta*grad(phi)) = RHS )
    compute_b_coeff( u, v, w, ep, ro);
    
    // solve div(b*grad(phi)) = diveu/c
    solve_for_phi();

    project_velocity( u, v, w, ro, -c);

    for ( int lev=0; lev <= m_amrcore -> finestLevel() ; ++lev )
	set_velocity_bcs( lev, u, v, w );

    // Print out infos
    if (verbose)
    {
	compute_diveu(u, v, w, ep, 1.0);
	
        Print() << " >> After projection\n";  
	for (int lev=0; lev <= m_amrcore -> finestLevel(); ++lev )
	{
	    MultiFab::Saxpy( *m_diveu[lev],  dt,  *ep[lev], 0, 0, 1, m_diveu[lev]->nGrow() ); 
	    MultiFab::Saxpy( *m_diveu[lev], -dt, *epo[lev], 0, 0, 1, m_diveu[lev]->nGrow() );
	    
	    Print() << "  * On level "<< lev
		    << " max(abs(dep/dt + diveu)) = " << m_diveu[lev] -> norm0() << "\n";
	}
    }
}




//
// Computes the following decomposition:
// 
//    u + c*grad(phi)/ro = u*  with  div(ep*u) = 0
//
// Inputs:
// 
//   lev    = the AMR level
//   u,v,w  = the MAC velocity field to be projected
//   ep     = the cell-centered volume fraction
//   ro     = the cell-centered density
//   c      = a real constant (tipically c = dt)
//
// Outputs:
//
//  phi     = the projection auxiliary function
//  u,v,w   = the PROJECTED MAC velocity field
//
// Notes:
//
//  phi is computed by solving
//
//       div(ep*grad(phi)/ro) = div(ep * u*)/c
//
//  WARNING: this method returns the MAC velocity with up-to-date BCs in place
// 
void 
MacProjection::apply_projection ( Vector< std::unique_ptr<MultiFab> >& u, 
				  Vector< std::unique_ptr<MultiFab> >& v,
				  Vector< std::unique_ptr<MultiFab> >& w,
				  Vector< std::unique_ptr<MultiFab> >& ep,
				  const Vector< std::unique_ptr<MultiFab> >& ro,
				  amrex::Real c )
{
    BL_PROFILE("MacProjection::apply_projection()");

    if (verbose)
    	Print() << "MAC Projection:\n";

    // Check that everything is consistent with amrcore
    update_internals();
    
    // Compute right hand side, AKA div(ep * u)
    compute_diveu(u, v, w, ep, c);

    // Print out infos
    if (verbose)
    {
        Print() << " >> Before projection\n" ; 
	for (int lev=0; lev <= m_amrcore -> finestLevel(); ++lev )
	{
	    // Multiply by c because diveu contains div(eu)/c
	    Print() << "  * On level "<< lev
		    << " max(abs(diveu)) = " << (m_diveu[lev] -> norm0())*c << "\n";
	}
    }
    
    // Compute beta coefficients ( div(beta*grad(phi)) = RHS )
    compute_b_coeff( u, v, w, ep, ro);
    
    // solve div(b*grad(phi)) = diveu/c
    solve_for_phi();

    project_velocity( u, v, w, ro, -c);

    for ( int lev=0; lev <= m_amrcore -> finestLevel() ; ++lev )
	set_velocity_bcs( lev, u, v, w );

    // Print out infos
    if (verbose)
    {
	compute_diveu(u, v, w, ep, 1.0);
	
        Print() << " >> After projection\n";  
	for (int lev=0; lev <= m_amrcore -> finestLevel(); ++lev )
	    Print() << "  * On level "<< lev
		    << " max(abs(diveu)) = " << m_diveu[lev] -> norm0() << "\n";
       
    }
}




//
//  Evaluate div(eu)/c at cell centers
//
//  u,v,w  : staggered edge velocities
//  ep     : cell centered volume fraction
//  c      : scaling factor 
// 
void
MacProjection::compute_diveu (  Vector< std::unique_ptr<MultiFab> >& u,
				Vector< std::unique_ptr<MultiFab> >& v,
				Vector< std::unique_ptr<MultiFab> >& w,
				Vector< std::unique_ptr<MultiFab> >& ep,
				const Real c = 1.0 )
{
  

    for ( int lev=0; lev <= m_amrcore->finestLevel() ; ++lev )
    {
        // Initialize to 0 since we are multiplying by (1/c) before we set bc's
        m_diveu[lev]->setVal(0.);

	set_velocity_bcs( lev, u, v, w );
	set_ccmf_bcs( lev, *ep[lev] );
    
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(*m_diveu[lev],true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();
	    
	    compute_mac_diveu ( BL_TO_FORTRAN_BOX(bx),
				BL_TO_FORTRAN_ANYD((*m_diveu[lev])[mfi]),
				(*ep[lev])[mfi].dataPtr (),
				BL_TO_FORTRAN_ANYD((*u[lev])[mfi]),
				BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
				BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
				m_amrcore -> Geom(lev).CellSize() );  
	}

	m_diveu[lev] -> mult ( 1.0/c, 1 );
    }
}


//
// Set the BCs for velocity only
// 
void
MacProjection::set_velocity_bcs ( int lev,
				      Vector< std::unique_ptr<MultiFab> >& u,
				      Vector< std::unique_ptr<MultiFab> >& v,
				      Vector< std::unique_ptr<MultiFab> >& w )
{
    BL_PROFILE("MacProjection::set_velocity_bcs()");

    u[lev] -> FillBoundary( m_amrcore -> Geom(lev).periodicity() );
    v[lev] -> FillBoundary( m_amrcore -> Geom(lev).periodicity() );
    w[lev] -> FillBoundary( m_amrcore -> Geom(lev).periodicity() );
     
    Box domain( m_amrcore->Geom(lev).Domain() );
	
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*m_diveu[lev], true); mfi.isValid(); ++mfi)
    {
	const Box& bx = (*m_diveu[lev])[mfi].box();
	
	set_mac_velocity_bcs ( bx.loVect(), bx.hiVect(),
			       BL_TO_FORTRAN_ANYD((*u[lev])[mfi]),
			       BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
			       BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
			       m_bc_ilo->dataPtr(), m_bc_ihi->dataPtr(),
			       m_bc_jlo->dataPtr(), m_bc_jhi->dataPtr(),
			       m_bc_klo->dataPtr(), m_bc_khi->dataPtr(),
			       domain.loVect(), domain.hiVect(),
			       &m_nghost );
    }	
    
}

//
// Set Cell-Centered-Multifab BCs
// 
void
MacProjection::set_ccmf_bcs ( int lev, MultiFab& mf )
{

    Box domain( m_amrcore->Geom(lev).Domain() );

    if(!mf.boxArray().ixType().cellCentered())
        amrex::Error("MacProjection::set_ccmf_bcs() can only be used for cell-centered arrays!");

    // Impose periodic bc's at domain boundaries and fine-fine copies in the
    // interior It is essential that we do this before the call to fill_bc0
    // below since fill_bc0 can extrapolate out to fill ghost cells outside the
    // domain after we have filled ghost cells inside the domain, but doing
    // this call after fill_bc0 can't fill ghost cells from ghost cells.
    mf.FillBoundary(m_amrcore -> Geom(lev).periodicity());

    // Fill all cell-centered arrays with first-order extrapolation at domain
    // boundaries
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
        const Box& sbx = mf[mfi].box();

	fill_bc0( mf[mfi].dataPtr(), sbx.loVect(), sbx.hiVect(),
		  m_bc_ilo->dataPtr(), m_bc_ihi->dataPtr(),
		  m_bc_jlo->dataPtr(), m_bc_jhi->dataPtr(),
		  m_bc_klo->dataPtr(), m_bc_khi->dataPtr(),
		  domain.loVect(),  domain.hiVect(),
		  &m_nghost );
    }

    // Impose periodic bc's at domain boundaries and fine-fine copies in the
    // interior It's not 100% clear whether we need this call or not.  Worth
    // testing.
    mf.FillBoundary( m_amrcore -> Geom(lev).periodicity() );
}


//
// Compute phi by solving the poisson equation :
//
//                  div( b * grad(phi) ) = rhs  
// 
void
MacProjection::solve_for_phi () 
{
    BL_PROFILE("MacProjection::compute_phi()");
    
    // 
    // First define the matrix (operator).
    // Class MLABecLaplacian describes the following operator:
    //
    //       (alpha * a - beta * (del dot b grad)) phi
    //
    LPInfo                       info;
    MLABecLaplacian              matrix( m_amrcore -> Geom(),
					 m_amrcore -> boxArray(),
					 m_amrcore -> DistributionMap(),
					 info );

    // It is essential that we set MaxOrder of the solver to 2
    // if we want to use the standard phi(i)-phi(i-1) approximation
    // for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the
    // gradient at a Dirichlet boundary.
    matrix.setMaxOrder( 2 );
    matrix.setDomainBC( m_lobc, m_hibc );
    matrix.setScalars ( 0.0, -1.0 );

    for (int lev=0; lev <= m_amrcore->finestLevel(); ++lev)
    {
	std::array<MultiFab const*,AMREX_SPACEDIM>   b;

	// Copy the PPE coefficient into the proper data strutcure
	b[0] = GetVecOfConstPtrs( m_b[lev] )[0]; 
	b[1] = GetVecOfConstPtrs( m_b[lev] )[1];
	b[2] = GetVecOfConstPtrs( m_b[lev] )[2];
	
	matrix.setBCoeffs( lev, b );

	// By this point we must have filled the Dirichlet values of phi stored in the ghost cells
	m_phi[lev] -> setVal(0.0);
	matrix.setLevelBC ( lev, GetVecOfConstPtrs(m_phi)[lev] );
    }
    
    // 
    // Then setup the solver ----------------------
    //
    MLMG  solver(matrix);

    solver.setMaxIter(m_mg_max_iter);
    solver.setMaxFmgIter(m_mg_max_fmg_iter);
    solver.setVerbose(m_mg_verbose);
    solver.setCGVerbose(m_mg_cg_verbose);
    solver.setCGMaxIter(m_mg_cg_maxiter);

    // This ensures that ghost cells of phi are correctly filled when returned from the solver
    solver.setFinalFillBC(true);

    solver.solve( GetVecOfPtrs(m_phi), GetVecOfConstPtrs(m_diveu), m_mg_rtol, m_mg_atol );

 
    for (int lev=0; lev <= m_amrcore->finestLevel(); ++lev)
	m_phi[lev] -> FillBoundary(m_amrcore -> Geom(lev).periodicity());

    // Get the gradient of phi
    solver.getGradSolution( GetVecOfArrOfPtrs(m_gradphi) );

}


//
// Computes the staggered Poisson's operator coefficients:
//
//      bcoeff = ep/ro
//
// Values are edge-centered.
// 
void
MacProjection::compute_b_coeff ( const Vector< std::unique_ptr<MultiFab> >& u,
				 const Vector< std::unique_ptr<MultiFab> >& v,
				 const Vector< std::unique_ptr<MultiFab> >& w,
				 const Vector< std::unique_ptr<MultiFab> >& ep,
				 const Vector< std::unique_ptr<MultiFab> >& ro )
{
    BL_PROFILE("MacProjection::compute_b_coeff");

    // Directions
    int xdir = 1;
    int ydir = 2;
    int zdir = 3;

    for ( int lev=0; lev <= m_amrcore->finestLevel(); ++lev )
    {
    
#ifdef _OPENMP
#pragma omp parallel 
#endif
	for (MFIter mfi(*ep[lev],true); mfi.isValid(); ++mfi)
	{
	    // Boxes for staggered components
	    Box ubx = mfi.tilebox (e_x);
	    Box vbx = mfi.tilebox (e_y);
	    Box wbx = mfi.tilebox (e_z);

	    // X direction
	    compute_bcoeff_mac (BL_TO_FORTRAN_BOX(ubx),
				BL_TO_FORTRAN_ANYD((*(m_b[lev][0]))[mfi]),
				BL_TO_FORTRAN_ANYD((*u[lev])[mfi]),
				BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
				(*ep[lev])[mfi].dataPtr(), &xdir );

	    // Y direction
	    compute_bcoeff_mac (BL_TO_FORTRAN_BOX(vbx),
				BL_TO_FORTRAN_ANYD((*(m_b[lev][1]))[mfi]),
				BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
				BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
				(*ep[lev])[mfi].dataPtr(), &ydir );

	    // Z direction
	    compute_bcoeff_mac (BL_TO_FORTRAN_BOX(wbx),
				BL_TO_FORTRAN_ANYD((*(m_b[lev][2]))[mfi]),
				BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
				BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
				(*ep[lev])[mfi].dataPtr(), &zdir );
	}

	m_b[lev][0] -> FillBoundary( m_amrcore -> Geom(lev).periodicity() );
	m_b[lev][1] -> FillBoundary( m_amrcore -> Geom(lev).periodicity() );
	m_b[lev][2] -> FillBoundary( m_amrcore -> Geom(lev).periodicity() );
    }
}



// 
// Compute u = u* - c * grad(phi) / ro
// 
void
MacProjection::project_velocity ( Vector< std::unique_ptr<MultiFab> >& u,
				  Vector< std::unique_ptr<MultiFab> >& v,
				  Vector< std::unique_ptr<MultiFab> >& w,
				  const Vector< std::unique_ptr<MultiFab> >& ro,
				  const Real c )
{
    BL_PROFILE("MacProjection::project_velocity()");

    int xdir = 1;
    int ydir = 2;
    int zdir = 3;

    for (int lev=0; lev <= m_amrcore->finestLevel(); ++lev )
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(*m_phi[lev],true); mfi.isValid(); ++mfi)
	{
	    // Boxes for staggered components
	    Box ubx = mfi.tilebox (e_x);
	    Box vbx = mfi.tilebox (e_y);
	    Box wbx = mfi.tilebox (e_z);

	    project_mac_velocity(
		BL_TO_FORTRAN_BOX(ubx),
		BL_TO_FORTRAN_ANYD((*u[lev])[mfi]),
		BL_TO_FORTRAN_ANYD((*m_gradphi[lev][0])[mfi]),
		BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
		&c, &xdir );

	    project_mac_velocity(
		BL_TO_FORTRAN_BOX(vbx),
		BL_TO_FORTRAN_ANYD((*v[lev])[mfi]),
		BL_TO_FORTRAN_ANYD((*m_gradphi[lev][1])[mfi]),
		BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
		&c, &ydir );

	    
	    project_mac_velocity(
		BL_TO_FORTRAN_BOX(wbx),
		BL_TO_FORTRAN_ANYD((*w[lev])[mfi]),
		BL_TO_FORTRAN_ANYD((*m_gradphi[lev][2])[mfi]),
		BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
		&c, &zdir );	
	    
	}
    }
}


