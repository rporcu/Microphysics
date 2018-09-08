#include <AMReX_ParmParse.H>

#include <mfix_proj_F.H>
#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

// For multigrid
#include <AMReX_MLMG.H>
// #include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLNodeLaplacian.H>

//
// Computes the following decomposition:
//
//    u + grad(phi)/ro_g = u*,     where div(eps*u) = 0
//
// where u* is a non-div-free velocity field, stored
// by components in u_g, v_g, and w_g. The resulting div-free
// velocity field, u, overwrites the value of u* in u_g, v_g, and w_g.
//
// phi is an auxiliary function related to the pressure p_g by the relation:
//
//     new p_g  = old p_g + phi
void 
mfix_level::mfix_apply_projection ( int lev, amrex::Real scaling_factor, bool proj_2 )
{
    BL_PROFILE("mfix_level::mfix_apply_projection");

    // Swap ghost cells and apply BCs to velocity
    mfix_set_velocity_bcs (lev,0);

    // Print info about predictor step
    {
        amrex::Print() << "Before projection \n";
        mfix_print_max_vel (lev);
        mfix_compute_diveu (lev);
        amrex::Print() << "max(abs(diveu)) = " << mfix_norm0(diveu, lev, 0) << "\n";
    }

    // Here we add the (1/rho gradp) back to ustar (note the +dt)
    // We leave the (-1/rho gradp0) term in vel_g
    if (proj_2)
    {
       mfix_add_grad_phi ( lev,  scaling_factor, (*p_g[lev]));
       mfix_set_velocity_bcs (lev,0);
    }

    // Compute right hand side, AKA div(ep_g* u) / dt
    mfix_compute_diveu (lev);

    //
    // NOTE: WE MULTIPLY BY NEGATIVE 1/DT to account for negative in solver
    //
    diveu[lev] -> mult (-1.0/scaling_factor, diveu[lev]->nGrow() );

    // Compute the PPE coefficients
    mfix_compute_bcoeff_ppe ( lev );

    // Set BCs for Poisson's solver
    int bc_lo[3], bc_hi[3];
    Box domain(geom[lev].Domain());
    
    set_ppe_bc (bc_lo, bc_hi,
		domain.loVect(), domain.hiVect(),
		&nghost, 
		bc_ilo.dataPtr(), bc_ihi.dataPtr(),
		bc_jlo.dataPtr(), bc_jhi.dataPtr(),
		bc_klo.dataPtr(), bc_khi.dataPtr());

    // Initialize phi to zero (any non-zero bc's are stored in p0)
    phi[lev] -> setVal(0.);
 
    // Solve PPE
    MultiFab fluxes(vel_g[lev]->boxArray(), vel_g[lev]->DistributionMap(),
                    vel_g[lev]->nComp(), vel_g[lev]->nGrow(), MFInfo(),
                    *ebfactory[lev]);

    solve_poisson_equation( lev, bcoeff, phi, diveu, bc_lo, bc_hi, fluxes );
 
    // Correct the velocity field
    MultiFab::Divide( fluxes, *ep_g[lev], 0, 0, 1, 0 );
    MultiFab::Divide( fluxes, *ep_g[lev], 0, 1, 1, 0 );
    MultiFab::Divide( fluxes, *ep_g[lev], 0, 2, 1, 0 );

    //
    // NOTE: THE SIGN OF DT (scaling_factor) IS CORRECT HERE
    //
    amrex::Print() << "Multiplying fluxes by dt " << scaling_factor << std::endl;
    fluxes.mult ( -scaling_factor, fluxes.nGrow() );
    MultiFab::Add( *vel_g[lev], fluxes, 0, 0, 1, 0);

    if (proj_2)
    {
       // p := phi
       MultiFab::Copy (*p_g[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
    }
    else
    {
       // p := p + phi
       MultiFab::Add (*p_g[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
    }

    // Swap ghost cells and apply BCs to velocity
    mfix_set_velocity_bcs (lev,0);

    // Print info about predictor step
    {
        amrex::Print() << "After  projection \n";
        mfix_print_max_vel (lev);
        mfix_compute_diveu (lev);
        amrex::Print() << "max(abs(diveu)) = " <<  mfix_norm0(vel_g, lev, 0) << "\n";
    }
}

//
// Solve PPE:
//
//                  div( eps_g/rho * grad(phi) ) = div(eps_g*u) 
// 
void
mfix_level::solve_poisson_equation (  int lev,
				      Vector< Vector< std::unique_ptr<MultiFab> > >& b,
				      Vector< std::unique_ptr<MultiFab> >& this_phi,
				      Vector< std::unique_ptr<MultiFab> >& rhs,
				      int bc_lo[], int bc_hi[],
                                      MultiFab& fluxes )
{
    BL_PROFILE("mfix_level::solve_poisson_equation");
    
    if (nodal_pressure == 1)
    {

       // 
       // First define the matrix (operator).
       //
       //        (del dot b sigma grad)) phi
       //
       LPInfo                       info;
       MLNodeLaplacian              matrix(geom, grids, dmap, info);

       matrix.setGaussSeidel(true);
       matrix.setHarmonicAverage(false);
       matrix.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
   			    {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );

       matrix.setSigma(0, *(b[lev][0]));
       matrix.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::Sigma);

       // By this point we must have filled the Dirichlet values of phi stored in the ghost cells
       this_phi[lev]->setVal(0.);
       matrix.setLevelBC ( lev, GetVecOfConstPtrs(this_phi)[lev] );

       // 
       // Then setup the solver ----------------------
       //
       MLMG  solver(matrix);
	
       solver.setMaxIter (mg_max_iter);
       solver.setMaxFmgIter (mg_max_fmg_iter);
       solver.setVerbose (mg_verbose);
       solver.setCGVerbose (mg_cg_verbose);
       solver.setCGMaxIter (mg_cg_maxiter);

       // 
       // Finally, solve the system
       //
       solver.solve ( GetVecOfPtrs(this_phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol );

       this_phi[lev] -> FillBoundary (geom[lev].periodicity());

    } else {

       // 
       // First define the matrix (operator).
       // Class MLABecLaplacian describes the following operator:
       //
       //       (alpha * a - beta * (del dot b grad)) phi
       //
       LPInfo                       info;
       MLEBABecLap                  matrix(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));
       Vector<const MultiFab*>      tmp;
       array<MultiFab const*,AMREX_SPACEDIM>   b_tmp;

       // Copy the PPE coefficient into the proper data strutcure
       tmp = amrex::GetVecOfConstPtrs ( b[lev] ) ;
       b_tmp[0] = tmp[0];
       b_tmp[1] = tmp[1];
       b_tmp[2] = tmp[2];

       // It is essential that we set MaxOrder of the solver to 2
       // if we want to use the standard phi(i)-phi(i-1) approximation
       // for the gradient at Dirichlet boundaries.
       // The solver's default order is 3 and this uses three points for the
       // gradient at a Dirichlet boundary.
       matrix.setMaxOrder(2);

       // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
       matrix.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
   			    {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );
    
       matrix.setScalars ( 0.0, -1.0 );    
       matrix.setBCoeffs ( lev, b_tmp );

       // By this point we must have filled the Dirichlet values of phi stored in the ghost cells
       this_phi[lev]->setVal(0.);
       matrix.setLevelBC ( lev, GetVecOfConstPtrs(this_phi)[lev] );

       // 
       // Then setup the solver ----------------------
       //
       MLMG  solver(matrix);
	
       solver.setMaxIter (mg_max_iter);
       solver.setMaxFmgIter (mg_max_fmg_iter);
       solver.setVerbose (mg_verbose);
       solver.setCGVerbose (mg_cg_verbose);
       solver.setCGMaxIter (mg_cg_maxiter);

       // This ensures that ghost cells of phi are correctly filled when returned from the solver
       solver.setFinalFillBC(true);

       // 
       // Finally, solve the system
       //
       solver.solve ( GetVecOfPtrs(this_phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol );
       solver.getFluxes( {&fluxes} );

       this_phi[lev] -> FillBoundary (geom[lev].periodicity());
    }
}

//
// Computes bcoeff = ep_g/ro_g at the faces of the scalar cells
// 
void
mfix_level::mfix_compute_bcoeff_ppe (int lev)
{
    BL_PROFILE("mfix_level::mfix_compute_bcoeff_ppe");

    // Directions
    int xdir = 1;
    int ydir = 2;
    int zdir = 3;
    
    if (nodal_pressure == 1)
    {
#ifdef _OPENMP
#pragma omp parallel 
#endif
       for (MFIter mfi(*ro_g[lev],true); mfi.isValid(); ++mfi)
       {
	   // Cell-centered tilebox
	   Box bx = mfi.tilebox();

	   // X direction
	   compute_bcoeff_nd (BL_TO_FORTRAN_BOX(bx),
		              BL_TO_FORTRAN_ANYD((*(bcoeff[lev][0]))[mfi]),
		              BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		              (*ep_g[lev])[mfi].dataPtr(), &xdir );

	   // Y direction
	   compute_bcoeff_nd (BL_TO_FORTRAN_BOX(bx),
		              BL_TO_FORTRAN_ANYD((*(bcoeff[lev][1]))[mfi]),
		              BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		              (*ep_g[lev])[mfi].dataPtr(), &ydir );

	   // Z direction
	   compute_bcoeff_nd (BL_TO_FORTRAN_BOX(bx),
		              BL_TO_FORTRAN_ANYD((*(bcoeff[lev][2]))[mfi]),
		              BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		              (*ep_g[lev])[mfi].dataPtr(), &zdir );
	
       }
    } else {
#ifdef _OPENMP
#pragma omp parallel 
#endif
       for (MFIter mfi(*ro_g[lev],true); mfi.isValid(); ++mfi)
       {
	   // Tileboxes for staggered components
	   Box ubx = mfi.tilebox (e_x);
	   Box vbx = mfi.tilebox (e_y);
	   Box wbx = mfi.tilebox (e_z);

	   // X direction
	   compute_bcoeff_cc (BL_TO_FORTRAN_BOX(ubx),
		              BL_TO_FORTRAN_ANYD((*(bcoeff[lev][0]))[mfi]),
  		              BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		              (*ep_g[lev])[mfi].dataPtr(), &xdir );

	   // Y direction
	   compute_bcoeff_cc (BL_TO_FORTRAN_BOX(vbx),
		              BL_TO_FORTRAN_ANYD((*(bcoeff[lev][1]))[mfi]),
		              BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		              (*ep_g[lev])[mfi].dataPtr(), &ydir );

	   // Z direction
	   compute_bcoeff_cc (BL_TO_FORTRAN_BOX(wbx),
		              BL_TO_FORTRAN_ANYD((*(bcoeff[lev][2]))[mfi]),
		              BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
		              (*ep_g[lev])[mfi].dataPtr(), &zdir );
	
       }
    }

    bcoeff[lev][0] -> FillBoundary(geom[lev].periodicity());
    bcoeff[lev][1] -> FillBoundary(geom[lev].periodicity());
    bcoeff[lev][2] -> FillBoundary(geom[lev].periodicity());
}
