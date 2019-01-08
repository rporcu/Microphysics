#include <AMReX_ParmParse.H>

#include <mfix_proj_F.H>
#include <mfix_mac_F.H>
#include <mfix_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

// For multigrid
#include <AMReX_MLMG.H>
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
//     new p_g  = phi
//
//     except in the initial iterations when
//
//     new p_g  = old p_g + phi
void
mfix::mfix_apply_projection ( amrex::Real time, amrex::Real scaling_factor, bool proj_2 )
{
    BL_PROFILE("mfix::mfix_apply_projection");


    // Set domain BCs for Poisson's solver
    // The domain BCs refer to level 0 only
    int bc_lo[3], bc_hi[3];
    Box domain(geom[0].Domain());

    set_ppe_bc(bc_lo, bc_hi,
               domain.loVect(), domain.hiVect(),
               &nghost,
               bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
               bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
               bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    // Swap ghost cells and apply BCs to velocity -- we need to do this to make sure
    //      that the divergence operator can see inflow values
    mfix_set_velocity_bcs(time, 0);

    mfix_compute_diveu(time);

    // Print info about predictor step
    for (int lev = 0; lev < nlev; lev++)
    {
        amrex::Print() << "AT LEVEL " << lev << " BEFORE PROJECTION: \n";
        mfix_print_max_vel(lev);
        mfix_print_max_gp (lev);
        amrex::Print() << "Min and Max of ep_g " << ep_g[lev]->min(0) << " "
                                                 << ep_g[lev]->max(0) << std::endl;
    }

    for (int lev = 0; lev < nlev; lev++)
    {
        amrex::Print() << "max(abs(diveu)) = " << mfix_norm0(diveu, lev, 0) << "\n";

        // Here we add (dt * (1/rho gradp)) to ustar
        if (proj_2)
        {
            // Convert velocities to momenta
            for (int n = 0; n < 3; n++)
                MultiFab::Multiply(*vel_g[lev],(*ro_g[lev]),0,n,1,vel_g[lev]->nGrow());

            MultiFab::Saxpy (*vel_g[lev], scaling_factor, *gp[lev], 0, 0, 3, vel_g[lev]->nGrow());

            // Convert momenta back to velocities
            for (int n = 0; n < 3; n++)
                MultiFab::Divide(*vel_g[lev],(*ro_g[lev]),0,n,1,vel_g[lev]->nGrow());

        }
    }

    mfix_set_velocity_bcs (time, 0);

    // Compute right hand side, AKA div(ep_g* u) / dt

    mfix_compute_diveu(time);

    for (int lev = 0; lev < nlev; lev++)
       amrex::Print() << "At level " << lev << ": max(abs(dive(u+gp/rho))) = " << mfix_norm0(diveu, lev, 0) << "\n";

    // Initialize phi to zero (any non-zero bc's are stored in p0)
    for (int lev = 0; lev < nlev; lev++)
        phi[lev] -> setVal(0.);

    // Compute the PPE coefficients = (ep_g / rho)
    mfix_compute_bcoeff_ppe( );

    Vector<std::unique_ptr<MultiFab> > fluxes;

    fluxes.resize(nlev);
    for (int lev = 0; lev < nlev; lev++)
    {
        // We don't need more than 1 ghost cell in fluxes so no need to make it bigger
        fluxes[lev].reset(new MultiFab(vel_g[lev]->boxArray(), vel_g[lev]->DistributionMap(),
                                       vel_g[lev]->nComp(), 1, MFInfo(),
                                       *ebfactory[lev]));

        fluxes[lev]->setVal(1.e200);
    }

    // Solve PPE
    solve_poisson_equation( bcoeff, phi, diveu, fluxes, bc_lo, bc_hi);

    for (int lev = 0; lev < nlev; lev++)
    {
        // The fluxes currently hold MINUS (dt) * (ep_g/rho) * grad(phi) so we divide by ep_g
        MultiFab::Divide( *fluxes[lev], *ep_g[lev], 0, 0, 1, 0 );
        MultiFab::Divide( *fluxes[lev], *ep_g[lev], 0, 1, 1, 0 );
        MultiFab::Divide( *fluxes[lev], *ep_g[lev], 0, 2, 1, 0 );

        // Now we correct the velocity with MINUS (dt) * (1/rho) * grad(phi),
        MultiFab::Add( *vel_g[lev], *fluxes[lev], 0, 0, 3, 0);

        // The fluxes currently hold MINUS (dt) * (1/rho) * grad(phi),
        // so now we multiply by rho and divide by (-dt) to get grad(phi)
        fluxes[lev]->mult ( -1/scaling_factor, fluxes[lev]->nGrow() );
        for (int n = 0; n < 3; n++)
            MultiFab::Multiply(*fluxes[lev],(*ro_g[lev]),0,n,1,fluxes[lev]->nGrow());

        // phi currently holds (dt) * phi so we divide by (dt) to get (phi)
        phi[lev]->mult ( 1/scaling_factor, phi[lev]->nGrow() );

        if (proj_2)
        {
            // p := phi
            MultiFab::Copy (*p_g[lev],    *phi[lev], 0, 0, 1,    phi[lev]->nGrow());
            MultiFab::Copy ( *gp[lev], *fluxes[lev], 0, 0, 3, fluxes[lev]->nGrow());
        }
        else
        {
            // p := p + phi
            MultiFab::Add (*p_g[lev],    *phi[lev], 0, 0, 1,    phi[lev]->nGrow());
            MultiFab::Add ( *gp[lev], *fluxes[lev], 0, 0, 3, fluxes[lev]->nGrow());
        }
    }

    for (int lev = nlev-1; lev > 0; lev--)
    {
        avgDown (lev-1, *vel_g[lev], *vel_g[lev-1]);
        avgDown (lev-1, *   gp[lev],    *gp[lev-1]);
    }

    // Swap ghost cells and apply BCs to velocity
    mfix_set_velocity_bcs (time, 0);

    mfix_compute_diveu(time);

    for (int lev = 0; lev < nlev; lev++)
    {
        amrex::Print() << "AT LEVEL " << lev << " AFTER PROJECTION: \n";
        mfix_print_max_vel(lev);
        mfix_print_max_gp (lev);
        amrex::Print() << "max(abs(diveu)) = " <<  mfix_norm0(diveu, lev, 0) << "\n";
    }
}

//
// Solve PPE:
//
//                  div( eps_g/rho * grad(phi) ) = div(eps_g*u)
//
void
mfix::solve_poisson_equation ( Vector< Vector< std::unique_ptr<MultiFab> > >& b,
             Vector< std::unique_ptr<MultiFab> >& this_phi,
             Vector< std::unique_ptr<MultiFab> >& rhs,
             Vector< std::unique_ptr<MultiFab> >& fluxes,
             int bc_lo[], int bc_hi[] )
{
    BL_PROFILE("mfix::solve_poisson_equation");

    //
    // First define the matrix (operator).
    //
    //        (del dot b sigma grad)) phi
    //
    LPInfo                       info;
    MLNodeLaplacian              matrix(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));

    matrix.setGaussSeidel(true);
    matrix.setHarmonicAverage(false);
    matrix.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                         {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );

    for (int lev = 0; lev < nlev; lev++)
      {
        matrix.setSigma(lev, *(b[lev][0]));

        // By this point we must have filled the Dirichlet values of phi stored in the ghost cells
        this_phi[lev]->setVal(0.);
        matrix.setLevelBC ( lev, GetVecOfConstPtrs(this_phi)[lev] );
      }

    //
    // Then setup the solver ----------------------
    //
    MLMG  solver(matrix);

    solver.setMaxIter (mg_max_iter);
    solver.setMaxFmgIter (mg_max_fmg_iter);
    solver.setVerbose (mg_verbose);
    solver.setCGVerbose (mg_cg_verbose);
    solver.setCGMaxIter (mg_cg_maxiter);

    // solver.setBottomSolver(MLMG::BottomSolver::bicgstab);
    // solver.setPreSmooth (20);
    // solver.setPostSmooth (20);

    //
    // Finally, solve the system
    //
    solver.solve ( GetVecOfPtrs(this_phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol );
    solver.getFluxes( amrex::GetVecOfPtrs(fluxes) );

    for (int lev = 0; lev < nlev; lev++)
      this_phi[lev] -> FillBoundary(geom[lev].periodicity());
}

//
// Computes bcoeff = ep_g/ro_g at the faces of the scalar cells
//
void
mfix::mfix_compute_bcoeff_ppe ()
{
    BL_PROFILE("mfix::mfix_compute_bcoeff_ppe");

    // Directions
    int xdir = 1;
    int ydir = 2;
    int zdir = 3;

    for (int lev = 0; lev < nlev; lev++)
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

      bcoeff[lev][0] -> FillBoundary(geom[lev].periodicity());
      bcoeff[lev][1] -> FillBoundary(geom[lev].periodicity());
      bcoeff[lev][2] -> FillBoundary(geom[lev].periodicity());
    }
}
