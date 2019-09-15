#include <mfix.H>

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

    // Swap ghost cells and apply BCs to velocity -- we need to do this to make sure
    //      that the divergence operator can see inflow values
    mfix_set_velocity_bcs(time, vel_g, 0);

    mfix_compute_diveu(time);

    // Print info about predictor step
    for (int lev = 0; lev < nlev; lev++)
    {
        amrex::Print() << "AT LEVEL " << lev << " BEFORE PROJECTION: \n";
        mfix_print_max_vel(lev);
        mfix_print_max_gp (lev);
        amrex::Print() << "Min and Max of ep_g " << mfix_min(*ep_g[lev],lev,0) << " "
                                                 << mfix_max(*ep_g[lev],lev,0) << std::endl;
    }

    for (int lev = 0; lev < nlev; lev++)
    {
        amrex::Print() << "   max(abs(diveu)) = " << mfix_norm0(diveu, lev, 0) << "\n";

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

    mfix_set_velocity_bcs (time, vel_g, 0);

    // Compute right hand side, AKA div(ep_g* u) / dt

    mfix_compute_diveu(time);

    for (int lev = 0; lev < nlev; lev++)
       amrex::Print() << "At level " << lev << ": max(abs(dive(u+gp/rho))) = " << mfix_norm0(diveu, lev, 0) << "\n";

    // Initialize phi to zero (any non-zero bc's are stored in p0)
    for (int lev = 0; lev < nlev; lev++)
        phi_nd[lev] -> setVal(0.);

    // Compute the PPE coefficients = (ep_g / rho)
    for (int lev = 0; lev < nlev; lev++)
    {
         MultiFab::Copy  (*bcoeff_nd[lev],*ep_g[lev],0,0,1,0);
         MultiFab::Divide(*bcoeff_nd[lev],*ro_g[lev],0,0,1,0);
    }

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
    mfix_setup_nodal_solver ();
    mfix_solve_poisson_equation (phi_nd, diveu, bcoeff_nd, fluxes);

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
        phi_nd[lev]->mult ( 1/scaling_factor, phi_nd[lev]->nGrow() );

        if (proj_2)
        {
            // p := phi
            MultiFab::Copy (*p_g[lev], *phi_nd[lev], 0, 0, 1, phi_nd[lev]->nGrow());
            MultiFab::Copy ( *gp[lev], *fluxes[lev], 0, 0, 3, fluxes[lev]->nGrow());
        }
        else
        {
            // p := p + phi
            MultiFab::Add (*p_g[lev], *phi_nd[lev], 0, 0, 1, phi_nd[lev]->nGrow());
            MultiFab::Add ( *gp[lev], *fluxes[lev], 0, 0, 3, fluxes[lev]->nGrow());
        }
    }

    for (int lev = nlev-1; lev > 0; lev--)
    {
        avgDown (lev-1, *vel_g[lev], *vel_g[lev-1]);
        avgDown (lev-1, *   gp[lev],    *gp[lev-1]);
    }

    // Swap ghost cells and apply BCs to velocity
    mfix_set_velocity_bcs (time, vel_g, 0);

    mfix_compute_diveu(time);

    for (int lev = 0; lev < nlev; lev++)
    {
        amrex::Print() << "AT LEVEL " << lev << " AFTER PROJECTION: \n";
        mfix_print_max_vel(lev);
        mfix_print_max_gp (lev);
        amrex::Print() << "   max(abs(diveu)) = " <<  mfix_norm0(diveu, lev, 0) << "\n";
    }
}

void
mfix::mfix_solve_poisson_equation ( Vector< std::unique_ptr<MultiFab> >& this_phi,
                                    Vector< std::unique_ptr<MultiFab> >& rhs,
                                    Vector< std::unique_ptr<MultiFab> >& b,
                                    Vector< std::unique_ptr<MultiFab> >& fluxes)
{
    BL_PROFILE("mfix::mfix_solve_poisson_equation");
    //
    // When we created nodal_matrix we didn't call setSigma so the coefficients are initially 0
    // Here we set them to the values we will actually use
    //
    for (int lev = 0; lev < nlev; lev++)
       nodal_matrix->setSigma(lev, *b[lev]);

    // We must call setLevelBC regardless of the boundary conditions.
    // If Dirichlet then the values of this_phi in the ghost cells outside the domain
    //    boundary will be used; if Neumann then the default is homogeneous
    for (int lev = 0; lev < nlev; lev++)
       nodal_matrix->setLevelBC(lev, GetVecOfConstPtrs(this_phi)[lev]);

    if (nodal_bottom_solver_type == "smoother")
    {
       nodal_solver->setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if (nodal_bottom_solver_type == "bicg")
    {
       nodal_solver->setBottomSolver(MLMG::BottomSolver::bicgstab);
    }
    else if (nodal_bottom_solver_type == "cg")
    {
       nodal_solver->setBottomSolver(MLMG::BottomSolver::cg);
    }
    else if (nodal_bottom_solver_type == "bicgcg")
    {
       nodal_solver->setBottomSolver(MLMG::BottomSolver::bicgcg);
    }
    else if (nodal_bottom_solver_type == "cgbicg")
    {
       nodal_solver->setBottomSolver(MLMG::BottomSolver::cgbicg);
    }
    else if (nodal_bottom_solver_type == "hypre")
    {
       nodal_solver->setBottomSolver(MLMG::BottomSolver::hypre);
    }

    //
    // Solve div( eps_g/rho * grad(phi) ) = div(eps_g*u)
    //
    nodal_solver->solve    ( GetVecOfPtrs(this_phi), GetVecOfConstPtrs(rhs), nodal_mg_rtol, nodal_mg_atol);
    nodal_solver->getFluxes( GetVecOfPtrs(fluxes) );

    for (int lev = 0; lev < nlev; lev++)
      this_phi[lev] -> FillBoundary(geom[lev].periodicity());
}
