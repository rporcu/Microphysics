#include <AMReX_ParmParse.H>
#include <mfix_mac_F.H>
#include <mfix_proj_F.H>
#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>

void
mfix_level::EvolveFluidProjection(int lev, int nstep, int steady_state, Real& dt,  Real& time, Real stop_time )
{
    BL_PROFILE_REGION_START("mfix::EvolveFluidProjection");
    BL_PROFILE("mfix::EvolveFluidProjection");

    amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";
    // Extrapolate boundary values for density and volume fraction
    // The subsequent call to mfix_set_scalar_bcs will only overwrite
    // rop_g and ep_g ghost values for PINF and POUT
    fill_mf_bc ( lev, *rop_g[lev] );
    fill_mf_bc ( lev, *ep_g[lev] );
    fill_mf_bc ( lev, *mu_g[lev] );

    // Fill ghost nodes and reimpose boundary conditions
    mfix_set_scalar_bcs (lev);
    mfix_set_velocity_bcs (lev,0);
    
    //
    // Start loop: if we are not seeking a steady state solution,
    // the loop will execute only once
    //
    int keep_looping = 1;
    int iter = 1;

    do
    {
        mfix_compute_dt(lev, time, stop_time, steady_state, dt);

        if (steady_state)
        {
           amrex::Print() << "\n   Iteration " << iter << " with dt = " << dt << "\n" << std::endl;
        } else {
           amrex::Print() << "\n   Step " << nstep+1 << ": from old_time " \
                          << time << " to new time " << time+dt
                          << " with dt = " << dt << "\n" << std::endl;
        }
 
        // Back up field
        // Backup field variable to old
	MultiFab::Copy (*ep_go[lev],  *ep_g[lev],  0, 0,  ep_g[lev]->nComp(),  ep_go[lev]->nGrow());
	MultiFab::Copy ( *p_go[lev],   *p_g[lev],  0, 0,   p_g[lev]->nComp(),   p_go[lev]->nGrow());
	MultiFab::Copy (*ro_go[lev],  *ro_g[lev],  0, 0,  ro_g[lev]->nComp(),  ro_go[lev]->nGrow());
	MultiFab::Copy (*rop_go[lev], *rop_g[lev], 0, 0, rop_g[lev]->nComp(), rop_go[lev]->nGrow());
	MultiFab::Copy (*vel_go[lev], *vel_g[lev], 0, 0, vel_g[lev]->nComp(), vel_go[lev]->nGrow());

        // User hooks
        for (MFIter mfi(*ep_g[lev], true); mfi.isValid(); ++mfi)
    	   mfix_usr2();

	//
	// Time integration step
	//
	
	// Calculate drag coefficient
	if (solve_dem)
	    mfix_calc_drag_fluid(lev);

	// Create temporary multifabs to hold the old-time conv and divtau
	//    so we don't have to re-compute them in the corrector
        MultiFab   conv_old(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
        MultiFab divtau_old(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
	
        // Predictor step 
        bool proj_2 = true;
        mfix_apply_predictor ( lev, conv_old, divtau_old, dt, proj_2 );

	// Print info about predictor step
        {
	    amrex::Print() << "\nAfter predictor step:\n";
    	    mfix_print_max_vel (lev);
    	    mfix_compute_diveu (lev);
	    amrex::Print() << "max(abs(diveu)) = " << mfix_norm0(diveu, lev, 0) << "\n";
        }

	// Calculate drag coefficient
	if (solve_dem)
	    mfix_calc_drag_fluid(lev);
	
	// Corrector step 
        proj_2 = true;
	mfix_apply_corrector ( lev, conv_old, divtau_old, dt, proj_2 );

	// Print info about corrector step
        {
	    amrex::Print() << "\nAfter corrector step:\n";
    	    mfix_print_max_vel (lev);
    	    mfix_compute_diveu (lev);
	    amrex::Print() << "max(abs(diveu)) = " << mfix_norm0(diveu, lev, 0) << "\n";
        }
	    
	// 
        // Check whether to exit the loop or not
	// 
	if (steady_state) {
	    keep_looping = !steady_state_reached ( lev, dt ); 
        } else {
	    keep_looping = 0;
	}

	// Update interations count
	++iter;
    }
    while ( keep_looping );

    BL_PROFILE_REGION_STOP("mfix::EvolveFluidProjection");
}

void
mfix_level::mfix_project_velocity (int lev)
{
    // Project velocity field to make sure initial velocity is divergence-free
    Real dummy_dt = 1.0;

    amrex::Print() << "Initial projection:\n";

    // Need to add this call here so that the MACProjection internal arrays
    //  are allocated so that the cell-centered projection can use the MAC
    //  data structures and set_velocity_bcs routine
    mac_projection->update_internals();

    bool proj_2 = true;
    mfix_apply_projection ( lev, dummy_dt, proj_2 );

   // We initialize p_g back to zero (p0_g may still be still non-zero)
   p_g[lev]->setVal(0.0);
}

void
mfix_level::mfix_initial_iterations (int lev, Real dt, Real stop_time, int steady_state)
{
    // Fill ghost cells
    mfix_set_scalar_bcs (lev);
    mfix_set_velocity_bcs (lev,0);

    // Copy vel_g into vel_go
    MultiFab::Copy (*vel_go[lev], *vel_g[lev],   0, 0, vel_g[lev]->nComp(), vel_go[lev]->nGrow());

    Real time = 0.0;
    mfix_compute_dt(lev,time,stop_time,steady_state,dt);

    amrex::Print() << "Doing initial pressure iterations with dt = " << dt << std::endl;

    //  Create temporary multifabs to hold conv and divtau
    MultiFab   conv(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
    MultiFab divtau(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);

    for (int iter = 0; iter < 3; ++iter)
    {
       amrex::Print() << "In initial_iterations: iter = " << iter <<  "\n";

       bool proj_2 = false;
       mfix_apply_predictor (lev, conv, divtau, dt, proj_2);

       // Replace vel_g by the original values 
       MultiFab::Copy (*vel_g[lev], *vel_go[lev], 0, 0, vel_g[lev]->nComp(), vel_g[lev]->nGrow());

       //exit(0);
    }
}

//
// Compute predictor:
//
//  1. Compute
//
//     vel_g = vel_go + dt * R_u^n + dt * divtau*(1/rop_g)
//
//  2. Add explicit forcing term ( AKA gravity, lagged pressure gradient,
//     and explicit part of particles momentum exchange )
//
//     vel_g = vel_g + dt * ( g - grad(p_g+p0)/ro_g + drag_u/rop_g )
//
//  3. Add implicit forcing term ( AKA implicit part of particles
//     momentum exchange )
//
//     vel_g = vel_g / ( 1 + dt * f_gds/rop_g )
//
//  4. Solve for phi
//
//     div( ep_g * grad(phi) / ro_g ) = div( ep_g * vel_g / dt + grad(p_g)/ro_g )
//
//  5. Compute
//
//     vel_g = vel_g -  dt * grad(phi) / ro_g
//
//  6. Define
//
//     p_g = phi
//
void
mfix_level::mfix_apply_predictor (int lev, MultiFab& conv_old, MultiFab& divtau_old,
                                  amrex::Real dt, bool proj_2)
{
    // Compute the explicit advective term R_u^n
    mfix_compute_ugradu_predictor( lev, conv_old, vel_go );

    // If explicit_diffusion == true  then we compute the full diffusive terms here
    // If explicit_diffusion == false then we compute only the off-diagonal terms here
    mfix_compute_divtau( lev, divtau_old, vel_go );
    
    // First add the convective term
    MultiFab::Saxpy (*vel_g[lev], dt,   conv_old, 0, 0, 3, 0);

    // Add the diffusion terms (either all if explicit_diffusion == true or just the
    //    off-diagonal terms if explicit_diffusion == false)
    MultiFab::Saxpy (*vel_g[lev], dt, divtau_old, 0, 0, 3, 0);

    // Add the forcing terms
    mfix_apply_forcing_terms ( lev, dt, vel_g);

    // Convert velocities to momenta
    for (int n = 0; n < 3; n++)
       MultiFab::Multiply(*vel_g[lev],(*ro_g[lev]),0,n,1,vel_g[lev]->nGrow());

    // Add (-dt grad p to momenta)
    MultiFab::Saxpy (*vel_g[lev], -dt,  *gp[lev], 0, 0, 3, vel_g[lev]->nGrow());
    MultiFab::Saxpy (*vel_g[lev], -dt, *gp0[lev], 0, 0, 3, vel_g[lev]->nGrow());

    // Convert momenta back to velocities
    for (int n = 0; n < 3; n++)
       MultiFab::Divide(*vel_g[lev],(*ro_g[lev]),0,n,1,vel_g[lev]->nGrow());

    // If doing implicit diffusion, solve here for u^*
    if (!explicit_diffusion)
       mfix_diffuse_velocity(lev,dt);

    // Add the drag term implicitly
    if (solve_dem)
       mfix_compute_intermediate_velocity ( lev, dt );
 
    // Project velocity field
    mfix_apply_projection ( lev, dt, proj_2 );
}

//
// Compute corrector:
//
//  1. Compute
//
//     vel_g = vel_go + dt * (R_u^* + R_u^n) / 2 + dt * divtau*(1/rop_g)
//
//     where the starred variables are computed using "predictor-step" variables.
//
//  2. Add explicit forcing term ( AKA gravity, lagged pressure gradient,
//     and explicit part of particles momentum exchange )
//
//     vel_g = vel_g + dt * ( g - grad(p_g+p0)/ro_g + drag_u/rop_g )
//
//  3. Add implicit forcing term ( AKA implicit part of particles
//     momentum exchange )
//
//     vel_g = vel_g / ( 1 + dt * f_gds/rop_g )
//
//  4. Solve for phi
//
//     div( ep_g * grad(phi) / ro_g ) = div( ep_g * vel_g / dt + grad(p_g)/ro_g )
//
//  5. Compute
//
//     vel_g = vel_g -  dt * grad(phi) / ro_g
//
//  6. Define
//
//     p_g = phi
//
void
mfix_level::mfix_apply_corrector (int lev, MultiFab& conv_old, MultiFab& divtau_old,
                                  amrex::Real dt, bool proj_2)
{
    BL_PROFILE("mfix_level::mfix_apply_corrector");

    MultiFab   conv(grids[lev], dmap[lev], 3, 0 );
    MultiFab divtau(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);

    // Compute the explicit advective term R_u^*
    mfix_compute_ugradu_corrector( lev,   conv, vel_g );

    // If explicit_diffusion == true  then we compute the full diffusive terms here
    // If explicit_diffusion == false then we compute only the off-diagonal terms here
    mfix_compute_divtau ( lev, divtau, vel_g);
        
    // Define u_g = u_go + dt/2 (R_u^* + R_u^n) 
    MultiFab::LinComb (*vel_g[lev], 1.0, *vel_go[lev], 0, dt/2.0, conv    , 0, 0, 3, 0); 
    MultiFab::Saxpy   (*vel_g[lev],                       dt/2.0, conv_old, 0, 0, 3, 0);

    // Add the diffusion terms (either all if explicit_diffusion == true or just the
    //    off-diagonal terms if explicit_diffusion == false)
    MultiFab::Saxpy (*vel_g[lev], dt/2.0, divtau    , 0, 0, 3, 0);
    MultiFab::Saxpy (*vel_g[lev], dt/2.0, divtau_old, 0, 0, 3, 0);

    // Add forcing terms
    mfix_apply_forcing_terms ( lev, dt, vel_g);

    // Convert velocities to momenta
    for (int n = 0; n < 3; n++)
       MultiFab::Multiply(*vel_g[lev],(*ro_g[lev]),0,n,1,vel_g[lev]->nGrow());

    // Add (-dt grad p to momenta)
    MultiFab::Saxpy (*vel_g[lev], -dt,  *gp[lev], 0, 0, 3, vel_g[lev]->nGrow());
    MultiFab::Saxpy (*vel_g[lev], -dt, *gp0[lev], 0, 0, 3, vel_g[lev]->nGrow());

    // Convert momenta back to velocities
    for (int n = 0; n < 3; n++)
       MultiFab::Divide(*vel_g[lev],(*ro_g[lev]),0,n,1,vel_g[lev]->nGrow());

    // If doing implicit diffusion, solve here for u^*
    if (!explicit_diffusion)
       mfix_diffuse_velocity(lev,dt);
 
    // Compute intermediate velocity if drag terms present
    if (solve_dem)
       mfix_compute_intermediate_velocity ( lev, dt );
 
    // Apply projection
    mfix_apply_projection ( lev, dt, proj_2 );
}

void
mfix_level::mfix_add_grad_phi (int lev, amrex::Real coeff, MultiFab& this_phi)
{
    BL_PROFILE("mfix_level::mfix_add_grad_phi");
    
    if (nodal_pressure == 1)
    {
#ifdef _OPENMP
#pragma omp parallel 
#endif
       for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
       {
   	   // Tilebox
   	   Box bx = mfi.tilebox();

	   add_grad_phind (
	       BL_TO_FORTRAN_BOX(bx),  
	       BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
	       BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
	       BL_TO_FORTRAN_ANYD(this_phi[mfi]),
	       geom[lev].CellSize(), &coeff);
       }
    } else {
#ifdef _OPENMP
#pragma omp parallel 
#endif
       for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
       {
   	   // Tilebox
   	   Box bx = mfi.tilebox();

	   add_grad_phicc (
	       BL_TO_FORTRAN_BOX(bx),  
	       BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
	       BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
	       this_phi[mfi].dataPtr(),
	       geom[lev].CellSize(), &coeff);
       }
    }
    
}


void
mfix_level::mfix_apply_forcing_terms (int lev, amrex::Real dt,
				      Vector< std::unique_ptr<MultiFab> >& vel) 

{
    BL_PROFILE("mfix_level::mfix_apply_forcing_terms");

    Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
    {
	// Tilebox
	Box bx = mfi.tilebox ();

	add_forcing ( BL_TO_FORTRAN_BOX(bx),  
		      BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*drag[lev])[mfi]),
		      BL_TO_FORTRAN_ANYD((*ro_g[lev])[mfi]),
                      (*rop_g[lev])[mfi].dataPtr(),
		      domain.loVect (), domain.hiVect (),
		      geom[lev].CellSize (), &dt);
    }
}

//
// Implicit solve for the intermediate velocity.
// Currently this means accounting for the implicit part of the fluid/particle
// momentum exchange
// 
void
mfix_level::mfix_compute_intermediate_velocity ( int lev, amrex::Real dt )

{
    BL_PROFILE("mfix_level::mfix_compute_intermediate_velocity");

    // Whole domain
    Box domain(geom[lev].Domain());
    
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi) {
	
	// Tilebox
	Box bx = mfi.tilebox();

	compute_intermediate_velocity ( BL_TO_FORTRAN_BOX(bx),  
					BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*f_gds[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
					&dt );
    }
}


//
// Compute div(ep_g * u)
// 
void
mfix_level::mfix_compute_diveu (int lev)
{
    Box domain(geom[lev].Domain());

    if (nodal_pressure == 1)
    {

       // Create a temporary multifab to hold (ep_g * vel_g)
       MultiFab vec(vel_g[lev]->boxArray(), vel_g[lev]->DistributionMap(), vel_g[lev]->nComp(), vel_g[lev]->nGrow());

       // Fill it with (ep_g * vel_g)
       vec.copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),vel_g[lev]->nGrow(),vel_g[lev]->nGrow());
       for (int n = 0; n < 3; n++)
          MultiFab::Multiply(vec,(*ep_g[lev]),0,n,1,1);

#ifdef _OPENMP
#pragma omp parallel
#endif
      // Extrapolate Dirichlet values to ghost cells -- but do it differently in that 
      //  no-slip walls are treated exactly like slip walls -- this is only relevant
      //  when going into the projection
      for (MFIter mfi(vec, true); mfi.isValid(); ++mfi)
        {
          set_vec_bcs ( BL_TO_FORTRAN_ANYD(vec[mfi]),
                        bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                        bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                        bc_klo.dataPtr(), bc_khi.dataPtr(),
                        domain.loVect(), domain.hiVect(),
                        &nghost);
        }

        vec.FillBoundary (geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIter mfi(*diveu[lev],true); mfi.isValid(); ++mfi)
       {
           // Note: this box is nodal!
	   const Box& bx = mfi.tilebox();
	
	   compute_diveund ( BL_TO_FORTRAN_BOX(bx),
			     BL_TO_FORTRAN_ANYD((*diveu[lev])[mfi]),
			     BL_TO_FORTRAN_ANYD(vec[mfi]),
	                     geom[lev].CellSize());
       }

    }
    else
    {

       int extrap_dir_bcs = 1;
       mfix_set_velocity_bcs (lev,extrap_dir_bcs);
       vel_g[lev]->FillBoundary (geom[lev].periodicity());     

       // Create face centered multifabs for ep_g and vel_g
       MultiFab epu( vel_g[lev]->boxArray(),  vel_g[lev]-> DistributionMap(),
                     vel_g[lev]->nComp(), nghost, MFInfo(), *ebfactory[lev]);

       MultiFab::Copy( epu, *vel_g[lev], 0, 0, 3, vel_g[lev]->nGrow() );

       for (int n = 0; n < 3; n++)
          MultiFab::Multiply( epu, *ep_g[lev], 0, n, 1, vel_g[lev]->nGrow() );

       Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> epu_fc;
       mfix_average_cc_to_fc( lev, epu, epu_fc );

       // This does not need to have correct ghost values in place
       EB_computeDivergence( *diveu[lev], GetArrOfConstPtrs(epu_fc), geom[lev] );

    }

    // Restore velocities to carry Dirichlet values on faces
    int extrap_dir_bcs = 0;
    mfix_set_velocity_bcs (lev,extrap_dir_bcs);
}

//
// Check if steady state has been reached by verifying that
// 
//      max(abs( u^(n+1) - u^(n) )) < tol * dt
//      max(abs( v^(n+1) - v^(n) )) < tol * dt
//      max(abs( w^(n+1) - w^(n) )) < tol * dt
// 

int
mfix_level::steady_state_reached (int lev, Real dt)
{
    //
    // Count number of access 
    //
    static int naccess = 0;

    //
    // Make sure velocity is up to date
    // 
    mfix_set_velocity_bcs (lev,0);
    
    // 
    // Use temporaries to store the difference
    // between current and previous solution
    // 
    MultiFab temp_vel(vel_g[lev]->boxArray(), vel_g[lev]->DistributionMap(),3,0);
    MultiFab::LinComb (temp_vel, 1.0, *vel_g[lev], 0, -1.0, *vel_go[lev], 0, 0, 3, 0);

    MultiFab tmp;
    
    if (nodal_pressure)
    {
       const BoxArray & nd_grid = amrex::convert(grids[lev], IntVect{1,1,1});
       tmp.define( nd_grid, dmap[lev], 1, 0 );     
    }
    else
    {
	tmp.define( grids[lev], dmap[lev], 1, 0 );
    }

    MultiFab::LinComb (tmp, 1.0, *p_g[lev], 0, -1.0, *p_go[lev], 0, 0, 1, 0);
    
    Real delta_u = mfix_norm0(temp_vel,lev,0);
    Real delta_v = mfix_norm0(temp_vel,lev,1);
    Real delta_w = mfix_norm0(temp_vel,lev,2);
    Real delta_p = mfix_norm0(tmp,lev,0);
    
    Real tol = steady_state_tol; 
    
    int condition1 = (delta_u < tol*dt) && (delta_v < tol*dt ) && (delta_w < tol*dt);

    //
    // Second stop condition
    //
    Real du_n1 = mfix_norm1(temp_vel, lev, 0); 
    Real dv_n1 = mfix_norm1(temp_vel, lev, 1); 
    Real dw_n1 = mfix_norm1(temp_vel, lev, 2); 
    Real dp_n1 = mfix_norm1(tmp, lev, 0); 
    Real uo_n1 = mfix_norm1(vel_go, lev, 0);
    Real vo_n1 = mfix_norm1(vel_go, lev, 1);
    Real wo_n1 = mfix_norm1(vel_go, lev, 2);
    Real po_n1 = mfix_norm1(p_go, lev, 0);
    
    Real tmp1, tmp2, tmp3, tmp4;

    Real local_tol = 1.0e-8;
    
    if ( uo_n1 < local_tol ) {
    	tmp1 = 0.0;
    } else {
    	tmp1 = du_n1 / uo_n1;
    };

    if ( vo_n1 < local_tol ) {
    	tmp2 = 0.0;
    } else {
    	tmp2 = dv_n1 / vo_n1;
    };
    
    if ( wo_n1 < local_tol ) {
    	tmp3 = 0.0;
    } else {
    	tmp3 = dw_n1 / wo_n1;
    };

    if ( po_n1 < local_tol ) {
    	tmp4 = 0.0;
    } else {
    	tmp4 = dp_n1 / po_n1;
    };

    int condition2 = (tmp1 < tol) && (tmp2 < tol) && (tmp3 < tol); // && (tmp4 < tol);

    //
    // Print out info on steady state checks
    //
    amrex::Print() << "\nSteady state check:\n";
    amrex::Print() << "||u-uo||/||uo|| , du/dt  = " << tmp1 <<" , "<< delta_u/dt << "\n";
    amrex::Print() << "||v-vo||/||vo|| , dv/dt  = " << tmp2 <<" , "<< delta_v/dt << "\n";
    amrex::Print() << "||w-wo||/||wo|| , dw/dt  = " << tmp3 <<" , "<< delta_w/dt << "\n";
    amrex::Print() << "||p-po||/||po|| , dp/dt  = " << tmp4 <<" , "<< delta_p/dt << "\n";

    // Count # access
    naccess++;

    // 
    //  Always return negative to first access. This way
    //  initial zero velocity field do not test for false positive
    //
    if ( naccess == 1 ) {
	return 0;
    } else {
	return condition1 || condition2;
    };
}

//
// Set the BCs for all the variables EXCEPT pressure or velocity.
//
void
mfix_level::mfix_set_scalar_bcs (int lev)
{
  BL_PROFILE("mfix_level::mfix_set_scalar_bcs()");

  Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(*ep_g[lev], true); mfi.isValid(); ++mfi)
    {
      set_scalar_bcs ( BL_TO_FORTRAN_ANYD((*ep_g[lev])[mfi]),
                      (*ro_g[lev])[mfi].dataPtr (),
                      (*rop_g[lev])[mfi].dataPtr (),
                      (*mu_g[lev])[mfi].dataPtr (),
                      (*lambda_g[lev])[mfi].dataPtr (),
                      bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                      bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                      bc_klo.dataPtr(), bc_khi.dataPtr(),
                      domain.loVect(), domain.hiVect(),
                      &nghost );
    }
        ep_g[lev] -> FillBoundary (geom[lev].periodicity());
        ro_g[lev] -> FillBoundary (geom[lev].periodicity());
       rop_g[lev] -> FillBoundary (geom[lev].periodicity());
        mu_g[lev] -> FillBoundary (geom[lev].periodicity());
    lambda_g[lev] -> FillBoundary (geom[lev].periodicity());
}

//
// Set the BCs for velocity only
//
void
mfix_level::mfix_set_velocity_bcs (int lev, int extrap_dir_bcs)
{
  BL_PROFILE("mfix_level::mfix_set_velocity_bcs()");

  vel_g[lev] -> FillBoundary (geom[lev].periodicity());

  Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(*vel_g[lev], true); mfi.isValid(); ++mfi)
    {
      set_velocity_bcs ( BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
                         bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                         bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                         bc_klo.dataPtr(), bc_khi.dataPtr(),
                         domain.loVect(), domain.hiVect(),
                         &nghost, &extrap_dir_bcs );
    }
}

//
// Fills ghost cell values of pressure appropriately for the BC type
//
void
mfix_level::mfix_extrap_pressure (int lev, std::unique_ptr<amrex::MultiFab>& p)
{
    BL_PROFILE("mfix_level::mfix_extrap_pressure()");
    if (nodal_pressure == 1) return;
 
    Box domain(geom[lev].Domain());
 
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for (MFIter mfi(*p, true); mfi.isValid(); ++mfi) {
 
        extrap_pressure_to_ghost_cells (
            BL_TO_FORTRAN_ANYD((*p)[mfi]),
            bc_ilo.dataPtr(), bc_ihi.dataPtr(),
            bc_jlo.dataPtr(), bc_jhi.dataPtr(),
            bc_klo.dataPtr(), bc_khi.dataPtr(),
            domain.loVect(), domain.hiVect(),
            &nghost);
    }
}

void
mfix_level::check_for_nans (int lev)
{
    bool ug_has_nans = vel_g[lev] -> contains_nan (0);
    bool vg_has_nans = vel_g[lev] -> contains_nan (1);
    bool wg_has_nans = vel_g[lev] -> contains_nan (2);
    bool pg_has_nans =   p_g[lev] -> contains_nan (0);
    bool ropg_has_nans = rop_g[lev] -> contains_nan (0);

    if (ug_has_nans)
	amrex::Print() << "WARNING: u_g contains NaNs!!!";

    if (vg_has_nans)
	amrex::Print() << "WARNING: v_g contains NaNs!!!";

    if (wg_has_nans)
	amrex::Print() << "WARNING: w_g contains NaNs!!!";

    if (pg_has_nans)
	amrex::Print() << "WARNING: p_g contains NaNs!!!";

    if (ropg_has_nans)
	amrex::Print() << "WARNING: rop_g contains NaNs!!!";

}

//
// Print the maximum values of the velocity components
//
void
mfix_level::mfix_print_max_vel(int lev)
{
    amrex::Print() << "max(abs(u/v/w/p))  = " << 
       mfix_norm0(vel_g, lev, 0) << "  " <<
       mfix_norm0(vel_g, lev, 1) << "  " <<
       mfix_norm0(vel_g, lev, 2) << "  " <<
       mfix_norm0(p_g,   lev, 0) << "  " << std::endl;
}



//
// This subroutines averages component by component
// The assumption is that cc is multicomponent
// 
void
mfix_level::mfix_average_cc_to_fc ( int lev, const MultiFab& cc,
                                    Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& fc )
{
   AMREX_ASSERT(cc.nComp()==AMREX_SPACEDIM);
   AMREX_ASSERT(AMREX_SPACEDIM==3);
   
   // 
   // First allocate fc
   //
   BoxArray x_ba = cc.boxArray();
   x_ba.surroundingNodes(0);
   fc[0].reset(new MultiFab(x_ba,cc.DistributionMap(),1,nghost, MFInfo(), *ebfactory[lev]));
   fc[0]->setVal(1.e200);

   BoxArray y_ba = cc.boxArray();
   y_ba.surroundingNodes(1);
   fc[1].reset(new MultiFab(y_ba,cc.DistributionMap(),1,nghost, MFInfo(), *ebfactory[lev]));
   fc[1]->setVal(1.e200);

   BoxArray z_ba = cc.boxArray();
   z_ba.surroundingNodes(2);
   fc[2].reset(new MultiFab(z_ba,cc.DistributionMap(),1,nghost, MFInfo(), *ebfactory[lev]));
   fc[2]->setVal(1.e200);

   //
   // Average
   // We do not care about EB because faces in covered regions
   // should never get used so we can set them to whatever values
   // we like
   //
#ifdef _OPENMP
#pragma omp parallel 
#endif
   for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
   {
      // Boxes for staggered components
      Box bx = mfi.tilebox();

      
      average_cc_to_fc( BL_TO_FORTRAN_BOX(bx),
                        BL_TO_FORTRAN_ANYD((*fc[0])[mfi]),
                        BL_TO_FORTRAN_ANYD((*fc[1])[mfi]),
                        BL_TO_FORTRAN_ANYD((*fc[2])[mfi]),
                        BL_TO_FORTRAN_ANYD(cc[mfi]));

   }

   // We do not fill BCs and halo regions in this routine    
} 





