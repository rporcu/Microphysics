#include <mfix_mac_F.H>
#include <mfix_proj_F.H>
#include <mfix_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>

void
mfix::EvolveFluid( int nstep, Real& dt,  Real& time, Real stop_time )
{
    BL_PROFILE_REGION_START("mfix::EvolveFluid");
    BL_PROFILE("mfix::EvolveFluid");

    amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";

    // Extrapolate boundary values for density and volume fraction
    // The subsequent call to mfix_set_scalar_bcs will only overwrite
    // rop_g and ep_g ghost values for PINF and POUT
    for (int lev = 0; lev < nlev; lev++)
    {
       fill_mf_bc ( lev, *rop_g[lev] );
       fill_mf_bc ( lev, *ep_g[lev] );
       fill_mf_bc ( lev, *mu_g[lev] );
    }

    // Fill ghost nodes and reimpose boundary conditions
    mfix_set_velocity_bcs (time, 0);
    mfix_set_scalar_bcs ();

    //
    // Start loop: if we are not seeking a steady state solution,
    // the loop will execute only once
    //
    int keep_looping = 1;
    int iter = 1;

    // Create temporary multifabs to hold the old-time conv and divtau
    //    so we don't have to re-compute them in the corrector
    Vector<std::unique_ptr<MultiFab> > conv_old;
    Vector<std::unique_ptr<MultiFab> > divtau_old;

      conv_old.resize(nlev);
    divtau_old.resize(nlev);

    for (int lev = 0; lev < nlev; lev++)
    {
         conv_old[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
       divtau_old[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
    }

    do
    {
        mfix_compute_dt(time, stop_time, dt);

        // Set new and old time to correctly use in fillpatching
        for (int lev = 0; lev < nlev; lev++)
        {
            t_old[lev] = time;
            t_new[lev] = time+dt;
        }

        if (steady_state)
        {
           amrex::Print() << "\n   Iteration " << iter << " with dt = " << dt << "\n" << std::endl;
        } else {
           amrex::Print() << "\n   Step " << nstep+1 << ": from old_time " \
                          << time << " to new time " << time+dt
                          << " with dt = " << dt << "\n" << std::endl;
        }

        for (int lev = 0; lev < nlev; lev++)
        {
           // Back up field variables to old
          MultiFab::Copy (*ep_go[lev],  *ep_g[lev],  0, 0,  ep_g[lev]->nComp(),  ep_go[lev]->nGrow());
          MultiFab::Copy ( *p_go[lev],   *p_g[lev],  0, 0,   p_g[lev]->nComp(),   p_go[lev]->nGrow());
          MultiFab::Copy (*ro_go[lev],  *ro_g[lev],  0, 0,  ro_g[lev]->nComp(),  ro_go[lev]->nGrow());
          MultiFab::Copy (*rop_go[lev], *rop_g[lev], 0, 0, rop_g[lev]->nComp(), rop_go[lev]->nGrow());
          MultiFab::Copy (*vel_go[lev], *vel_g[lev], 0, 0, vel_g[lev]->nComp(), vel_go[lev]->nGrow());

           // User hooks
           for (MFIter mfi(*ep_g[lev], true); mfi.isValid(); ++mfi)
              mfix_usr2();
        }

        //
        // Time integration step
        //
        Real new_time = time+dt;

        // Calculate drag coefficient
        if (solve_dem)
          mfix_calc_drag_fluid(time);

        // Predictor step
        bool proj_2_pred = true;
        mfix_apply_predictor ( conv_old, divtau_old, time, dt, proj_2_pred );

        // Print info about predictor step
        amrex::Print() << "\nAfter predictor step at time " << new_time << std::endl;
        for (int lev = 0; lev < nlev; lev++)
          mfix_print_max_vel (lev);

        mfix_compute_diveu(new_time);

        for (int lev = 0; lev < nlev; lev++)
          amrex::Print() << "max(abs(diveu)) = " << mfix_norm0(diveu, lev, 0) << "\n";

        // Calculate drag coefficient
        if (solve_dem)
          mfix_calc_drag_fluid(new_time);

        bool proj_2_corr = true;
        // Corrector step
        if (!steady_state)
          mfix_apply_corrector ( conv_old, divtau_old, time, dt, proj_2_corr );

        // Print info about corrector step
        amrex::Print() << "\nAfter corrector step at time " << new_time << std::endl;
        for (int lev = 0; lev < nlev; lev++)
          mfix_print_max_vel (lev);

        mfix_compute_diveu (new_time);

        for (int lev = 0; lev < nlev; lev++)
          amrex::Print() << "  max(abs(diveu)) = " << mfix_norm0(diveu, lev, 0) << "\n";

        //
        // Check whether to exit the loop or not
        //
        if (steady_state) {
          keep_looping = !steady_state_reached ( dt, iter);
        } else {
          keep_looping = 0;
        }

        // Update interations count
        ++iter;
    }
    while ( keep_looping );

    BL_PROFILE_REGION_STOP("mfix::EvolveFluid");
}

void
mfix::mfix_project_velocity ()
{
    // Project velocity field to make sure initial velocity is divergence-free
    Real dummy_dt = 1.0;

    amrex::Print() << "Initial projection:\n";

    // Need to add this call here so that the MACProjection internal arrays
    //  are allocated so that the cell-centered projection can use the MAC
    //  data structures and set_velocity_bcs routine
    mac_projection->update_internals();

    bool proj_2 = true;
    Real time = 0.0;
    mfix_apply_projection ( time, dummy_dt, proj_2 );

   // We initialize p_g and gp back to zero (p0_g may still be still non-zero)
   for (int lev = 0; lev < nlev; lev++)
   {
      p_g[lev]->setVal(0.0);
       gp[lev]->setVal(0.0);
   }
}

void
mfix::mfix_initial_iterations (Real dt, Real stop_time)
{

   Real time = 0.0;
   mfix_compute_dt(time,stop_time,dt);

   amrex::Print() << "Doing initial pressure iterations with dt = " << dt << std::endl;

   // Fill ghost cells
   mfix_set_scalar_bcs ();
   mfix_set_velocity_bcs (time, 0);

   // Copy vel_g into vel_go
   for (int lev = 0; lev < nlev; lev++)
      MultiFab::Copy (*vel_go[lev], *vel_g[lev],   0, 0, vel_g[lev]->nComp(), vel_go[lev]->nGrow());

    // Create temporary multifabs to hold conv and divtau
    Vector<std::unique_ptr<MultiFab> > conv;
    Vector<std::unique_ptr<MultiFab> > divtau;

      conv.resize(nlev);
    divtau.resize(nlev);

    for (int lev = 0; lev < nlev; lev++)
    {
         conv[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
       divtau[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
    }

   for (int iter = 0; iter < initial_iterations; ++iter)
   {
       amrex::Print() << " " << std::endl;
       amrex::Print() << "In initial_iterations: iter = " << iter <<  "\n";

       bool proj_2 = false;

       mfix_apply_predictor (conv, divtau, time, dt, proj_2);

       // Replace vel_g by the original values
       for (int lev = 0; lev < nlev; lev++)
          MultiFab::Copy (*vel_g[lev], *vel_go[lev], 0, 0, vel_g[lev]->nComp(), vel_g[lev]->nGrow());

       // Reset the boundary values (necessary if they are time-dependent)
       mfix_set_velocity_bcs (time, 0);
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
//     vel_g = vel_g + dt * ( g - grad(p_g+p0)/ro_g )
//
//  3. Add implicit forcing term ( AKA implicit part of particles
//     momentum exchange )
//
//     vel_g = (vel_g + drag_u/rop_g) / ( 1 + dt * f_gds/rop_g )
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
mfix::mfix_apply_predictor (Vector< std::unique_ptr<MultiFab> >& conv_old,
                      Vector< std::unique_ptr<MultiFab> >& divtau_old,
                            Real time, Real dt, bool proj_2)
{
    // We use the new-time value for things computed on the "*" state
    Real new_time = time + dt;

    // Compute the explicit advective term R_u^n
    mfix_compute_ugradu_predictor( conv_old, vel_go, time );

    for (int lev = 0; lev < nlev; lev++)
    {
        // If explicit_diffusion == true  then we compute the full diffusive terms here
        // If explicit_diffusion == false then we compute only the off-diagonal terms here
        mfix_compute_divtau( lev, *divtau_old[lev], vel_go );

        // First add the convective term
        MultiFab::Saxpy (*vel_g[lev], dt, *conv_old[lev], 0, 0, 3, 0);

        // Add the diffusion terms (either all if explicit_diffusion == true or just the
        //    off-diagonal terms if explicit_diffusion == false)
        MultiFab::Saxpy (*vel_g[lev], dt, *divtau_old[lev], 0, 0, 3, 0);

        // Add the gravitational forcing
        mfix_add_gravity_and_gp ( lev, dt);
    }

    // Add the drag term implicitly
    if (solve_dem)
        mfix_add_drag_terms ( dt );

    // If doing implicit diffusion, solve here for u^*
    if (!explicit_diffusion)
        mfix_diffuse_velocity(new_time,dt);

    // Project velocity field
    mfix_apply_projection ( new_time, dt, proj_2 );

    mfix_set_velocity_bcs (new_time, 0);
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
//     vel_g = vel_g + dt * ( g - grad(p_g+p0)/ro_g )
//
//  3. Add implicit forcing term ( AKA implicit part of particles
//     momentum exchange )
//
//     vel_g = (vel_g + drag_u/rop_g) / ( 1 + dt * f_gds/rop_g )
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
mfix::mfix_apply_corrector (Vector< std::unique_ptr<MultiFab> >& conv_old,
                      Vector< std::unique_ptr<MultiFab> >& divtau_old,
                            Real time, Real dt, bool proj_2)
{
    BL_PROFILE("mfix::mfix_apply_corrector");

    // We use the new-time value for things computed on the "*" state
    Real new_time = time + dt;

    // Create temporary multifabs to hold the old-time conv and divtau
    //    so we don't have to re-compute them in the corrector
    Vector<std::unique_ptr<MultiFab> > conv;
    Vector<std::unique_ptr<MultiFab> > divtau;

      conv.resize(nlev);
    divtau.resize(nlev);

    for (int lev = 0; lev < nlev; lev++)
    {
         conv[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
       divtau[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
    }

    // Compute the explicit advective term R_u^*
    mfix_compute_ugradu_corrector( conv, vel_g, new_time );

    for (int lev = 0; lev < nlev; lev++)
    {
        // If explicit_diffusion == true  then we compute the full diffusive terms here
        // If explicit_diffusion == false then we compute only the off-diagonal terms here
        mfix_compute_divtau( lev, *divtau[lev], vel_g);

        // Define u_g = u_go + dt/2 (R_u^* + R_u^n)
        MultiFab::LinComb (*vel_g[lev], 1.0, *vel_go[lev], 0, dt/2.0, *conv[lev]    , 0, 0, 3, 0);
        MultiFab::Saxpy   (*vel_g[lev],                       dt/2.0, *conv_old[lev], 0, 0, 3, 0);

        // Add the diffusion terms (either all if explicit_diffusion == true or just the
        //    off-diagonal terms if explicit_diffusion == false)
        MultiFab::Saxpy (*vel_g[lev], dt/2.0, *divtau[lev]    , 0, 0, 3, 0);
        MultiFab::Saxpy (*vel_g[lev], dt/2.0, *divtau_old[lev], 0, 0, 3, 0);

        // Add the gravitational forcing
        mfix_add_gravity_and_gp ( lev, dt);
    }

    // Compute intermediate velocity if drag terms present
    if (solve_dem)
        mfix_add_drag_terms (dt);

    // If doing implicit diffusion, solve here for u^*
    if (!explicit_diffusion)
       mfix_diffuse_velocity(new_time,dt);

    // Apply projection
    mfix_apply_projection (new_time, dt, proj_2);

    mfix_set_velocity_bcs (new_time, 0);
}

void
mfix::mfix_add_gravity_and_gp (int lev, amrex::Real dt) 
{
    BL_PROFILE("mfix::mfix_add_gravity_and_gp");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
    {
      // Tilebox
      Box bx = mfi.tilebox ();

      const auto& vel_fab = vel_g[lev]->array(mfi);

      for (int n = 0; n < 3; n++)
       for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); k++)
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++)
          for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++)
          {
             vel_fab(i,j,k,n) += dt * gravity[n];
          }
    }

    for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
    {
      // Grown tilebox
      Box bx = mfi.growntilebox ();

      const auto& vel_fab = vel_g[lev]->array(mfi);
      const auto&  gp_fab =    gp[lev]->array(mfi);
      const auto& den_fab =  ro_g[lev]->array(mfi);

      for (int n = 0; n < 3; n++)
       for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); k++)
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++)
          for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++)
          {
             vel_fab(i,j,k,n) -= dt * gp_fab(i,j,k,n)/den_fab(i,j,k);
          }

      for (int n = 0; n < 3; n++)
       for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); k++)
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++)
          for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++)
          {
             vel_fab(i,j,k,n) -= dt * gp00[n]/den_fab(i,j,k);
          }
    }
}

//
// Implicit solve for the intermediate velocity.
// Currently this means accounting for the implicit part of the fluid/particle
// momentum exchange
//
void
mfix::mfix_add_drag_terms ( amrex::Real dt )

{
  BL_PROFILE("mfix::mfix_add_drag");

  for (int lev = 0; lev < nlev; lev++)
    {
      // The volume fraction of each fluid cell (1 if uncovered, 0 if covered)
      const amrex::MultiFab* volfrac = &(ebfactory[lev] -> getVolFrac());

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
        {
          // Tilebox
          Box bx = mfi.tilebox();

          add_drag_terms ( BL_TO_FORTRAN_BOX(bx),
                           BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
                           BL_TO_FORTRAN_ANYD((*f_gds[lev])[mfi]),
                           BL_TO_FORTRAN_ANYD((*drag[lev])[mfi]),
                           BL_TO_FORTRAN_ANYD((*rop_g[lev])[mfi]),
                           BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                           &dt );
        }
    }
}

//
// Check if steady state has been reached by verifying that
//
//      max(abs( u^(n+1) - u^(n) )) < tol * dt
//      max(abs( v^(n+1) - v^(n) )) < tol * dt
//      max(abs( w^(n+1) - w^(n) )) < tol * dt
//

int
mfix::steady_state_reached (Real dt, int iter)
{
    //
    // Count number of access
    //
    static int naccess = 0;

    int condition1[nlev];
    int condition2[nlev];

    Real time = 0.;
    mfix_set_velocity_bcs (time, 0);

    //
    // Make sure velocity is up to date
    //
    for (int lev = 0; lev < nlev; lev++)
    {

       //
       // Use temporaries to store the difference
       // between current and previous solution
       //
       MultiFab temp_vel(vel_g[lev]->boxArray(), vel_g[lev]->DistributionMap(),3,0);
       MultiFab::LinComb (temp_vel, 1.0, *vel_g[lev], 0, -1.0, *vel_go[lev], 0, 0, 3, 0);

       MultiFab tmp;

       const BoxArray & nd_grid = amrex::convert(grids[lev], IntVect{1,1,1});
       tmp.define( nd_grid, dmap[lev], 1, 0 );

       MultiFab::LinComb (tmp, 1.0, *p_g[lev], 0, -1.0, *p_go[lev], 0, 0, 1, 0);

       Real delta_u = mfix_norm0(temp_vel,lev,0);
       Real delta_v = mfix_norm0(temp_vel,lev,1);
       Real delta_w = mfix_norm0(temp_vel,lev,2);
       Real delta_p = mfix_norm0(tmp,lev,0);

       Real tol = steady_state_tol;

       condition1[lev] = (delta_u < tol*dt) && (delta_v < tol*dt ) && (delta_w < tol*dt);

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

       condition2[lev] = (tmp1 < tol) && (tmp2 < tol) && (tmp3 < tol); // && (tmp4 < tol);

       //
       // Print out info on steady state checks
       //
       amrex::Print() << "\nSteady state check at level " << lev << ":\n";
       amrex::Print() << "||u-uo||/||uo|| , du/dt  = " << tmp1 <<" , "<< delta_u/dt << "\n";
       amrex::Print() << "||v-vo||/||vo|| , dv/dt  = " << tmp2 <<" , "<< delta_v/dt << "\n";
       amrex::Print() << "||w-wo||/||wo|| , dw/dt  = " << tmp3 <<" , "<< delta_w/dt << "\n";
       amrex::Print() << "||p-po||/||po|| , dp/dt  = " << tmp4 <<" , "<< delta_p/dt << "\n";
    }

    int reached = 1;
    for (int lev = 0; lev < nlev; lev++)
    {
       reached = reached && (condition1[lev] || condition2[lev]);
    }

    reached = reached || (iter >= steady_state_max_iter);

    // Count # access
    naccess++;

    //
    //  Always return negative to first access. This way
    //  initial zero velocity field do not test for false positive
    //
    if ( naccess == 1 ) {
       return 0;
    } else {
       return reached;
    };
}
