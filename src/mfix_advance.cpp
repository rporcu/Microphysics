#include <mfix_proj_F.H>
#include <mfix_F.H>
#include <mfix.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

void
mfix::EvolveFluid( int nstep, Real& dt,  Real& time, Real stop_time, Real coupling_timing )
{
    BL_PROFILE_REGION_START("mfix::EvolveFluid");
    BL_PROFILE("mfix::EvolveFluid");

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "EvolveFluid Start";
            MemProfiler::report(ss.str());
        }
#endif

    amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";

    // Extrapolate boundary values for ro_g, tracer, ep_g and mu_g
    // The subsequent call to mfix_set_scalar_bcs will only overwrite
    // ep_g ghost values for PINF and POUT
    for (int lev = 0; lev < nlev; lev++)
    {
        ro_g[lev]->FillBoundary(geom[lev].periodicity());
        trac[lev]->FillBoundary(geom[lev].periodicity());
        ep_g[lev]->FillBoundary(geom[lev].periodicity());
        mu_g[lev]->FillBoundary(geom[lev].periodicity());
    }

    // Fill ghost nodes and reimpose boundary conditions
    mfix_set_velocity_bcs (time, vel_g, 0);
    mfix_set_scalar_bcs   (time, ro_g, trac, ep_g, mu_g);

    //
    // Start loop: if we are not seeking a steady state solution,
    // the loop will execute only once
    //
    int keep_looping = 1;
    int iter = 1;

    // Create temporary multifabs to hold the old-time conv and divtau
    //    so we don't have to re-compute them in the corrector
    Vector<std::unique_ptr<MultiFab> > conv_u_old;
    Vector<std::unique_ptr<MultiFab> > conv_s_old;
    Vector<std::unique_ptr<MultiFab> > divtau_old;

    conv_u_old.resize(nlev);
    conv_s_old.resize(nlev);
    divtau_old.resize(nlev);

    for (int lev = 0; lev < nlev; lev++)
    {
       conv_u_old[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), *ebfactory[lev]));
       conv_s_old[lev].reset(new MultiFab(grids[lev], dmap[lev], 2, 0, MFInfo(), *ebfactory[lev]));
       divtau_old[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), *ebfactory[lev]));
    }

    do
    {
        mfix_compute_dt(nstep, time, stop_time, dt);

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
          MultiFab::Copy ( *ep_go[lev],  *ep_g[lev],  0, 0,  ep_g[lev]->nComp(),  ep_go[lev]->nGrow());
          MultiFab::Copy (  *p_go[lev],   *p_g[lev],  0, 0,   p_g[lev]->nComp(),   p_go[lev]->nGrow());
          MultiFab::Copy ( *ro_go[lev],  *ro_g[lev],  0, 0,  ro_g[lev]->nComp(),  ro_go[lev]->nGrow());
          MultiFab::Copy (*trac_o[lev],  *trac[lev],  0, 0,  trac[lev]->nComp(), trac_o[lev]->nGrow());
          MultiFab::Copy (*vel_go[lev], *vel_g[lev], 0, 0, vel_g[lev]->nComp(),  vel_go[lev]->nGrow());

           // User hooks
           for (MFIter mfi(*ep_g[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
              mfix_usr2();
        }

        //
        // Time integration step
        //
        Real new_time = time+dt;

        // Calculate drag coefficient
        if (solve_dem) {
          Real start_drag = ParallelDescriptor::second();
          mfix_calc_drag_fluid(time);
          coupling_timing += ParallelDescriptor::second() - start_drag;
        }

        // Predictor step
        bool proj_2_pred = true;
        mfix_apply_predictor ( conv_u_old, conv_s_old, divtau_old, time, dt, proj_2_pred );

        // Calculate drag coefficient
        if (solve_dem)
        {
          Real start_drag = ParallelDescriptor::second();
          amrex::Print() << "\nRecalculating drag ..." << std::endl;
          mfix_calc_drag_fluid(new_time);
          coupling_timing += ParallelDescriptor::second() - start_drag;
        }

        bool proj_2_corr = true;
        // Corrector step
        if (!steady_state)
           mfix_apply_corrector ( conv_u_old, conv_s_old, divtau_old, time, dt, proj_2_corr );

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

    if (test_tracer_conservation)
       amrex::Print() << "Sum tracer volume wgt = " << volWgtSum(0,*trac[0],0) << " " << volEpsWgtSum(0,*trac[0],0) << std::endl;

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "EvolveFluid Stop";
            MemProfiler::report(ss.str());
        }
#endif

    BL_PROFILE_REGION_STOP("mfix::EvolveFluid");
}

void
mfix::mfix_project_velocity ()
{
    // Project velocity field to make sure initial velocity is divergence-free
    Real dummy_dt = 1.0;

    amrex::Print() << "Initial projection:\n";

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
    int nstep = 0;
    mfix_compute_dt(nstep,time,stop_time,dt);

    amrex::Print() << "Doing initial pressure iterations with dt = " << dt << std::endl;

    // Fill ghost nodes and reimpose boundary conditions
    mfix_set_velocity_bcs (time, vel_g, 0);
    mfix_set_scalar_bcs   (time, ro_g, trac, ep_g, mu_g);

    // Copy vel_g into vel_go
    for (int lev = 0; lev < nlev; lev++)
       MultiFab::Copy (*vel_go[lev], *vel_g[lev],   0, 0, vel_g[lev]->nComp(), vel_go[lev]->nGrow());

    if (solve_dem)
       mfix_calc_drag_fluid(time);

    // Create temporary multifabs to hold conv and divtau
    Vector<std::unique_ptr<MultiFab> > conv_u;
    Vector<std::unique_ptr<MultiFab> > conv_s;
    Vector<std::unique_ptr<MultiFab> > divtau;

    conv_u.resize(nlev);
    conv_s.resize(nlev);
    divtau.resize(nlev);

    for (int lev = 0; lev < nlev; lev++)
    {
       conv_u[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), *ebfactory[lev]));
       conv_s[lev].reset(new MultiFab(grids[lev], dmap[lev], 2             , 0, MFInfo(), *ebfactory[lev]));
       divtau[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), *ebfactory[lev]));
    }

   for (int iter = 0; iter < initial_iterations; ++iter)
   {
       amrex::Print() << " " << std::endl;
       amrex::Print() << "In initial_iterations: iter = " << iter <<  "\n";

       bool proj_2 = false;

       mfix_apply_predictor (conv_u, conv_s, divtau, time, dt, proj_2);

       // Reset any quantities which might have been updated
       for (int lev = 0; lev < nlev; lev++)
          MultiFab::Copy (*vel_g[lev], *vel_go[lev], 0, 0, vel_g[lev]->nComp(), vel_g[lev]->nGrow());

       if (advect_density)
          for (int lev = 0; lev < nlev; lev++)
             MultiFab::Copy (*ro_g[lev], *ro_go[lev], 0, 0, ro_g[lev]->nComp(), ro_g[lev]->nGrow());

       if (advect_tracer)
          for (int lev = 0; lev < nlev; lev++)
             MultiFab::Copy (*trac[lev], *trac_o[lev], 0, 0, trac[lev]->nComp(), trac[lev]->nGrow());

       // Reset the boundary values (necessary if they are time-dependent)
       mfix_set_velocity_bcs (time, vel_g, 0);
       mfix_set_scalar_bcs   (time, ro_g, trac, ep_g, mu_g);
   }
}

//
// Compute predictor:
//
//  1. Compute
//
//     vel_g = vel_go + dt * R_u^n + dt * divtau*(1/(ro_g*ep_g))
//
//  2. Add explicit forcing term ( AKA gravity, lagged pressure gradient,
//     and explicit part of particles momentum exchange )
//
//     vel_g = vel_g + dt * ( g - grad(p_g+p0)/ro_g )
//
//  3. Add implicit forcing term ( AKA implicit part of particles
//     momentum exchange )
//
//     drag_coeff = drag(3)
//     drag_coeff*velp = drag(0:2)
//
//     vel_g = (vel_g + (drag_coeff*velp)/(ro_g*ep_g) / ( 1 + dt * drag_coeff/(ro_g*ep_g)
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
mfix::mfix_apply_predictor (Vector< std::unique_ptr<MultiFab> >& conv_u_old,
                            Vector< std::unique_ptr<MultiFab> >& conv_s_old,
                            Vector< std::unique_ptr<MultiFab> >& divtau_old,
                            Real time, Real dt, bool proj_2)
{
    // We use the new-time value for things computed on the "*" state
    Real new_time = time + dt;

    // Compute the explicit advective term R_u^n
    mfix_compute_convective_term( conv_u_old, conv_s_old, vel_go, ep_g, ro_go, trac_o ,time );

    // FOR NOW WE STILL DIVIDE BY EP_G BUT WE DON'T WANT TO KEEP DOING THIS!
    for (int lev = 0; lev < nlev; lev++)
       for (int i = 0; i < 3; i++)
          MultiFab::Divide  (*conv_u_old[lev],*ep_g[lev],0,i,1,0);

    for (int lev = 0; lev < nlev; lev++)
       for (int i = 0; i < 2; i++)
          MultiFab::Divide  (*conv_s_old[lev],*ep_g[lev],0,i,1,0);

    int explicit_diffusion_pred = 1;

    if (explicit_diffusion_pred == 1)
       mfix_compute_divtau( divtau_old, vel_go, time);
    else
       for (int lev = 0; lev < nlev; lev++)
          divtau_old[lev]->setVal(0.);

    for (int lev = 0; lev < nlev; lev++)
    {
        EB_set_covered(*divtau_old[lev], 0, divtau_old[lev]->nComp(), divtau_old[lev]->nGrow(), 0.0);

        // First add the convective term
        MultiFab::Saxpy (*vel_g[lev], dt, *conv_u_old[lev], 0, 0, 3, 0);

        // Make sure to do this multiply before we update density!
        if (advect_tracer)
        {
           int conv_comp = 1;
           MultiFab::Multiply(*trac[lev],*ro_go[lev],0,0,1,0);
           MultiFab::Saxpy   (*trac[lev], dt, *conv_s_old[lev], conv_comp, 0, 1, 0);
        }

        if (advect_density)
        {
           int conv_comp = 0;
           MultiFab::Saxpy (*ro_g[lev], dt, *conv_s_old[lev], conv_comp, 0, 1, 0);
        }

        // Make sure to do this divide after we update density!
        if (advect_tracer)
           MultiFab::Divide  (*trac[lev],*ro_g[lev],0,0,1,0);

        // Add the explicit diffusion terms
        if (explicit_diffusion_pred == 1)
           MultiFab::Saxpy (*vel_g[lev], dt, *divtau_old[lev], 0, 0, 3, 0);
    }

    // Add source terms
    mfix_add_gravity_and_gp(dt);

    // Add the drag term implicitly
    if (solve_dem)
        mfix_add_drag_implicit (dt);

    // If doing implicit diffusion, solve here for u^*
    if (explicit_diffusion_pred == 0)
        mfix_diffuse_velocity_tensor(new_time,dt);

    // Project velocity field
    mfix_apply_projection ( new_time, dt, proj_2 );

    mfix_set_velocity_bcs (new_time, vel_g, 0);
}

//
// Compute corrector:
//
//  1. Compute
//
//     vel_g = vel_go + dt * (R_u^* + R_u^n) / 2 + dt * divtau*(1/(ro_g*ep_g))
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
//     vel_g = (vel_g + (drag_coeff*velp)/(ro_g*ep_g) / ( 1 + dt * drag_coeff/(ro_g*ep_g)
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
mfix::mfix_apply_corrector (Vector< std::unique_ptr<MultiFab> >& conv_u_old,
                            Vector< std::unique_ptr<MultiFab> >& conv_s_old,
                            Vector< std::unique_ptr<MultiFab> >& divtau_old,
                            Real time, Real dt, bool proj_2)
{
    BL_PROFILE("mfix::mfix_apply_corrector");

    // We use the new-time value for things computed on the "*" state
    Real new_time = time + dt;

    // Create temporary multifabs to hold the new-time conv and divtau
    Vector<std::unique_ptr<MultiFab> > conv_u;
    Vector<std::unique_ptr<MultiFab> > conv_s;
    Vector<std::unique_ptr<MultiFab> > divtau;

    conv_u.resize(nlev);
    conv_s.resize(nlev);
    divtau.resize(nlev);

    for (int lev = 0; lev < nlev; lev++)
    {
       conv_u[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), *ebfactory[lev]));
       conv_s[lev].reset(new MultiFab(grids[lev], dmap[lev], 2             , 0, MFInfo(), *ebfactory[lev]));
       divtau[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), *ebfactory[lev]));
    }

    // Compute the explicit advective term R_u^*
    mfix_compute_convective_term( conv_u, conv_s, vel_g, ep_g, ro_g, trac, new_time );

    // FOR NOW WE STILL DIVIDE BY EP_G BUT WE DON'T WANT TO KEEP DOING THIS!
    for (int lev = 0; lev < nlev; lev++)
       for (int i = 0; i < 3; i++)
          MultiFab::Divide  (*conv_u[lev],*ep_g[lev],0,i,1,0);

    for (int lev = 0; lev < nlev; lev++)
    {
        for (int i = 0; i < 2; i++)
          MultiFab::Divide  (*conv_s[lev],*ep_g[lev],0,i,1,0);

        // Make sure to do this multiply before we update density!
        if (advect_tracer)
        {
          int conv_comp = 1;
          MultiFab::Multiply(*trac_o[lev],*ro_go[lev],0,0,1,0);

          // Add the convective terms so trac = trac_o + dt/2 (R_s^* + R_s^n)
          MultiFab::LinComb (*trac[lev], 1.0, *trac_o[lev], 0, dt/2.0, *conv_s[lev]    , conv_comp, 0, 1, 0);
          MultiFab::Saxpy   (*trac[lev],                       dt/2.0, *conv_s_old[lev], conv_comp, 0, 1, 0);
        }

        if (advect_density)
        {
           int conv_comp = 0;
           // Add the convective terms so trac = trac_o + dt/2 (R_s^* + R_s^n)
           MultiFab::LinComb (*ro_g[lev], 1.0, *ro_go[lev], 0, dt/2.0, *conv_s[lev]    , conv_comp, 0, 1, 0);
           MultiFab::Saxpy   (*ro_g[lev],                      dt/2.0, *conv_s_old[lev], conv_comp, 0, 1, 0);
        }

        // Make sure to do this divide after we update density!
        if (advect_tracer)
           MultiFab::Divide  (*trac[lev],*ro_g[lev],0,0,1,0);

        // Add the convective terms so u_g = u_go + dt/2 (R_u^* + R_u^n)
        MultiFab::LinComb (*vel_g[lev], 1.0, *vel_go[lev], 0, dt/2.0, *conv_u[lev]    , 0, 0, 3, 0);
        MultiFab::Saxpy   (*vel_g[lev],                       dt/2.0, *conv_u_old[lev], 0, 0, 3, 0);
    }

    // Add the explicit diffusion term so u_g = u_go + dt/2 (R_u^* + R_u^n) + dt/2 (Lu)^n
    for (int lev = 0; lev < nlev; lev++)
        MultiFab::Saxpy (*vel_g[lev], dt/2.0, *divtau_old[lev], 0, 0, 3, 0);

    // Add source terms
    mfix_add_gravity_and_gp(dt);

    // Add the drag term implicitly
    if (solve_dem)
        mfix_add_drag_implicit(dt);

    // Solve for u^star s.t. u^star = u_go + dt/2 (R_u^* + R_u^n) + dt/2 (Lu)^n + dt/2 (Lu)^star
    mfix_diffuse_velocity_tensor(new_time,.5*dt);

    // Apply projection
    mfix_apply_projection (new_time, dt, proj_2);

    mfix_set_velocity_bcs (new_time, vel_g, 0);
}

void
mfix::mfix_add_gravity_and_gp (Real dt)
{
    BL_PROFILE("mfix::mfix_add_gravity_and_gp");
    for (int lev = 0; lev < nlev; lev++)
    {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(*vel_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
         // Tilebox
         Box bx = mfi.tilebox ();

         const auto& vel_fab = vel_g[lev]->array(mfi);
         const auto&  gp_fab =    gp[lev]->array(mfi);
         const auto& den_fab =  ro_g[lev]->array(mfi);

         const auto grav_loc = gravity;
         const auto  gp0_loc = gp0;

         AMREX_FOR_3D(bx, i, j, k,
         {
             Real inv_dens = 1.0 / den_fab(i,j,k);
             vel_fab(i,j,k,0) += dt * ( grav_loc[0]-(gp_fab(i,j,k,0)+gp0_loc[0])*inv_dens );
             vel_fab(i,j,k,1) += dt * ( grav_loc[1]-(gp_fab(i,j,k,1)+gp0_loc[1])*inv_dens );
             vel_fab(i,j,k,2) += dt * ( grav_loc[2]-(gp_fab(i,j,k,2)+gp0_loc[2])*inv_dens );
         });

         // NOTE: here we do not need host-device synchronization since it is
         // already included in the MFIter destructor
       }
    }
}

//
// Explicit solve for the intermediate velocity.
// Currently this means accounting for the implicit part of the fluid/particle
// momentum exchange
//
void
mfix::mfix_add_drag_explicit (Real dt)
{
  /*
     This adds both components of the drag term
     So the drag term we add is beta * (particle_velocity - fluid_velocity)
                              = drag(0:2) - drag(3) * fluid_velocity
  */

  BL_PROFILE("mfix::mfix_add_drag_explicit");

  for (int lev = 0; lev < nlev; lev++)
  {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vel_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      // Tilebox
      Box bx = mfi.tilebox ();

      const auto&  vel_fab = vel_g[lev]->array(mfi);
      const auto& drag_fab =  drag[lev]->array(mfi);
      const auto&   ro_fab =  ro_g[lev]->array(mfi);
      const auto&   ep_fab =  ep_g[lev]->array(mfi);

      AMREX_FOR_3D(bx, i, j, k,
      {
          Real orop  = dt / (ro_fab(i,j,k) * ep_fab(i,j,k));

          Real drag_0 = (drag_fab(i,j,k,0) - drag_fab(i,j,k,3)*vel_fab(i,j,k,0)) * orop;
          Real drag_1 = (drag_fab(i,j,k,1) - drag_fab(i,j,k,3)*vel_fab(i,j,k,1)) * orop;
          Real drag_2 = (drag_fab(i,j,k,2) - drag_fab(i,j,k,3)*vel_fab(i,j,k,2)) * orop;

          vel_fab(i,j,k,0) += drag_0;
          vel_fab(i,j,k,1) += drag_1;
          vel_fab(i,j,k,2) += drag_2;
      });

      // NOTE: here we do not need host-device synchronization since it is
      // already included in the MFIter destructor
    }
  }
}

//
// Implicit solve for the intermediate velocity.
// Currently this means accounting for the implicit part of the fluid/particle
// momentum exchange
//
void
mfix::mfix_add_drag_implicit (Real dt)
{
  /*
     This adds both components of the drag term
     So the drag term we add is beta * (particle_velocity - fluid_velocity)
                              = drag(0:2) - drag(3) * fluid_velocity
  */

  BL_PROFILE("mfix::mfix_add_drag_implicit");

  for (int lev = 0; lev < nlev; lev++)
  {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vel_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      // Tilebox
      Box bx = mfi.tilebox ();

      const auto&  vel_fab = vel_g[lev]->array(mfi);
      const auto& drag_fab =  drag[lev]->array(mfi);
      const auto&   ro_fab =  ro_g[lev]->array(mfi);
      const auto&   ep_fab =  ep_g[lev]->array(mfi);

      AMREX_FOR_3D(bx, i, j, k,
      {
          Real orop  = dt / (ro_fab(i,j,k) * ep_fab(i,j,k));
          Real denom = 1.0 / (1.0 + drag_fab(i,j,k,3) * orop);

          vel_fab(i,j,k,0) = (vel_fab(i,j,k,0) + drag_fab(i,j,k,0) * orop) * denom;
          vel_fab(i,j,k,1) = (vel_fab(i,j,k,1) + drag_fab(i,j,k,1) * orop) * denom;
          vel_fab(i,j,k,2) = (vel_fab(i,j,k,2) + drag_fab(i,j,k,2) * orop) * denom;
      });
 
      // NOTE: here we do not need host-device synchronization since it is
      // already included in the MFIter destructor
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
    mfix_set_velocity_bcs (time, vel_g, 0);

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

    reached = reached || (iter >= steady_state_maxiter);

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
