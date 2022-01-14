#include <mfix.H>

#include <AMReX_VisMF.H>
#include <mfix_mf_helpers.H>
#include <mfix_dem_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_eb_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_solvers.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif


using namespace Solvers;


void
mfix::mfix_project_velocity ()
{
    // Project velocity field to make sure initial velocity is divergence-free
    Real dummy_dt = 1.0;

    amrex::Print() << "Initial projection:\n";

    bool proj_2 = true;
    Real time = 0.0;

    // Apply projection -- depdt=0 for now
    Vector< MultiFab* > depdt(finest_level+1);
    for (int lev(0); lev <= finest_level; ++lev)
      depdt[lev] = MFHelpers::createFrom(*(m_leveldata[lev]->ep_g), 0.0, 1).release();

    mfix_apply_nodal_projection(depdt, time, dummy_dt, dummy_dt, proj_2,
                                get_vel_g_old(), get_vel_g(), get_p_g(), get_gp(),
                                get_ep_g(), get_txfr(), get_ro_g_const());

    for (int lev(0); lev <= finest_level; ++lev)
      delete depdt[lev];

    // We initialize p_g and gp back to zero (p0_g may still be still non-zero)
    for (int lev = 0; lev <= finest_level; lev++) {
      m_leveldata[lev]->p_g->setVal(0);
      m_leveldata[lev]->gp->setVal(0);
    }

    // Call the initial redistribution after the initial projection.
    InitialRedistribution(time);
}

void
mfix::mfix_initial_iterations (Real dt, Real stop_time)
{
  Real time = 0.0;
  int nstep = 0;

  mfix_compute_dt(nstep, time, stop_time, dt, dt);

  amrex::Print() << "Doing initial pressure iterations with dt = " << dt << "\n";

  // Fill ghost nodes and reimpose boundary conditions
  mfix_set_velocity_bcs(time, get_vel_g(), 0);
  mfix_set_density_bcs(time, get_ro_g());
  mfix_set_tracer_bcs(time, get_trac());

  if (advect_enthalpy) {
    mfix_set_temperature_bcs(time, get_T_g());
    mfix_set_enthalpy_bcs(time, get_h_g());
  }

  if (advect_enthalpy && EB::fix_temperature)
    mfix_set_eb_temperature_bcs(get_T_g_on_eb());

  if (solve_species)
    mfix_set_species_bcs(time, get_X_gk());

  // Copy vel_g and p_g into vel_go and p_go
  for (int lev = 0; lev <= finest_level; lev++) {
    MultiFab::Copy(*m_leveldata[lev]->vel_go, *m_leveldata[lev]->vel_g, 0, 0,
                   m_leveldata[lev]->vel_g->nComp(), m_leveldata[lev]->vel_g->nGrow());

    MultiFab::Copy(*m_leveldata[lev]->p_go, *m_leveldata[lev]->p_g, 0, 0,
                   m_leveldata[lev]->p_g->nComp(), m_leveldata[lev]->p_g->nGrow());
  }

  if (DEM::solve || PIC::solve) {
    mfix_calc_txfr_fluid(get_txfr(), get_ep_g(), get_ro_g(), get_vel_g(),
                         get_T_g(), get_X_gk(), time);

    if (reactions.solve) {
      mfix_calc_chem_txfr(get_chem_txfr(), get_ep_g(), get_ro_g(), get_vel_g(),
                          get_p_g(), get_T_g(), get_X_gk(), time);
    }
  }

  // Create temporary multifabs to hold conv and vel_RHS
  Vector< MultiFab* > conv_u(finest_level+1, nullptr);
  Vector< MultiFab* > conv_s(finest_level+1, nullptr);
  Vector< MultiFab* > conv_X(finest_level+1, nullptr);
  Vector< MultiFab* > ro_RHS_old(finest_level+1, nullptr);
  Vector< MultiFab* > lap_trac_old(finest_level+1, nullptr);
  Vector< MultiFab* > enthalpy_RHS_old(finest_level+1, nullptr);
  Vector< MultiFab* > lap_T_old(finest_level+1, nullptr);
  Vector< MultiFab* > species_RHS_old(finest_level+1, nullptr);
  Vector< MultiFab* > vel_RHS_old(finest_level+1, nullptr);
  Vector< MultiFab* > div_J_old(finest_level+1, nullptr);
  Vector< MultiFab* > div_hJ_old(finest_level+1, nullptr);
  Vector< Real > rhs_pressure_g_old(finest_level+1, 0.);
  Vector< Real > rhs_pressure_g(finest_level+1, 0.);

  for (int lev = 0; lev <= finest_level; lev++)
  {
    conv_u[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
    conv_u[lev]->setVal(0.0);

    // 2+ntrac: one for density, ntrac for tracers and one for enthalpy
    conv_s[lev] = new MultiFab(grids[lev], dmap[lev], 2+ntrac, 0, MFInfo(), *ebfactory[lev]);
    conv_s[lev]->setVal(0.0);

    if (advect_density) {
      ro_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
      ro_RHS_old[lev]->setVal(0.0);
    }

    if (advect_tracer) {
      lap_trac_old[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), *ebfactory[lev]);
      lap_trac_old[lev]->setVal(0.0);
    }

    if (advect_enthalpy) {
      enthalpy_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
      enthalpy_RHS_old[lev]->setVal(0.0);
      lap_T_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
      lap_T_old[lev]->setVal(0.0);
    }

    if (solve_species) {
      div_J_old[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies, 0, MFInfo(), *ebfactory[lev]);
      div_J_old[lev]->setVal(0.0);
    }

    if (advect_enthalpy && solve_species) {
      div_hJ_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
      div_hJ_old[lev]->setVal(0.0);
    }

    if (solve_species) {
      conv_X[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies, 0, MFInfo(), *ebfactory[lev]);
      conv_X[lev]->setVal(0.0);

      species_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies, 0, MFInfo(), *ebfactory[lev]);
      species_RHS_old[lev]->setVal(0.0);
    }

    if (reactions.solve) {
      vel_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
    }
  }

  for (int iter = 0; iter < initial_iterations; ++iter)
  {
    amrex::Print() << " " << std::endl;
    amrex::Print() << "In initial_iterations: iter = " << iter << "\n";

    bool proj_2 = false;

    auto dt_copy = dt;

    Real coupling_timing(0);

    mfix_apply_predictor(conv_u, conv_s, conv_X, ro_RHS_old, lap_trac_old,
        enthalpy_RHS_old, lap_T_old, species_RHS_old, vel_RHS_old,
        div_J_old, div_hJ_old, rhs_pressure_g_old, rhs_pressure_g, time, dt,
        dt_copy, proj_2, coupling_timing);

    // Reset any quantities which might have been updated
    for (int lev = 0; lev <= finest_level; lev++) {
      MultiFab::Copy(*m_leveldata[lev]->vel_g, *m_leveldata[lev]->vel_go, 0, 0,
                     m_leveldata[lev]->vel_g->nComp(), m_leveldata[lev]->vel_g->nGrow());

      if (advect_density)
        MultiFab::Copy(*m_leveldata[lev]->ro_g, *m_leveldata[lev]->ro_go, 0, 0,
                       m_leveldata[lev]->ro_g->nComp(), m_leveldata[lev]->ro_g->nGrow());

      if (advect_enthalpy) {
        MultiFab::Copy(*m_leveldata[lev]->T_g, *m_leveldata[lev]->T_go, 0, 0,
                       m_leveldata[lev]->T_g->nComp(), m_leveldata[lev]->T_g->nGrow());
        MultiFab::Copy(*m_leveldata[lev]->h_g, *m_leveldata[lev]->h_go, 0, 0,
                       m_leveldata[lev]->h_g->nComp(), m_leveldata[lev]->h_g->nGrow());
      }

      if (advect_tracer)
        MultiFab::Copy(*m_leveldata[lev]->trac, *m_leveldata[lev]->trac_o, 0, 0,
                       m_leveldata[lev]->trac->nComp(), m_leveldata[lev]->trac->nGrow());

      if (solve_species) {
        MultiFab::Copy(*m_leveldata[lev]->X_gk, *m_leveldata[lev]->X_gko, 0, 0,
                       m_leveldata[lev]->X_gk->nComp(), m_leveldata[lev]->X_gk->nGrow());
      }

//      if ((m_constraint_type == ConstraintType::IdealGasOpenSystem ||
//           m_constraint_type == ConstraintType::IdealGasClosedSystem) && advect_enthalpy) {
//        MultiFab::Copy(*m_leveldata[lev]->pressure_g, *m_leveldata[lev]->pressure_go, 0, 0,
//                        m_leveldata[lev]->pressure_g->nComp(), m_leveldata[lev]->pressure_g->nGrow());
//      }
    }

    // Reset the boundary values (necessary if they are time-dependent)
    mfix_set_velocity_bcs(time, get_vel_g(), 0);
    mfix_set_density_bcs(time, get_ro_g());
    mfix_set_tracer_bcs(time, get_trac());

    if (advect_enthalpy)
      mfix_set_temperature_bcs(time, get_T_g());

    if (advect_enthalpy)
      mfix_set_enthalpy_bcs(time, get_h_g());

    if (solve_species)
      mfix_set_species_bcs(time, get_X_gk());
  }

  for (int lev = 0; lev <= finest_level; lev++)
  {
     delete conv_u[lev];
     delete conv_s[lev];

     if (advect_density)
       delete ro_RHS_old[lev];

     if (advect_tracer)
       delete lap_trac_old[lev];

     if (advect_enthalpy) {
       delete enthalpy_RHS_old[lev];
       delete lap_T_old[lev];
     }

     if (advect_enthalpy && solve_species) {
       delete div_hJ_old[lev];
     }

     if (solve_species) {
       delete conv_X[lev];
       delete species_RHS_old[lev];
       delete div_J_old[lev];
     }

     if (reactions.solve)
       delete vel_RHS_old[lev];
  }
}

//
// Explicit solve for the intermediate velocity.
// Currently this means accounting for the explicit part of the fluid/particle
// momentum exchange
//
void
mfix::mfix_add_txfr_explicit (Real dt,
                              Vector<MultiFab*      > const& vel_in,
                              Vector<MultiFab*      > const& h_g_in,
                              Vector<MultiFab*      > const& T_g_in,
                              Vector<MultiFab const*> const& X_gk_in,
                              Vector<MultiFab const*> const& txfr_in,
                              Vector<MultiFab const*> const& rho_in,
                              Vector<MultiFab const*> const& ep_g_in)
{
  /*
     This adds both components of the drag term
     So the drag term we add is beta * (particle_velocity - fluid_velocity)
                              = drag(0:2) - drag(3) * fluid_velocity
  */

  BL_PROFILE("mfix::mfix_add_txfr_explicit");

  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

  auto& fluid_parms = *fluid.parameters;

  for (int lev = 0; lev <= finest_level; lev++) {

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ep_g_in[lev]->Factory());
    const auto& flags = factory.getMultiEBCellFlagFab();
    const auto& volfrac = factory.getVolFrac();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      // Tilebox
      Box bx = mfi.tilebox();

      Array4<Real      > const&  vel_array = vel_in[lev]->array(mfi);
      Array4<Real const> const& txfr_array = txfr_in[lev]->array(mfi);
      Array4<Real const> const&   ro_array = rho_in[lev]->array(mfi);
      Array4<Real const> const&   ep_array = ep_g_in[lev]->array(mfi);

      amrex::ParallelFor(bx,[dt,vel_array,txfr_array,ro_array,ep_array]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const Real orop  = dt / (ro_array(i,j,k) * ep_array(i,j,k));

        const Real beta = txfr_array(i,j,k,Transfer::beta);

        const Real vel_x = vel_array(i,j,k,0);
        const Real vel_y = vel_array(i,j,k,1);
        const Real vel_z = vel_array(i,j,k,2);

        const Real drag_0 = (txfr_array(i,j,k,Transfer::velx) - beta*vel_x) * orop;
        const Real drag_1 = (txfr_array(i,j,k,Transfer::vely) - beta*vel_y) * orop;
        const Real drag_2 = (txfr_array(i,j,k,Transfer::velz) - beta*vel_z) * orop;

        vel_array(i,j,k,0) = vel_x + drag_0;
        vel_array(i,j,k,1) = vel_y + drag_1;
        vel_array(i,j,k,2) = vel_z + drag_2;
      });

      if(advect_enthalpy) {

        Array4<Real      > const& hg_array  = h_g_in[lev]->array(mfi);
        Array4<Real      > const& Tg_array  = T_g_in[lev]->array(mfi);
        Array4<Real const> const& Xgk_array = X_gk_in[lev]->const_array(mfi);

        auto const& flags_arr = flags.const_array(mfi);
        auto const& volfrac_arr = volfrac.const_array(mfi);

        const int nspecies_g = fluid.nspecies;
        const int fluid_is_a_mixture = fluid.is_a_mixture;

        amrex::ParallelFor(bx,[dt,hg_array,Tg_array,txfr_array,ro_array,ep_array,
            fluid_parms,Xgk_array,nspecies_g,fluid_is_a_mixture,flags_arr,
            volfrac_arr,run_on_device]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

          if (!cell_is_covered) {
            const Real epg_loc = ep_array(i,j,k);
            const Real vfrac   = volfrac_arr(i,j,k);

            const Real orop  = dt / (ro_array(i,j,k) * epg_loc);

            const Real Tg_old = Tg_array(i,j,k);
            const Real Ts     = txfr_array(i,j,k,Transfer::gammaTp);
            const Real gamma  = txfr_array(i,j,k,Transfer::gamma);

            const Real hg = hg_array(i,j,k) + (Ts - gamma * Tg_old) * orop;
            hg_array(i,j,k) = hg;

            // ************************************************************
            // Newton-Raphson solver for solving implicit equation for
            // temperature
            // ************************************************************
            // Residual computation
            auto R = [&] AMREX_GPU_DEVICE (Real Tg_arg)
            {
              Real hg_loc(0);

              if (!fluid_is_a_mixture) {

                hg_loc = run_on_device ?
                  fluid_parms.calc_h_g<RunOn::Device>(Tg_arg, cell_is_covered) :
                  fluid_parms.calc_h_g<RunOn::Host>(Tg_arg, cell_is_covered);

              } else {

                for (int n(0); n < nspecies_g; ++n) {
                  const Real h_gk = run_on_device ?
                    fluid_parms.calc_h_gk<RunOn::Device>(Tg_arg, n, cell_is_covered) :
                    fluid_parms.calc_h_gk<RunOn::Host>(Tg_arg, n, cell_is_covered);

                  hg_loc += Xgk_array(i,j,k,n)*h_gk;
                }
              }

              return hg_loc - hg;
            };

            // Partial derivative computation
            auto partial_R = [&] AMREX_GPU_DEVICE (Real Tg_arg)
            {
              Real gradient(0);

              if (!fluid_is_a_mixture) {

                gradient = run_on_device ?
                  fluid_parms.calc_partial_h_g<RunOn::Device>(Tg_arg) :
                  fluid_parms.calc_partial_h_g<RunOn::Host>(Tg_arg);

              } else {

                for (int n(0); n < nspecies_g; ++n) {
                  const Real h_gk = run_on_device ?
                    fluid_parms.calc_partial_h_gk<RunOn::Device>(Tg_arg,n) :
                    fluid_parms.calc_partial_h_gk<RunOn::Host>(Tg_arg,n);

                  gradient += Xgk_array(i,j,k,n)*h_gk;
                }
              }

              return gradient;
            };

            Real Tg_new(Tg_old);

            int solver_iterations(0);

            {
              DampedNewton::DampingFactor damping_factor(0., 0.);
              solver_iterations = 
                DampedNewton::solve(Tg_new, R, partial_R, damping_factor(epg_loc, vfrac),
                                    1.e-8, 1.e-8, 500);

            } if (solver_iterations == 500) {

              DampedNewton::DampingFactor damping_factor(1., 0.);
              solver_iterations =
                DampedNewton::solve(Tg_new, R, partial_R, damping_factor(epg_loc, vfrac),
                                    1.e-7, 1.e-7, 500);

            } if (solver_iterations == 500) {

              DampedNewton::DampingFactor damping_factor(1., 1.);
              solver_iterations =
                DampedNewton::solve(Tg_new, R, partial_R, damping_factor(epg_loc, vfrac),
                                    1.e-6, 1.e-6, 500);

            } if (solver_iterations == 500) {
              amrex::Abort("Damped-Newton solver did not converge");
            }


            Tg_array(i,j,k) = Tg_new;
          }
        });
      }
    }
  }
}

//
// Implicit solve for the intermediate velocity.
// Currently this means accounting for the implicit part of the fluid/particle
// momentum exchange
//
void
mfix::mfix_add_txfr_implicit (Real dt,
                              Vector<MultiFab*      > const& vel_in,
                              Vector<MultiFab*      > const& h_g_in,
                              Vector<MultiFab*      > const& T_g_in,
                              Vector<MultiFab const*> const& X_gk_in,
                              Vector<MultiFab const*> const& txfr_in,
                              Vector<MultiFab const*> const& rho_in,
                              Vector<MultiFab const*> const& ep_g_in)
{
  /*
     This adds both components of the drag term
     So the drag term we add is beta * (particle_velocity - fluid_velocity)
                              = drag(0:2) - drag(3) * fluid_velocity
  */

  BL_PROFILE("mfix::mfix_add_txfr_implicit");

  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

  auto& fluid_parms = *fluid.parameters;

  for (int lev = 0; lev <= finest_level; lev++) {
    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ep_g_in[lev]->Factory());
    const auto& flags = factory.getMultiEBCellFlagFab();

    const auto& volfrac = factory.getVolFrac();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      // Tilebox
      Box bx = mfi.tilebox();

      Array4<Real      > const& vel_array  = vel_in[lev]->array(mfi);
      Array4<Real const> const& txfr_array = txfr_in[lev]->const_array(mfi);
      Array4<Real const> const& ro_array   = rho_in[lev]->const_array(mfi);
      Array4<Real const> const& ep_array   = ep_g_in[lev]->const_array(mfi);

      amrex::ParallelFor(bx,[dt,vel_array,txfr_array,ro_array,ep_array]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real orop  = dt / (ro_array(i,j,k) * ep_array(i,j,k));
        Real denom = 1.0 / (1.0 + txfr_array(i,j,k,3) * orop);

        vel_array(i,j,k,0) = (vel_array(i,j,k,0) + txfr_array(i,j,k,Transfer::velx) * orop) * denom;
        vel_array(i,j,k,1) = (vel_array(i,j,k,1) + txfr_array(i,j,k,Transfer::vely) * orop) * denom;
        vel_array(i,j,k,2) = (vel_array(i,j,k,2) + txfr_array(i,j,k,Transfer::velz) * orop) * denom;
      });

      if (advect_enthalpy) {

        const int fluid_is_a_mixture = fluid.is_a_mixture;

        Array4<Real const> dummy_arr;

        Array4<Real      > const& hg_array  = h_g_in[lev]->array(mfi);
        Array4<Real      > const& Tg_array  = T_g_in[lev]->array(mfi);
        Array4<Real const> const& Xgk_array = fluid_is_a_mixture ? X_gk_in[lev]->const_array(mfi) : dummy_arr;

        auto const& flags_arr = flags.const_array(mfi);
        auto const& volfrac_arr = volfrac.const_array(mfi);

        const int nspecies_g = fluid.nspecies;

        amrex::ParallelFor(bx,[dt,hg_array,Tg_array,txfr_array,ro_array,ep_array,
            fluid_parms,Xgk_array,nspecies_g,fluid_is_a_mixture,flags_arr,
            volfrac_arr,run_on_device]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

          if (!cell_is_covered) {

            const Real hg = hg_array(i,j,k);

            const Real gammaTp = txfr_array(i,j,k,Transfer::gammaTp);
            const Real gamma = txfr_array(i,j,k,Transfer::gamma);

            const Real epg_loc = ep_array(i,j,k);
            const Real vfrac   = volfrac_arr(i,j,k);

            const Real ep_ro_g = epg_loc*ro_array(i,j,k);

            // ************************************************************
            // Newton-Raphson solver for solving implicit equation for
            // temperature
            // ************************************************************
            // Residual computation
            auto R = [&] AMREX_GPU_DEVICE (Real Tg_arg)
            {
              Real hg_loc(0);

              if (!fluid_is_a_mixture) {

                hg_loc = run_on_device ?
                  fluid_parms.calc_h_g<RunOn::Device>(Tg_arg, cell_is_covered) :
                  fluid_parms.calc_h_g<RunOn::Host>(Tg_arg, cell_is_covered);

              } else {

                for (int n(0); n < nspecies_g; ++n) {
                  const Real h_gk = run_on_device ?
                    fluid_parms.calc_h_gk<RunOn::Device>(Tg_arg, n, cell_is_covered) :
                    fluid_parms.calc_h_gk<RunOn::Host>(Tg_arg, n, cell_is_covered);

                  hg_loc += Xgk_array(i,j,k,n)*h_gk;
                }
              }

              return ep_ro_g*(hg_loc - hg) + dt*gamma*Tg_arg - dt*gammaTp;
            };

            // Partial derivative computation
            auto partial_R = [&] AMREX_GPU_DEVICE (Real Tg_arg)
            {
              Real gradient(0);

              if (!fluid_is_a_mixture) {

                gradient = run_on_device ?
                  fluid_parms.calc_partial_h_g<RunOn::Device>(Tg_arg) :
                  fluid_parms.calc_partial_h_g<RunOn::Host>(Tg_arg);

              } else {

                for (int n(0); n < nspecies_g; ++n) {
                  const Real h_gk = run_on_device ?
                    fluid_parms.calc_partial_h_gk<RunOn::Device>(Tg_arg,n) :
                    fluid_parms.calc_partial_h_gk<RunOn::Host>(Tg_arg,n);

                  gradient += Xgk_array(i,j,k,n)*h_gk;
                }
              }

              return ep_ro_g*gradient + dt*gamma;
            };

            Real Tg_old = Tg_array(i,j,k);

            Real Tg_new(Tg_old);

            int solver_iterations(0);

            {
              DampedNewton::DampingFactor damping_factor(0., 0.);
              solver_iterations = 
                DampedNewton::solve(Tg_new, R, partial_R, damping_factor(epg_loc, vfrac),
                                    1.e-8, 1.e-8, 500);

            } if (solver_iterations == 500) {

              DampedNewton::DampingFactor damping_factor(1., 0.);
              solver_iterations =
                DampedNewton::solve(Tg_new, R, partial_R, damping_factor(epg_loc, vfrac),
                                    1.e-7, 1.e-7, 500);

            } if (solver_iterations == 500) {

              DampedNewton::DampingFactor damping_factor(1., 1.);
              solver_iterations =
                DampedNewton::solve(Tg_new, R, partial_R, damping_factor(epg_loc, vfrac),
                                    1.e-6, 1.e-6, 500);

            } if (solver_iterations == 500) {
              amrex::Abort("Damped-Newton solver did not converge");
            }


            Tg_array(i,j,k) = Tg_new;

            Real hg_new(0.);

            if (!fluid_is_a_mixture) {

              hg_new = run_on_device ?
                fluid_parms.calc_h_g<RunOn::Device>(Tg_new, cell_is_covered) :
                fluid_parms.calc_h_g<RunOn::Host>(Tg_new, cell_is_covered);

            } else {

              for (int n(0); n < nspecies_g; ++n) {
                const Real h_gk = run_on_device ?
                  fluid_parms.calc_h_gk<RunOn::Device>(Tg_new, n, cell_is_covered) :
                  fluid_parms.calc_h_gk<RunOn::Host>(Tg_new, n, cell_is_covered);

                hg_new += Xgk_array(i,j,k,n)*h_gk;
              }
            }

            hg_array(i,j,k) = hg_new;
          }
        });
      }
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

    amrex::Vector<int> condition1(finest_level+1);
    amrex::Vector<int> condition2(finest_level+1);

    Real time = 0.;

    mfix_set_velocity_bcs(time, get_vel_g(), 0);

    //
    // Make sure velocity is up to date
    //
    for (int lev = 0; lev <= finest_level; lev++)
    {

       //
       // Use temporaries to store the difference
       // between current and previous solution
       //
       MultiFab temp_vel(m_leveldata[lev]->vel_g->boxArray(),
                         m_leveldata[lev]->vel_g->DistributionMap(), 3, 0);
       MultiFab::LinComb(temp_vel, 1.0, *m_leveldata[lev]->vel_g, 0, -1.0,
                         *m_leveldata[lev]->vel_go, 0, 0, 3, 0);

       MultiFab tmp;

       const BoxArray & nd_grid = amrex::convert(grids[lev], IntVect{1,1,1});
       tmp.define(nd_grid, dmap[lev], 1, 0);

       MultiFab::LinComb(tmp, 1.0, *m_leveldata[lev]->p_g, 0, -1.0,
                         *m_leveldata[lev]->p_go, 0, 0, 1, 0);

       Real delta_u = temp_vel.norm0(0,0,false,true);
       Real delta_v = temp_vel.norm0(1,0,false,true);
       Real delta_w = temp_vel.norm0(2,0,false,true);
       Real delta_p = temp_vel.norm0(0,0,false,true);

       Real tol = steady_state_tol;

       condition1[lev] = (delta_u < tol*dt) && (delta_v < tol*dt ) && (delta_w < tol*dt);

       //
       // Second stop condition
       //
       Periodicity period = geom[lev].periodicity();

       Real du_n1 = temp_vel.norm1(0,period,true);
       Real dv_n1 = temp_vel.norm1(1,period,true);
       Real dw_n1 = temp_vel.norm1(2,period,true);
       Real dp_n1 =      tmp.norm1(0,period,true);
       Real uo_n1 = m_leveldata[lev]->vel_go->norm1(0,period,true);
       Real vo_n1 = m_leveldata[lev]->vel_go->norm1(1,period,true);
       Real wo_n1 = m_leveldata[lev]->vel_go->norm1(2,period,true);
       Real po_n1 = m_leveldata[lev]->p_go->norm1(0,period,true);

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
    for (int lev = 0; lev <= finest_level; lev++)
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
