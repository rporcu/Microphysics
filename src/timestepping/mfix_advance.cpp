#include <mfix.H>

#include <AMReX_VisMF.H>
#include <mfix_mf_helpers.H>
#include <mfix_dem.H>
#include <mfix_fluid.H>
#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_eb.H>
#include <mfix_pic.H>
#include <mfix_solvers.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

namespace EnthalpyImplicitUpdate {

struct Residue
{
  AMREX_GPU_HOST_DEVICE
  Residue (const int& i,
           const int& j,
           const int& k,
           const int& fluid_is_a_mixture,
           const int& nspecies_g,
           const MFIXFluidParms& fluid_parms,
           const Array4<const Real>& Xgk_array,
           const Real& hg,
           const Real& ep_ro_g,
           const Real& gamma,
           const Real& gammaTp,
           const Real& dt)
    : m_i(i)
    , m_j(j)
    , m_k(k)
    , m_fluid_is_a_mixture(fluid_is_a_mixture)
    , m_nspecies_g(nspecies_g)
    , m_fluid_parms(fluid_parms)
    , m_Xgk_array(Xgk_array)
    , m_hg(hg)
    , m_ep_ro_g(ep_ro_g)
    , m_gamma(gamma)
    , m_gammaTp(gammaTp)
    , m_dt(dt)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real operator() (const Real& Tg_arg)
  {
    Real hg_loc(0);

    if (!m_fluid_is_a_mixture) {
      hg_loc = m_fluid_parms.calc_h_g<run_on>(Tg_arg);
    } else {
      for (int n_g(0); n_g < m_nspecies_g; ++n_g) {
        const Real h_gk = m_fluid_parms.calc_h_gk<run_on>(Tg_arg, n_g);
        hg_loc += m_Xgk_array(m_i,m_j,m_k,n_g)*h_gk;
      }
    }

    return m_ep_ro_g*(hg_loc - m_hg) + m_dt*m_gamma*Tg_arg - m_dt*m_gammaTp;
  }

  const int& m_i; const int& m_j; const int& m_k;
  const int& m_fluid_is_a_mixture;
  const int& m_nspecies_g;
  const MFIXFluidParms& m_fluid_parms;
  const Array4<const Real>& m_Xgk_array;
  const Real& m_hg;
  const Real& m_ep_ro_g;
  const Real& m_gamma;
  const Real& m_gammaTp;
  const Real& m_dt;
};

struct Gradient
{
  AMREX_GPU_HOST_DEVICE
  Gradient (const int& i,
            const int& j,
            const int& k,
            const int& fluid_is_a_mixture,
            const int& nspecies_g,
            const MFIXFluidParms& fluid_parms,
            const Array4<const Real>& Xgk_array,
            const Real& ep_ro_g,
            const Real& gamma,
            const Real& dt)
    : m_i(i)
    , m_j(j)
    , m_k(k)
    , m_fluid_is_a_mixture(fluid_is_a_mixture)
    , m_nspecies_g(nspecies_g)
    , m_fluid_parms(fluid_parms)
    , m_Xgk_array(Xgk_array)
    , m_ep_ro_g(ep_ro_g)
    , m_gamma(gamma)
    , m_dt(dt)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real operator() (const Real& Tg_arg)
  {
    Real gradient(0);

    if (!m_fluid_is_a_mixture) {
      gradient = m_fluid_parms.calc_partial_h_g<run_on>(Tg_arg);
    } else {
      for (int n_g(0); n_g < m_nspecies_g; ++n_g) {
        const Real h_gk = m_fluid_parms.calc_partial_h_gk<run_on>(Tg_arg,n_g);
        gradient += m_Xgk_array(m_i,m_j,m_k,n_g)*h_gk;
      }
    }

    return m_ep_ro_g*gradient + m_dt*m_gamma;
  }

  const int& m_i; const int& m_j; const int& m_k;
  const int& m_fluid_is_a_mixture;
  const int& m_nspecies_g;
  const MFIXFluidParms& m_fluid_parms;
  const Array4<const Real>& m_Xgk_array;
  const Real& m_ep_ro_g;
  const Real& m_gamma;
  const Real& m_dt;
};

} // end namespace EnthalpyImplicitUpdate

using namespace EnthalpyImplicitUpdate;
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
    Vector< MultiFab* > eb_flow_vel(finest_level+1, nullptr);
    Vector< MultiFab* > depdt(finest_level+1);

    for (int lev(0); lev <= finest_level; ++lev) {
      depdt[lev] = new MultiFab(grids[lev], dmap[lev], 1, 1, MFInfo(), EBFactory(lev));
      depdt[lev]->setVal(0.);

      if (m_embedded_boundaries.has_flow()) {
        eb_flow_vel[lev] = new MultiFab(grids[lev], dmap[lev], 3, nghost_state(), MFInfo(), *ebfactory[lev]);
        eb_flow_vel[lev]->setVal(0.0);
      }
    }

    if (m_embedded_boundaries.has_flow()) {

       m_boundary_conditions.set_eb_velocity_bcs(time, eb_flow_vel);
    }

    mfix_apply_nodal_projection(depdt, time, dummy_dt, dummy_dt, proj_2,
                                get_vel_g_old(), get_vel_g(), get_p_g(), get_gp(),
                                get_ep_g(), get_txfr(), get_ro_g_const(),
                                GetVecOfConstPtrs(eb_flow_vel));

    for (int lev(0); lev <= finest_level; ++lev) {
      delete depdt[lev];
      if (m_embedded_boundaries.has_flow()) {
        delete eb_flow_vel[lev];
      }
    }

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
  m_timer.dt() = dt;

  amrex::Print() << "Doing initial pressure iterations with dt = " << dt << "\n";

  // Fill ghost nodes and reimpose boundary conditions
  m_boundary_conditions.set_velocity_bcs(time, get_vel_g(), 0);
  m_boundary_conditions.set_density_bcs(time, get_ro_g());
  m_boundary_conditions.set_tracer_bcs(time, fluid, get_trac());

  if (fluid.solve_enthalpy()) {
    m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g());
    m_boundary_conditions.set_enthalpy_bcs(time, fluid,get_h_g());
  }

  if (fluid.solve_enthalpy() && m_embedded_boundaries.fix_temperature())
    m_boundary_conditions.set_eb_temperature_bcs(get_T_g_on_eb());

  if (fluid.solve_species())
    m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk());

  // Copy vel_g and p_g into vel_go and p_go
  for (int lev = 0; lev <= finest_level; lev++) {
    MultiFab::Copy(*m_leveldata[lev]->vel_go, *m_leveldata[lev]->vel_g, 0, 0,
                   m_leveldata[lev]->vel_g->nComp(), m_leveldata[lev]->vel_g->nGrow());

    MultiFab::Copy(*m_leveldata[lev]->p_go, *m_leveldata[lev]->p_g, 0, 0,
                   m_leveldata[lev]->p_g->nComp(), m_leveldata[lev]->p_g->nGrow());
  }

  if (m_dem.solve() || m_pic.solve()) {
    mfix_calc_txfr_fluid(get_txfr(), get_ep_g(), get_ro_g(),
                         get_vel_g(), get_T_g(), get_X_gk(), get_thermodynamic_p_g(),
                         time);
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

  Vector< MultiFab* > eb_flow_vel(finest_level+1, nullptr);
  Vector< MultiFab* > eb_flow_scalars(finest_level+1, nullptr);
  Vector< MultiFab* > eb_flow_species(finest_level+1, nullptr);

  for (int lev = 0; lev <= finest_level; lev++)
  {
    conv_u[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
    conv_u[lev]->setVal(0.0);

    // 2+ntrac: one for density, ntrac for tracers and one for enthalpy
    conv_s[lev] = new MultiFab(grids[lev], dmap[lev], 2+ntrac, 0, MFInfo(), *ebfactory[lev]);
    conv_s[lev]->setVal(0.0);

    if (fluid.solve_density()) {
      ro_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
      ro_RHS_old[lev]->setVal(0.0);
    }

    if (fluid.solve_tracer()) {
      lap_trac_old[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), *ebfactory[lev]);
      lap_trac_old[lev]->setVal(0.0);
    }

    if (fluid.solve_enthalpy()) {
      enthalpy_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
      enthalpy_RHS_old[lev]->setVal(0.0);
      lap_T_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
      lap_T_old[lev]->setVal(0.0);
    }

    if (fluid.solve_species()) {
      div_J_old[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies(), 0, MFInfo(), *ebfactory[lev]);
      div_J_old[lev]->setVal(0.0);
    }

    if (fluid.solve_enthalpy() && fluid.solve_species()) {
      div_hJ_old[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
      div_hJ_old[lev]->setVal(0.0);
    }

    if (fluid.solve_species()) {
      conv_X[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies(), 0, MFInfo(), *ebfactory[lev]);
      conv_X[lev]->setVal(0.0);

      species_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies(), 0, MFInfo(), *ebfactory[lev]);
      species_RHS_old[lev]->setVal(0.0);
    }

    if (reactions.solve()) {
      vel_RHS_old[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
    }

    if (m_embedded_boundaries.has_flow()) {

      eb_flow_vel[lev] = new MultiFab(grids[lev], dmap[lev], 3, nghost_state(), MFInfo(), *ebfactory[lev]);
      eb_flow_vel[lev]->setVal(0.0);

      eb_flow_scalars[lev] = new MultiFab(grids[lev], dmap[lev], 2+ntrac, nghost_state(), MFInfo(), *ebfactory[lev]);
      eb_flow_scalars[lev]->setVal(0.0);

      if (fluid.solve_species()) {
        eb_flow_species[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies(), nghost_state(), MFInfo(), *ebfactory[lev]);
        eb_flow_species[lev]->setVal(0.0);
      }
    }
  }//nlev


  if (m_embedded_boundaries.has_flow()) {

     m_boundary_conditions.set_eb_velocity_bcs(time, eb_flow_vel);

     m_boundary_conditions.set_eb_scalar_bcs(fluid, eb_flow_scalars, eb_flow_species);
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
        div_J_old, div_hJ_old, rhs_pressure_g_old, rhs_pressure_g,
        eb_flow_vel, eb_flow_scalars, eb_flow_species, time, dt,
        dt_copy, proj_2, coupling_timing);

    // Reset any quantities which might have been updated
    for (int lev = 0; lev <= finest_level; lev++) {
      MultiFab::Copy(*m_leveldata[lev]->vel_g, *m_leveldata[lev]->vel_go, 0, 0,
                     m_leveldata[lev]->vel_g->nComp(), m_leveldata[lev]->vel_g->nGrow());

      if (fluid.solve_density())
        MultiFab::Copy(*m_leveldata[lev]->ro_g, *m_leveldata[lev]->ro_go, 0, 0,
                       m_leveldata[lev]->ro_g->nComp(), m_leveldata[lev]->ro_g->nGrow());

      if (fluid.solve_enthalpy()) {
        MultiFab::Copy(*m_leveldata[lev]->T_g, *m_leveldata[lev]->T_go, 0, 0,
                       m_leveldata[lev]->T_g->nComp(), m_leveldata[lev]->T_g->nGrow());
        MultiFab::Copy(*m_leveldata[lev]->h_g, *m_leveldata[lev]->h_go, 0, 0,
                       m_leveldata[lev]->h_g->nComp(), m_leveldata[lev]->h_g->nGrow());
      }

      if (fluid.solve_tracer())
        MultiFab::Copy(*m_leveldata[lev]->trac, *m_leveldata[lev]->trac_o, 0, 0,
                       m_leveldata[lev]->trac->nComp(), m_leveldata[lev]->trac->nGrow());

      if (fluid.solve_species()) {
        MultiFab::Copy(*m_leveldata[lev]->X_gk, *m_leveldata[lev]->X_gko, 0, 0,
                       m_leveldata[lev]->X_gk->nComp(), m_leveldata[lev]->X_gk->nGrow());
      }

      if (fluid.constraint_type()== MFIXFluidPhase::ConstraintType::IdealGasClosedSystem &&
          fluid.solve_enthalpy()) {
        *m_leveldata[lev]->thermodynamic_p_g = *m_leveldata[lev]->thermodynamic_p_go;
      }
    }

    // Reset the boundary values (necessary if they are time-dependent)
    m_boundary_conditions.set_velocity_bcs(time, get_vel_g(), 0);
    m_boundary_conditions.set_density_bcs(time, get_ro_g());
    m_boundary_conditions.set_tracer_bcs(time, fluid, get_trac());

    if (fluid.solve_enthalpy()) {
      m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g());
      m_boundary_conditions.set_enthalpy_bcs(time, fluid,get_h_g());
    }

    if (fluid.solve_species())
      m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk());
  }

  for (int lev = 0; lev <= finest_level; lev++)
  {
     delete conv_u[lev];
     delete conv_s[lev];

     if (fluid.solve_density())
       delete ro_RHS_old[lev];

     if (fluid.solve_tracer())
       delete lap_trac_old[lev];

     if (fluid.solve_enthalpy()) {
       delete enthalpy_RHS_old[lev];
       delete lap_T_old[lev];
     }

     if (fluid.solve_enthalpy() && fluid.solve_species()) {
       delete div_hJ_old[lev];
     }

     if (fluid.solve_species()) {
       delete conv_X[lev];
       delete species_RHS_old[lev];
       delete div_J_old[lev];
     }

     if (reactions.solve())
       delete vel_RHS_old[lev];

     if (m_embedded_boundaries.has_flow()) {
       delete eb_flow_vel[lev];
       delete eb_flow_scalars[lev];
       if (fluid.solve_species()) {
         delete eb_flow_species[lev];
       }
    }
  }
}

//
// Explicit solve for the intermediate velocity.
// Currently this means accounting for the explicit part of the fluid/particle
// momentum exchange
//
void
mfix::mfix_add_vel_txfr_explicit (Real dt,
                                  Vector<MultiFab*      > const& vel_in,
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

  for (int lev = 0; lev <= finest_level; lev++) {

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

      InterphaseTxfrIndexes txfr_idxs(fluid.nspecies(), reactions.nreactions());

      const int idx_velx_txfr = txfr_idxs.vel+0;
      const int idx_vely_txfr = txfr_idxs.vel+1;
      const int idx_velz_txfr = txfr_idxs.vel+2;
      const int idx_drag_txfr = txfr_idxs.drag_coeff;

      amrex::ParallelFor(bx,[dt,vel_array,txfr_array,ro_array,ep_array,
         idx_velx_txfr, idx_vely_txfr, idx_velz_txfr, idx_drag_txfr ]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const Real orop  = dt / (ro_array(i,j,k) * ep_array(i,j,k));

        const Real beta = txfr_array(i,j,k,idx_drag_txfr);

        const Real vel_x = vel_array(i,j,k,0);
        const Real vel_y = vel_array(i,j,k,1);
        const Real vel_z = vel_array(i,j,k,2);

        const Real drag_0 = (txfr_array(i,j,k,idx_velx_txfr) - beta*vel_x) * orop;
        const Real drag_1 = (txfr_array(i,j,k,idx_vely_txfr) - beta*vel_y) * orop;
        const Real drag_2 = (txfr_array(i,j,k,idx_velz_txfr) - beta*vel_z) * orop;

        vel_array(i,j,k,0) = vel_x + drag_0;
        vel_array(i,j,k,1) = vel_y + drag_1;
        vel_array(i,j,k,2) = vel_z + drag_2;
      });
    }
  }
}

void
mfix::mfix_add_enthalpy_txfr_explicit (Real dt,
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

  BL_PROFILE("mfix::mfix_add_enthalpy_txfr_explicit");

  const auto& fluid_parms = fluid.parameters();

  for (int lev = 0; lev <= finest_level; lev++) {

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ep_g_in[lev]->Factory());
    const auto& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*h_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      // Tilebox
      Box bx = mfi.tilebox();

      Array4<Real const> const& txfr_array = txfr_in[lev]->array(mfi);
      Array4<Real const> const&   ro_array = rho_in[lev]->array(mfi);
      Array4<Real const> const&   ep_array = ep_g_in[lev]->array(mfi);

      Array4<Real      > const& hg_array  = h_g_in[lev]->array(mfi);
      Array4<Real      > const& Tg_array  = T_g_in[lev]->array(mfi);
      Array4<Real const> const& Xgk_array = X_gk_in[lev]->const_array(mfi);

      auto const& flags_arr = flags.const_array(mfi);

      const int nspecies_g = fluid.nspecies();
      const int fluid_is_a_mixture = fluid.isMixture();

      const int is_IOProc = int(ParallelDescriptor::IOProcessor());

      InterphaseTxfrIndexes txfr_idxs(fluid.nspecies(), reactions.nreactions());

      const int idx_gammaTp_txfr = txfr_idxs.gammaTp;
      const int idx_convection_coeff_txfr = txfr_idxs.convection_coeff;

      amrex::ParallelFor(bx,[dt,hg_array,Tg_array,txfr_array,ro_array,ep_array,
          fluid_parms,Xgk_array,nspecies_g,fluid_is_a_mixture,flags_arr,
          idx_gammaTp_txfr, idx_convection_coeff_txfr,
          is_IOProc,abstol=newton_abstol,reltol=newton_reltol,maxiter=newton_maxiter]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

        if (!cell_is_covered) {
          const Real epg_loc = ep_array(i,j,k);

          const Real orop  = dt / (ro_array(i,j,k) * epg_loc);

          const Real Tg_old = Tg_array(i,j,k);
          const Real Ts     = txfr_array(i,j,k,idx_gammaTp_txfr);
          const Real gamma  = txfr_array(i,j,k,idx_convection_coeff_txfr);

          const Real hg = hg_array(i,j,k) + (Ts - gamma * Tg_old) * orop;
          hg_array(i,j,k) = hg;

          // ************************************************************
          // Newton-Raphson solver for solving implicit equation for
          // temperature
          // ************************************************************
          Newton::FluidEnthalpy::Residue residue(i, j, k, fluid_is_a_mixture,
              cell_is_covered, nspecies_g, fluid_parms, Xgk_array, hg);

          Newton::FluidEnthalpy::Gradient gradient(i, j, k, fluid_is_a_mixture,
              nspecies_g, fluid_parms, Xgk_array);

          Real Tg_new(Tg_old);

          Newton::solve(Tg_new, residue, gradient, abstol, reltol, maxiter, is_IOProc);

          Tg_array(i,j,k) = Tg_new;

          Real hg_new(0.);

          if (!fluid_is_a_mixture) {

            hg_new = fluid_parms.calc_h_g<run_on>(Tg_new);

          } else {

            for (int n(0); n < nspecies_g; ++n) {
              const Real h_gk = fluid_parms.calc_h_gk<run_on>(Tg_new, n);

              hg_new += Xgk_array(i,j,k,n)*h_gk;
            }
          }

          hg_array(i,j,k) = hg_new;
        }
      });
    }
  }
}

//
// Implicit solve for the intermediate velocity.
// Currently this means accounting for the implicit part of the fluid/particle
// momentum exchange
//
void
mfix::mfix_add_vel_txfr_implicit (Real dt,
                                  Vector<MultiFab*      > const& vel_in,
                                  Vector<MultiFab const*> const& txfr_in,
                                  Vector<MultiFab const*> const& rho_in,
                                  Vector<MultiFab const*> const& ep_g_in)
{
  /*
     This adds both components of the drag term
     So the drag term we add is beta * (particle_velocity - fluid_velocity)
                              = drag(0:2) - drag(3) * fluid_velocity
  */

  BL_PROFILE("mfix::mfix_add_vel_txfr_implicit");

  InterphaseTxfrIndexes txfr_idxs(fluid.nspecies(), reactions.nreactions());

  const int idx_velx_txfr = txfr_idxs.vel+0;
  const int idx_vely_txfr = txfr_idxs.vel+1;
  const int idx_velz_txfr = txfr_idxs.vel+2;

  for (int lev = 0; lev <= finest_level; lev++) {

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


      amrex::ParallelFor(bx,[dt,vel_array,txfr_array,ro_array,ep_array,
       idx_velx_txfr, idx_vely_txfr, idx_velz_txfr ]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real orop  = dt / (ro_array(i,j,k) * ep_array(i,j,k));
        Real denom = 1.0 / (1.0 + txfr_array(i,j,k,3) * orop);

        vel_array(i,j,k,0) = (vel_array(i,j,k,0) + txfr_array(i,j,k,idx_velx_txfr) * orop) * denom;
        vel_array(i,j,k,1) = (vel_array(i,j,k,1) + txfr_array(i,j,k,idx_vely_txfr) * orop) * denom;
        vel_array(i,j,k,2) = (vel_array(i,j,k,2) + txfr_array(i,j,k,idx_velz_txfr) * orop) * denom;
      });
    }
  }
}

void
mfix::mfix_add_enthalpy_txfr_implicit (Real dt,
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

  BL_PROFILE("mfix::mfix_add_energy_txfr_implicit");

  const auto& fluid_parms = fluid.parameters();

  InterphaseTxfrIndexes txfr_idxs(fluid.nspecies(), reactions.nreactions());

  const int idx_gammaTp_txfr = txfr_idxs.gammaTp;
  const int idx_convection_coeff_txfr = txfr_idxs.convection_coeff;

  for (int lev = 0; lev <= finest_level; lev++) {
    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ep_g_in[lev]->Factory());
    const auto& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*h_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      // Tilebox
      Box bx = mfi.tilebox();

      Array4<Real const> const& txfr_array = txfr_in[lev]->const_array(mfi);
      Array4<Real const> const& ro_array   = rho_in[lev]->const_array(mfi);
      Array4<Real const> const& ep_array   = ep_g_in[lev]->const_array(mfi);

      const int fluid_is_a_mixture = fluid.isMixture();

      Array4<Real const> dummy_arr;

      Array4<Real      > const& hg_array  = h_g_in[lev]->array(mfi);
      Array4<Real      > const& Tg_array  = T_g_in[lev]->array(mfi);
      Array4<Real const> const& Xgk_array = fluid_is_a_mixture ? X_gk_in[lev]->const_array(mfi) : dummy_arr;

      auto const& flags_arr = flags.const_array(mfi);

      const int nspecies_g = fluid.nspecies();

      const int is_IOProc = int(ParallelDescriptor::IOProcessor());

      amrex::ParallelFor(bx,[dt,hg_array,Tg_array,txfr_array,ro_array,ep_array,
          fluid_parms,Xgk_array,nspecies_g,fluid_is_a_mixture,flags_arr,
          idx_gammaTp_txfr, idx_convection_coeff_txfr,
          is_IOProc,abstol=newton_abstol,reltol=newton_reltol,maxiter=newton_maxiter]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

        if (!cell_is_covered) {

          const Real hg = hg_array(i,j,k);

          const Real gammaTp = txfr_array(i,j,k,idx_gammaTp_txfr);
          const Real gamma = txfr_array(i,j,k,idx_convection_coeff_txfr);

          const Real epg_loc = ep_array(i,j,k);

          const Real ep_ro_g = epg_loc*ro_array(i,j,k);

          // ************************************************************
          // Newton-Raphson solver for solving implicit equation for
          // temperature
          // ************************************************************
          Residue residue(i, j, k, fluid_is_a_mixture, nspecies_g, fluid_parms,
              Xgk_array, hg, ep_ro_g, gamma, gammaTp, dt);
          
          Gradient gradient(i, j, k, fluid_is_a_mixture, nspecies_g, fluid_parms,
              Xgk_array, ep_ro_g, gamma, dt);

          Real Tg_old = Tg_array(i,j,k);

          Real Tg_new(Tg_old);

          Newton::solve(Tg_new, residue, gradient, abstol, reltol, maxiter, is_IOProc);

          Tg_array(i,j,k) = Tg_new;

          Real hg_new(0.);

          if (!fluid_is_a_mixture) {

            hg_new = fluid_parms.calc_h_g<run_on>(Tg_new);

          } else {

            for (int n(0); n < nspecies_g; ++n) {
              const Real h_gk = fluid_parms.calc_h_gk<run_on>(Tg_new, n);

              hg_new += Xgk_array(i,j,k,n)*h_gk;
            }
          }

          hg_array(i,j,k) = hg_new;
        }
      });
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

    m_boundary_conditions.set_velocity_bcs(time, get_vel_g(), 0);

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
       amrex::Print() << "\nSteady state check at level " << lev << ":"
         << "\n||u-uo||/||uo|| = " << std::setw(12) << std::scientific << std::setprecision(4) << tmp1
         << "      du/dt = " << std::setw(12) << std::scientific << std::setprecision(4) << delta_u/dt
         << "\n||v-vo||/||vo|| = " << std::setw(12) << std::scientific << std::setprecision(4) << tmp2
         << "      dv/dt = " << std::setw(12) << std::scientific << std::setprecision(4) << delta_v/dt
         << "\n||w-wo||/||wo|| = " << std::setw(12) << std::scientific << std::setprecision(4) << tmp3
         << "      dw/dt = " << std::setw(12) << std::scientific << std::setprecision(4) << delta_w/dt
         << "\n||p-po||/||po|| = " << std::setw(12) << std::scientific << std::setprecision(4) << tmp4
         << "      dp/dt = " << std::setw(12) << std::scientific << std::setprecision(4) << delta_p/dt
         << "\n\n";
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
