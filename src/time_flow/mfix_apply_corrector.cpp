#include <mfix_F.H>
#include <mfix.H>

#include <AMReX_VisMF.H>
#include <MFIX_MFHelpers.H>
#include <MFIX_DEM_Parms.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif
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
mfix::mfix_apply_corrector (Vector< MultiFab* >& conv_u_old,
                            Vector< MultiFab* >& conv_s_old,
                            Vector< MultiFab* >& divtau_old,
                            Vector< MultiFab* >&   laps_old,
                            Real time, Real dt, bool proj_2)
{
    BL_PROFILE("mfix::mfix_apply_corrector");

    // We use the new-time value for things computed on the "*" state
    Real new_time = time + dt;

    // Create temporary multifabs to hold the new-time conv and divtau
    Vector< MultiFab* > conv_u;
    Vector< MultiFab* > conv_s;
    Vector< MultiFab* > divtau;

    conv_u.resize(nlev);
    conv_s.resize(nlev);
    divtau.resize(nlev);

    for (int lev = 0; lev < nlev; lev++)
    {
       conv_u[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
       conv_s[lev] = new MultiFab(grids[lev], dmap[lev], 2, 0, MFInfo(), *ebfactory[lev]);
       divtau[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);

       conv_u[lev]->setVal(0.0);
       conv_s[lev]->setVal(0.0);
       divtau[lev]->setVal(0.0);
    }

    // Compute the explicit advective term R_u^*
    mfix_compute_convective_term(conv_u, conv_s, get_vel_g(), get_ep_g(),
                                 get_ro_g(), get_trac(), new_time);

    // FOR NOW WE STILL DIVIDE BY EP_G BUT WE DON'T WANT TO KEEP DOING THIS!
    for (int lev = 0; lev < nlev; lev++)
      for (int i = 0; i < 3; i++)
        MultiFab::Divide(*conv_u[lev], *(m_leveldata[lev]->ep_g), 0, i, 1, 0);

    for (int lev = 0; lev < nlev; lev++)
    {
      for (int i = 0; i < 2; i++)
        MultiFab::Divide(*conv_s[lev], *(m_leveldata[lev]->ep_g), 0, i, 1, 0);

      // Make sure to do this multiply before we update density!
      if (advect_tracer)
      {
        int conv_comp = 1;
        MultiFab::Multiply(*m_leveldata[lev]->trac_o,
                           *m_leveldata[lev]->ro_go, 0, 0, 1, 0);

        // Add the convective terms so trac = trac_o + dt/2 (R_s^* + R_s^n)
        MultiFab::LinComb(*m_leveldata[lev]->trac, 1.0,
                          *m_leveldata[lev]->trac_o, 0, dt/2.0, *conv_s[lev],
                          conv_comp, 0, 1, 0);
        MultiFab::Saxpy(*m_leveldata[lev]->trac, dt/2.0, *conv_s_old[lev],
                        conv_comp, 0, 1, 0);
      }

      if (advect_density)
      {
         int conv_comp = 0;
         // Add the convective terms so trac = trac_o + dt/2 (R_s^* + R_s^n)
         MultiFab::LinComb(*m_leveldata[lev]->ro_g, 1.0,
                           *m_leveldata[lev]->ro_go, 0, dt/2.0, *conv_s[lev],
                           conv_comp, 0, 1, 0);
         MultiFab::Saxpy(*m_leveldata[lev]->ro_g, dt/2.0, *conv_s_old[lev],
                         conv_comp, 0, 1, 0);
      }

      // Make sure to do this divide after we update density!
      if (advect_tracer)
         MultiFab::Divide(*m_leveldata[lev]->trac,
                          *m_leveldata[lev]->ro_g, 0, 0, 1, 0);

      // Add the convective terms so u_g = u_go + dt/2 (R_u^* + R_u^n)
      MultiFab::LinComb(*m_leveldata[lev]->vel_g, 1.0,
                        *m_leveldata[lev]->vel_go, 0, dt/2.0,
                        *conv_u[lev], 0, 0, 3, 0);
      MultiFab::Saxpy(*m_leveldata[lev]->vel_g, dt/2.0,
                      *conv_u_old[lev], 0, 0, 3, 0);
    }

    // Add the explicit diffusion term so u_g = u_go + dt/2 (R_u^* + R_u^n) + dt/2 (Lu)^n
    for (int lev = 0; lev < nlev; lev++)
        MultiFab::Saxpy(*m_leveldata[lev]->vel_g, dt/2.0, *divtau_old[lev], 0, 0, 3, 0);

    // Add source terms
    mfix_add_gravity_and_gp(dt);

    // Add the drag term implicitly
    if (DEM::solve)
        mfix_add_drag_implicit(dt);

    //
    // Solve for u^star s.t. u^star = u_go + dt/2 (R_u^* + R_u^n) + dt/2 (Lu)^n + dt/2 (Lu)^star
    // Note we multiply ep_g by ro_g so that we pass in a single array holding (ro_g * ep_g)
    //
    for (int lev = 0; lev < nlev; lev++)
        MultiFab::Multiply(*m_leveldata[lev]->ep_g,
                           *m_leveldata[lev]->ro_g, 0, 0, 1,
                            m_leveldata[lev]->ep_g->nGrow());

    diffusion_op->diffuse_velocity(get_vel_g(), get_ep_g(), get_mu_g(), 0.5*dt);

    // mfix_set_tracer_bcs (new_time, trac, 0);
    diffusion_op->diffuse_scalar(get_trac(), get_ep_g(), mu_s, dt);

    for (int lev = 0; lev < nlev; lev++)
        MultiFab::Divide(*m_leveldata[lev]->ep_g,
                         *m_leveldata[lev]->ro_g, 0, 0, 1,
                         m_leveldata[lev]->ep_g->nGrow());

    //
    // Apply projection -- depdt=0 for now
    //
    Vector< MultiFab* > depdt(nlev);
    for (int lev(0); lev < nlev; ++lev)
        depdt[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0, 1).release();

    mfix_apply_nodal_projection(depdt, new_time, dt, proj_2);
    mfix_correct_small_cells(get_vel_g());

    for (int lev(0); lev < nlev; ++lev)
      delete depdt[lev];

    //mfix_set_velocity_bcs(new_time, vel_g, 0);
    
    for (int lev = 0; lev < nlev; lev++)
    {
       delete conv_u[lev];
       delete conv_s[lev];
       delete divtau[lev];
    }
}
