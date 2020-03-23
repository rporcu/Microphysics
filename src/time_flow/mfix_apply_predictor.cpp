#include <mfix_F.H>
#include <mfix.H>

#include <AMReX_VisMF.H>
#include <MFIX_MFHelpers.H>
#include <MFIX_DEM_Parms.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

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
mfix::mfix_apply_predictor (Vector< MultiFab* >& conv_u_old,
                            Vector< MultiFab* >& conv_s_old,
                            Vector< MultiFab* >& divtau_old,
                            Vector< MultiFab* >&   laps_old,
                            Real time,
                            Real dt,
                            bool proj_2)
{
    // We use the new-time value for things computed on the "*" state
    Real new_time = time + dt;

    // Compute the explicit advective term R_u^n
    mfix_compute_convective_term(conv_u_old, conv_s_old, get_vel_g_old(),
                                 get_ep_g(), get_ro_g_old(), get_trac_old(), time);

    // FOR NOW WE STILL DIVIDE BY EP_G BUT WE DON'T WANT TO KEEP DOING THIS!
    for (int lev = 0; lev < nlev; lev++)
      for (int i = 0; i < 3; i++)
        MultiFab::Divide(*conv_u_old[lev], *(m_leveldata[lev]->ep_g), 0, i, 1, 0);

    for (int lev = 0; lev < nlev; lev++)
      for (int i = 0; i < 2; i++)
        MultiFab::Divide(*conv_s_old[lev], *(m_leveldata[lev]->ep_g), 0, i, 1, 0);

    int explicit_diffusion_pred = 1;

    if (explicit_diffusion_pred == 1)
    {
        //mfix_set_velocity_bcs(time, vel_go, 0);
        diffusion_op->ComputeDivTau(divtau_old, get_vel_g_old(), get_ro_g(),
                                    get_ep_g(), get_mu_g());

        // mfix_set_tracer_bcs (time, trac_o);
        diffusion_op->ComputeLapS(laps_old, get_trac_old(), get_ro_g(),
                                  get_ep_g(), mu_s);

    } else {
       for (int lev = 0; lev < nlev; lev++)
       {
          divtau_old[lev]->setVal(0.);
            laps_old[lev]->setVal(0.);
       }
    }

    for (int lev = 0; lev < nlev; lev++)
    {
        EB_set_covered(*divtau_old[lev], 0, divtau_old[lev]->nComp(), divtau_old[lev]->nGrow(), 0.0);

        // First add the convective term
        MultiFab::Saxpy(*m_leveldata[lev]->vel_g, dt, *conv_u_old[lev], 0, 0, 3, 0);

        // Make sure to do this multiply before we update density!
        if (advect_tracer)
        {
           int conv_comp = 1;
           MultiFab::Multiply(*m_leveldata[lev]->trac,
                              *m_leveldata[lev]->ro_go, 0, 0, 1, 0);
           MultiFab::Saxpy(*m_leveldata[lev]->trac, dt,
                           *conv_s_old[lev], conv_comp, 0, 1, 0);
        }

        if (advect_density)
        {
           int conv_comp = 0;
           MultiFab::Saxpy(*m_leveldata[lev]->ro_g, dt, *conv_s_old[lev],
                           conv_comp, 0, 1, 0);
        }

        // Make sure to do this divide after we update density!
        if (advect_tracer)
           MultiFab::Divide(*m_leveldata[lev]->trac, *m_leveldata[lev]->ro_g, 0, 0, 1, 0);

        // Add the explicit diffusion terms
        if (explicit_diffusion_pred == 1)
           MultiFab::Saxpy(*m_leveldata[lev]->vel_g, dt, *divtau_old[lev], 0, 0, 3, 0);
    }

    // Add source terms
    mfix_add_gravity_and_gp(dt);

    // Add the drag term implicitly
    if (DEM::solve)
        mfix_add_drag_implicit(dt);

    // If doing implicit diffusion, solve here for u^*
    // Note we multiply ep_g by ro_g so that we pass in a single array holding (ro_g * ep_g)
    if (explicit_diffusion_pred == 0)
    {
      mfix_set_density_bcs(time, get_ro_g());
      mfix_set_scalar_bcs(time, get_trac(), get_mu_g());

      for (int lev = 0; lev < nlev; lev++)
        MultiFab::Multiply(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                           0, 0, 1, m_leveldata[lev]->ep_g->nGrow());

      mfix_set_velocity_bcs(new_time, get_vel_g(), 0);
      diffusion_op->diffuse_velocity(get_vel_g(), get_ep_g(), get_mu_g(), dt);

      // mfix_set_tracer_bcs (new_time, trac, 0);
      diffusion_op->diffuse_scalar(get_trac(), get_ep_g(), mu_s, dt);

      for (int lev = 0; lev < nlev; lev++)
          MultiFab::Divide(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                           0, 0, 1, m_leveldata[lev]->ep_g->nGrow());
    }

    // Project velocity field -- depdt=0 for now
    Vector< MultiFab* > depdt(nlev);
    for (int lev(0); lev < nlev; ++lev)
      depdt[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0, 1).release();

    mfix_apply_nodal_projection(depdt, new_time, dt, proj_2);

    mfix_correct_small_cells (get_vel_g());

    for (int lev(0); lev < nlev; ++lev)
      delete depdt[lev];

    //mfix_set_velocity_bcs(new_time, vel_g, 0);
}
