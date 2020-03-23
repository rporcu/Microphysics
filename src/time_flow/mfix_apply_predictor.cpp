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
                            Real l_dt,
                            bool proj_2)
{
    // We use the new-time value for things computed on the "*" state
    Real new_time = time + l_dt;

    // *************************************************************************************
    // Allocate space for half-time density
    // *************************************************************************************
    Vector<MultiFab> density_nph;
    for (int lev = 0; lev <= finest_level; ++lev)
        density_nph.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(),  *ebfactory[lev]);

    // *************************************************************************************
    // Compute the explicit advective term R_u^n
    // *************************************************************************************
    mfix_compute_convective_term(conv_u_old, conv_s_old, get_vel_g_old(),
                                 get_ep_g(), get_ro_g_old(), get_trac_old(), time);

    // *************************************************************************************
    // Modify the update so they are the updates to u and s, not to (rho u) and (rho s)
    // *************************************************************************************
    // FOR NOW WE STILL DIVIDE BY EP_G BUT WE DON'T WANT TO KEEP DOING THIS!
    for (int lev = 0; lev < nlev; lev++)
      for (int i = 0; i < 3; i++)
        MultiFab::Divide(*conv_u_old[lev], *(m_leveldata[lev]->ep_g), 0, i, 1, 0);

    // *************************************************************************************
    // Compute explicit diffusive update
    // *************************************************************************************
    int explicit_diffusion_pred = 1;

    if (explicit_diffusion_pred == 1)
    {
        //mfix_set_velocity_bcs(time, vel_go, 0);
        diffusion_op->ComputeDivTau(divtau_old, get_vel_g_old(), get_ro_g(),
                                    get_ep_g(), get_mu_g());

        diffusion_op->ComputeLapS(laps_old, get_trac_old(), get_ro_g(),
                                  get_ep_g(), mu_s);

    } else {
       for (int lev = 0; lev < nlev; lev++)
       {
          divtau_old[lev]->setVal(0.);
            laps_old[lev]->setVal(0.);
       }
    }

    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    if (!advect_density)
    {
        for (int lev = 0; lev <= finest_level; lev++)
            MultiFab::Copy(density_nph[lev], *(m_leveldata[lev]->ro_go), 0, 0, 1, 1);

    } else {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*ld.vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real  const> const& rho_o  = ld.ro_go->const_array(mfi);
                Array4<Real> const& rho_new       = ld.ro_g->array(mfi);
                Array4<Real> const& rho_nph       = density_nph[lev].array(mfi);
                Array4<Real> const& epg           = ld.ep_g->array(mfi);
                Array4<Real const> const& drdt    = conv_s_old[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rho_new(i,j,k) = epg(i,j,k)*rho_o(i,j,k) + l_dt * drdt(i,j,k);
                    rho_new(i,j,k) = rho_new(i,j,k) / epg(i,j,k);

                    rho_nph(i,j,k) = 0.5 * (rho_o(i,j,k) + rho_new(i,j,k));
                });
            } // mfi
        } // lev

    } // not constant density

    // *************************************************************************************
    // Update tracer(s)
    // *************************************************************************************
    int l_ntrac = ntrac;
    if (advect_tracer)
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*ld.vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real const> const& tra_o  = ld.trac_o->const_array(mfi);
                Array4<Real      > const& tra_n  = ld.trac->array(mfi);
                Array4<Real const> const& rho_o  = ld.ro_go->const_array(mfi);
                Array4<Real const> const& rho_n  = ld.ro_g->const_array(mfi);
                Array4<Real> const& epg          = ld.ep_g->array(mfi);
                Array4<Real const> const& dtdt   = conv_s_old[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int n = 0; n < l_ntrac; ++n)
                    {
                        tra_n(i,j,k) = (rho_o(i,j,k)*epg(i,j,k))*tra_o(i,j,k,n) + l_dt * dtdt(i,j,k);
                        tra_n(i,j,k) = tra_n(i,j,k) / (rho_n(i,j,k)*epg(i,j,k));
                    }
                });
            } // mfi
        } // lev
    } // advect_tracer

    // *************************************************************************************
    // Update velocity
    // *************************************************************************************
    for (int lev = 0; lev < nlev; lev++)
    {
        EB_set_covered(*divtau_old[lev], 0, divtau_old[lev]->nComp(), divtau_old[lev]->nGrow(), 0.0);

        // First add the convective term
        MultiFab::Saxpy(*m_leveldata[lev]->vel_g, l_dt, *conv_u_old[lev], 0, 0, 3, 0);

        // Add the explicit diffusion terms
        if (explicit_diffusion_pred == 1)
           MultiFab::Saxpy(*m_leveldata[lev]->vel_g, l_dt, *divtau_old[lev], 0, 0, 3, 0);
    }

    // *************************************************************************************
    // Add source terms to velocity
    // *************************************************************************************
    mfix_add_gravity_and_gp(l_dt);

    // *************************************************************************************
    // Add the drag term implicitly
    // *************************************************************************************
    if (DEM::solve)
        mfix_add_drag_implicit(l_dt);

    // *************************************************************************************
    // If doing implicit diffusion, solve here for u^*
    // Note we multiply ep_g by ro_g so that we pass in a single array holding (ro_g * ep_g)
    // *************************************************************************************
    if (explicit_diffusion_pred == 0)
    {
      mfix_set_density_bcs(time, get_ro_g());
      mfix_set_scalar_bcs(time, get_trac(), get_mu_g());

      for (int lev = 0; lev < nlev; lev++)
        MultiFab::Multiply(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                           0, 0, 1, m_leveldata[lev]->ep_g->nGrow());

      mfix_set_velocity_bcs(new_time, get_vel_g(), 0);
      diffusion_op->diffuse_velocity(get_vel_g(), get_ep_g(), get_mu_g(), l_dt);

      // mfix_set_tracer_bcs (new_time, trac, 0);
      diffusion_op->diffuse_scalar(get_trac(), get_ep_g(), mu_s, l_dt);

      for (int lev = 0; lev < nlev; lev++)
          MultiFab::Divide(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                           0, 0, 1, m_leveldata[lev]->ep_g->nGrow());
    }

    // *************************************************************************************
    // Project velocity field -- depdt=0 for now
    // *************************************************************************************
    Vector< MultiFab* > depdt(nlev);
    for (int lev(0); lev < nlev; ++lev)
      depdt[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0, 1).release();

    mfix_apply_nodal_projection(depdt, new_time, l_dt, proj_2);

    // *************************************************************************************
    // Correct small cells
    // *************************************************************************************
    mfix_correct_small_cells (get_vel_g());

    for (int lev(0); lev < nlev; ++lev)
      delete depdt[lev];
}
