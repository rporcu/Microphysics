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
//                    + dt * ( g - grad(p_g+p0)/ro_g )
//
//  2. Add drag term implicitly ( AKA implicit part of particle smomentum exchange )
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
    // Compute the explicit advective terms
    // Note that "conv_u_old" returns update to (ep_g u)
    // Note that "conv_s_old" returns update to (ep_g rho) and (ep_g rho tracer)
    // *************************************************************************************
    mfix_compute_convective_term(conv_u_old, conv_s_old, get_vel_g_old(),
                                 get_ep_g(), get_ro_g_old(), get_trac_old(), time);

    // *************************************************************************************
    // Compute explicit diffusive update
    // *************************************************************************************
    bool explicit_diffusion_pred = true;

    if (explicit_diffusion_pred)
    {
        //mfix_set_velocity_bcs(time, vel_go, 0);
        diffusion_op->ComputeDivTau(divtau_old, get_vel_g_old(), get_ro_g(),
                                    get_ep_g(), get_mu_g());

        diffusion_op->ComputeLapS(laps_old, get_trac_old(), get_ro_g(),
                                  get_ep_g(), mu_s);

        for (int lev = 0; lev <= finest_level; lev++)
        {
            EB_set_covered(*divtau_old[lev], 0, divtau_old[lev]->nComp(), divtau_old[lev]->nGrow(), 0.0);
            EB_set_covered(  *laps_old[lev], 0,   laps_old[lev]->nComp(),   laps_old[lev]->nGrow(), 0.0);
        }

    } else {
       for (int lev = 0; lev <= finest_level; lev++)
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
                Array4<Real      > const& rho_new = ld.ro_g->array(mfi);
                Array4<Real      > const& rho_nph = density_nph[lev].array(mfi);
                Array4<Real const> const& rho_o   = ld.ro_go->const_array(mfi);
                Array4<Real const> const& epg     = ld.ep_g->const_array(mfi);
                Array4<Real const> const& drdt_o  = conv_s_old[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    int conv_comp = 0;
                    rho_new(i,j,k) = epg(i,j,k)*rho_o(i,j,k) + l_dt * drdt_o(i,j,k,conv_comp);
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
                Array4<Real      > const& tra_n  = ld.trac->array(mfi);
                Array4<Real const> const& tra_o  = ld.trac_o->const_array(mfi);
                Array4<Real const> const& rho_o  = ld.ro_go->const_array(mfi);
                Array4<Real const> const& rho_n  = ld.ro_g->const_array(mfi);
                Array4<Real const> const& laps_o = laps_old[lev]->const_array(mfi);
                Array4<Real const> const& epg    = ld.ep_g->const_array(mfi);
                Array4<Real const> const& dtdt_o = conv_s_old[lev]->const_array(mfi);

                if (explicit_diffusion_pred)
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        for (int n = 0; n < l_ntrac; ++n)
                        {
                            int conv_comp = 1+n;
                            tra_n(i,j,k,n) = (rho_o(i,j,k)*epg(i,j,k))*tra_o(i,j,k,n)
                                             + l_dt * (dtdt_o(i,j,k,conv_comp) + laps_o(i,j,k));
                            tra_n(i,j,k,n) = tra_n(i,j,k,n) / (rho_n(i,j,k)*epg(i,j,k));
                        }
                    });
                } else { // Fully implicit diffusion
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        for (int n = 0; n < l_ntrac; ++n)
                        {
                            int conv_comp = 1+n;
                            tra_n(i,j,k,n) = (rho_o(i,j,k)*epg(i,j,k))*tra_o(i,j,k,n) + l_dt * dtdt_o(i,j,k,conv_comp);
                            tra_n(i,j,k,n) = tra_n(i,j,k,n) / (rho_n(i,j,k)*epg(i,j,k));
                        }
                    });
                }
            } // mfi
        } // lev
    } // advect_tracer

    // *************************************************************************************
    // Update velocity with convective update, diffusive update, gp and gravity source terms
    // *************************************************************************************
    for (int lev = 0; lev <= finest_level; lev++)
    {

       auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(*m_leveldata[lev]->vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
         // Tilebox
         Box bx = mfi.tilebox ();

         Array4<Real      > const& vel_n    = ld.vel_g->array(mfi);
         Array4<Real const> const& vel_o    = ld.vel_go->const_array(mfi);
         Array4<Real const> const& divtau_o = divtau_old[lev]->const_array(mfi);
         Array4<Real const> const& dudt_o   = conv_u_old[lev]->const_array(mfi);
         Array4<Real const> const& gp       = ld.gp->const_array(mfi);
         Array4<Real const> const& rho_nph  = density_nph[lev].const_array(mfi);
         Array4<Real const> const& epg      = ld.ep_g->const_array(mfi);

         // We need this until we remove static attribute from mfix::gravity
         const RealVect gp0_dev(gp0);
         const RealVect gravity_dev(gravity);

         if (explicit_diffusion_pred)
         {
             amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                 vel_n(i,j,k,0) = epg(i,j,k)*vel_o(i,j,k,0) + l_dt * dudt_o(i,j,k,0);
                 vel_n(i,j,k,1) = epg(i,j,k)*vel_o(i,j,k,1) + l_dt * dudt_o(i,j,k,1);
                 vel_n(i,j,k,2) = epg(i,j,k)*vel_o(i,j,k,2) + l_dt * dudt_o(i,j,k,2);
    
                 vel_n(i,j,k,0) = vel_n(i,j,k,0) / epg(i,j,k);
                 vel_n(i,j,k,1) = vel_n(i,j,k,1) / epg(i,j,k);
                 vel_n(i,j,k,2) = vel_n(i,j,k,2) / epg(i,j,k);
    
                 vel_n(i,j,k,0) += l_dt * divtau_o(i,j,k,0);
                 vel_n(i,j,k,1) += l_dt * divtau_o(i,j,k,1);
                 vel_n(i,j,k,2) += l_dt * divtau_o(i,j,k,2);

                 Real inv_dens = 1.0 / rho_nph(i,j,k);
                 vel_n(i,j,k,0) += l_dt * (gravity_dev[0]-(gp(i,j,k,0)+gp0_dev[0])*inv_dens);
                 vel_n(i,j,k,1) += l_dt * (gravity_dev[1]-(gp(i,j,k,1)+gp0_dev[1])*inv_dens);
                 vel_n(i,j,k,2) += l_dt * (gravity_dev[2]-(gp(i,j,k,2)+gp0_dev[2])*inv_dens);
             });

         } else { // Fully implicit

             amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                 vel_n(i,j,k,0) = epg(i,j,k)*vel_o(i,j,k,0) + l_dt * dudt_o(i,j,k,0);
                 vel_n(i,j,k,1) = epg(i,j,k)*vel_o(i,j,k,1) + l_dt * dudt_o(i,j,k,1);
                 vel_n(i,j,k,2) = epg(i,j,k)*vel_o(i,j,k,2) + l_dt * dudt_o(i,j,k,2);
    
                 vel_n(i,j,k,0) = vel_n(i,j,k,0) / epg(i,j,k);
                 vel_n(i,j,k,1) = vel_n(i,j,k,1) / epg(i,j,k);
                 vel_n(i,j,k,2) = vel_n(i,j,k,2) / epg(i,j,k);

                 Real inv_dens = 1.0 / rho_nph(i,j,k);
                 vel_n(i,j,k,0) += l_dt * (gravity_dev[0]-(gp(i,j,k,0)+gp0_dev[0])*inv_dens);
                 vel_n(i,j,k,1) += l_dt * (gravity_dev[1]-(gp(i,j,k,1)+gp0_dev[1])*inv_dens);
                 vel_n(i,j,k,2) += l_dt * (gravity_dev[2]-(gp(i,j,k,2)+gp0_dev[2])*inv_dens);
             });
         }
       } // mfi
    } // lev

    // *************************************************************************************
    // Add the drag term implicitly
    // *************************************************************************************
    if (DEM::solve)
        mfix_add_drag_implicit(l_dt);

    // *************************************************************************************
    // If doing implicit diffusion, solve here for u^*
    // Note we multiply ep_g by ro_g so that we pass in a single array holding (ro_g * ep_g)
    // *************************************************************************************
    if (!explicit_diffusion_pred)
    {
      mfix_set_density_bcs(time, get_ro_g());
      mfix_set_scalar_bcs(time, get_trac(), get_mu_g());

      for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Multiply(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                           0, 0, 1, m_leveldata[lev]->ep_g->nGrow());

      mfix_set_velocity_bcs(new_time, get_vel_g(), 0);
      diffusion_op->diffuse_velocity(get_vel_g(), get_ep_g(), get_mu_g(), l_dt);

      // mfix_set_tracer_bcs (new_time, trac, 0);
      diffusion_op->diffuse_scalar(get_trac(), get_ep_g(), mu_s, l_dt);

      for (int lev = 0; lev <= finest_level; lev++)
          MultiFab::Divide(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                           0, 0, 1, m_leveldata[lev]->ep_g->nGrow());
    }

    // *************************************************************************************
    // Project velocity field -- depdt=0 for now
    // *************************************************************************************
    Vector< MultiFab* > depdt(finest_level+1);
    for (int lev(0); lev <= finest_level; ++lev)
      depdt[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0, 1).release();

    mfix_apply_nodal_projection(depdt, new_time, l_dt, proj_2);

    // *************************************************************************************
    // Correct small cells
    // *************************************************************************************
    mfix_correct_small_cells (get_vel_g());

    for (int lev(0); lev <= finest_level; ++lev)
      delete depdt[lev];
}
