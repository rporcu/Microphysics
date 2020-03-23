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
//                    + dt * ( g - grad(p_g+p0)/ro_g )
//
//     where the starred variables are computed using "predictor-step" variables.
//
//  2. Add implicit forcing term ( AKA implicit part of particles
//     momentum exchange )
//
//     vel_g = (vel_g + (drag_coeff*velp)/(ro_g*ep_g) / ( 1 + dt * drag_coeff/(ro_g*ep_g)
//
//  3. Solve for phi
//
//     div( ep_g * grad(phi) / ro_g ) = div( ep_g * vel_g / dt + grad(p_g)/ro_g )
//
//  4. Compute
//
//     vel_g = vel_g -  dt * grad(phi) / ro_g
//
//  5. Define
//
//     p_g = phi
//
void
mfix::mfix_apply_corrector (Vector< MultiFab* >& conv_u_old,
                            Vector< MultiFab* >& conv_s_old,
                            Vector< MultiFab* >& divtau_old,
                            Vector< MultiFab* >&   laps_old,
                            Real time, Real l_dt, bool proj_2)
{
    BL_PROFILE("mfix::mfix_apply_corrector");

    // We use the new-time value for things computed on the "*" state
    Real new_time = time + l_dt;

    // *************************************************************************************
    // Allocate space for half-time density
    // *************************************************************************************
    Vector<MultiFab> density_nph;
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        density_nph.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(),  *ebfactory[lev]);
    }

    // *************************************************************************************
    // Create temporary multifabs to hold the new-time conv and divtau
    // *************************************************************************************
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

    // *************************************************************************************
    // Compute the explicit advective term R_u^*
    // *************************************************************************************
    mfix_compute_convective_term(conv_u, conv_s, get_vel_g(), get_ep_g(),
                                 get_ro_g(), get_trac(), new_time);

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
                Array4<Real const> const& drdt_o  = conv_s_old[lev]->const_array(mfi);
                Array4<Real const> const& drdt    = conv_s[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rho_new(i,j,k) = epg(i,j,k)*rho_o(i,j,k) 
                                    + 0.5 * l_dt * (drdt_o(i,j,k,0) + drdt(i,j,k,0));
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
                Array4<Real const> const& dtdt_o = conv_s_old[lev]->const_array(mfi);
                Array4<Real const> const& dtdt   = conv_s[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int n = 0; n < l_ntrac; ++n)
                    {
                        int conv_comp = 1+n;
                        tra_n(i,j,k,n) = (rho_o(i,j,k)*epg(i,j,k))*tra_o(i,j,k,n) 
                                      + 0.5 * l_dt * (dtdt_o(i,j,k,conv_comp) + dtdt(i,j,k,conv_comp));
                        tra_n(i,j,k,n) = tra_n(i,j,k,n) / (rho_n(i,j,k)*epg(i,j,k));
                    }
                });
            } // mfi
        } // lev
    } // advect_tracer

    // *************************************************************************************
    // Update velocity with convective update, diffusive update, gp and gravity source terms
    // *************************************************************************************
    for (int lev = 0; lev < nlev; lev++)
    {
       auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(*m_leveldata[lev]->vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
         // Tilebox
         Box bx = mfi.tilebox ();

         Array4<Real const> const& vel_o   = ld.vel_go->const_array(mfi);
         Array4<Real      > const& vel_n   = ld.vel_g->array(mfi);
         Array4<Real const> const& lapu_o  = divtau_old[lev]->const_array(mfi);
         Array4<Real const> const& dudt_o  = conv_u_old[lev]->const_array(mfi);
         Array4<Real const> const& dudt    = conv_u[lev]->const_array(mfi);
         Array4<Real const> const& gp      = ld.gp->const_array(mfi);
         Array4<Real const> const& rho_nph = density_nph[lev].const_array(mfi);
         Array4<Real const> const& epg     = ld.ep_g->const_array(mfi);

         // We need this until we remove static attribute from mfix::gravity
         const RealVect gp0_dev(gp0);
         const RealVect gravity_dev(gravity);

         amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
             vel_n(i,j,k,0) = epg(i,j,k)*vel_o(i,j,k,0) + 0.5 * l_dt * (dudt_o(i,j,k,0)+dudt(i,j,k,0));
             vel_n(i,j,k,1) = epg(i,j,k)*vel_o(i,j,k,1) + 0.5 * l_dt * (dudt_o(i,j,k,1)+dudt(i,j,k,1));
             vel_n(i,j,k,2) = epg(i,j,k)*vel_o(i,j,k,2) + 0.5 * l_dt * (dudt_o(i,j,k,2)+dudt(i,j,k,2));

             vel_n(i,j,k,0) = vel_n(i,j,k,0) / epg(i,j,k);
             vel_n(i,j,k,1) = vel_n(i,j,k,1) / epg(i,j,k);
             vel_n(i,j,k,2) = vel_n(i,j,k,2) / epg(i,j,k);

             Real inv_dens = 1.0 / rho_nph(i,j,k);
             vel_n(i,j,k,0) += l_dt * (gravity_dev[0]-(gp(i,j,k,0)+gp0_dev[0])*inv_dens);
             vel_n(i,j,k,1) += l_dt * (gravity_dev[1]-(gp(i,j,k,1)+gp0_dev[1])*inv_dens);
             vel_n(i,j,k,2) += l_dt * (gravity_dev[2]-(gp(i,j,k,2)+gp0_dev[2])*inv_dens);
         });
       }
    }

    // *************************************************************************************
    // Add the explicit diffusion term so u_g = u_go + dt/2 (R_u^* + R_u^n) + dt/2 (Lu)^n
    // *************************************************************************************
    for (int lev = 0; lev < nlev; lev++)
        MultiFab::Saxpy(*m_leveldata[lev]->vel_g, l_dt/2.0, *divtau_old[lev], 0, 0, 3, 0);

    // *************************************************************************************
    // Add the drag term implicitly
    // *************************************************************************************
    if (DEM::solve)
        mfix_add_drag_implicit(l_dt);

    // *************************************************************************************
    // Solve for u^star s.t. u^star = u_go + dt/2 (R_u^* + R_u^n) + dt/2 (Lu)^n + dt/2 (Lu)^star
    // Note we multiply ep_g by ro_g so that we pass in a single array holding (ro_g * ep_g)
    // *************************************************************************************
    for (int lev = 0; lev < nlev; lev++)
        MultiFab::Multiply(*m_leveldata[lev]->ep_g,
                           *m_leveldata[lev]->ro_g, 0, 0, 1,
                            m_leveldata[lev]->ep_g->nGrow());

    diffusion_op->diffuse_velocity(get_vel_g(), get_ep_g(), get_mu_g(), 0.5*l_dt);

    // mfix_set_tracer_bcs (new_time, trac, 0);
    diffusion_op->diffuse_scalar(get_trac(), get_ep_g(), mu_s, l_dt);

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

    mfix_apply_nodal_projection(depdt, new_time, l_dt, proj_2);
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
