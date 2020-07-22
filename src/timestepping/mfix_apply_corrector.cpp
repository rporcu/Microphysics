#include <mfix.H>

#include <AMReX_VisMF.H>
#include <mfix_mf_helpers.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

//
// Compute corrector:
//
//  1. Compute
//
//     vel_g = vel_go + dt * (R_u^* + R_u^n) / 2 + dt * ( g - grad(p_g+p0)/ro_g )
//
//     where the starred variables are computed using "predictor-step" variables.
//
//  2. Add implicit forcing term ( AKA implicit part of particles
//     momentum exchange )
//
//     vel_g = (vel_g + (drag_coeff*velp)/(ro_g*ep_g) / ( 1 + dt * drag_coeff/(ro_g*ep_g)
//
//  3. Add Crank-Nicolson diffusion
//
//     vel_g^* = vel_g + dt/2 * ( divtau^n * (1/(ro_g*ep_g)) + divtau^* * (1/(ro_g*ep_g)) )
//
//  3. Solve for phi
//
//     div( ep_g * grad(phi) / ro_g ) = div( ep_g * vel_g^* / dt + grad(p_g)/ro_g )
//
//  4. Compute
//
//     vel_g^{n+1} = vel_g^* -  dt * grad(phi) / ro_g
//
//  5. Define
//
//     p_g = phi
//
void
mfix::mfix_apply_corrector (Vector< MultiFab* >& conv_u_old,
                            Vector< MultiFab* >& conv_s_old,
                            Vector< MultiFab* >& conv_X_old,
                            Vector< MultiFab* >& divtau_old,
                            Vector< MultiFab* >&   laps_old,
                            Vector< MultiFab* >&   lapT_old,
                            Vector< MultiFab* >&   lapX_old,
                            Real time,
                            Real l_dt,
                            Real l_prev_dt,
                            bool proj_2)
{
    BL_PROFILE("mfix::mfix_apply_corrector");

    // We use the new-time value for things computed on the "*" state
    Real new_time = time + l_dt;

    // *************************************************************************************
    // Allocate space for half-time density and convective terms
    // *************************************************************************************
    Vector<MultiFab> density_nph;
    Vector<MultiFab*> conv_u;
    Vector<MultiFab*> conv_s;
    Vector<MultiFab*> conv_X;

    conv_u.resize(finest_level+1);
    conv_s.resize(finest_level+1);
    conv_X.resize(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        density_nph.emplace_back(grids[lev], dmap[lev],       1, 1, MFInfo(),  *ebfactory[lev]);

        conv_u[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
        // 3 components: one for density, one for tracer and one for enthalpy
        conv_s[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);

        conv_u[lev]->setVal(0.0);
        conv_s[lev]->setVal(0.0);

      if (advect_fluid_species) {
        conv_X[lev] = new MultiFab(grids[lev], dmap[lev], FLUID::nspecies_g, 0,
            MFInfo(), *ebfactory[lev]);

        conv_X[lev]->setVal(0.0);
      }
    }

    // *************************************************************************************
    // Compute the explicit advective term R_u^*
    // *************************************************************************************
    mfix_compute_convective_term(conv_u, conv_s, conv_X, get_vel_g(),
        get_ep_g(), get_ro_g(), get_h_g(), get_trac(), get_X_g(), new_time);

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
                  int conv_comp = 0;

                  const Real epg_loc = epg(i,j,k);
                  const Real rho_o_loc = rho_o(i,j,k);

                  Real rho = epg_loc*rho_o_loc + .5*l_dt*(drdt_o(i,j,k,conv_comp) + drdt(i,j,k,conv_comp));
                  rho /= epg_loc;

                  rho_new(i,j,k) = rho;

                  rho_nph(i,j,k) = 0.5 * (rho_o_loc + rho);
                });
            } // mfi
        } // lev

    } // not constant density

    // *************************************************************************************
    // Update enthalpy and temperature
    // *************************************************************************************
    if (advect_enthalpy)
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
                Array4<Real const> const& h_g_o  = ld.h_go->const_array(mfi);
                Array4<Real      > const& h_g_n  = ld.h_g->array(mfi);
                Array4<Real      > const& T_g_n  = ld.T_g->array(mfi);
                Array4<Real const> const& rho_o  = ld.ro_go->const_array(mfi);
                Array4<Real const> const& rho_n  = ld.ro_g->const_array(mfi);
                Array4<Real const> const& epg    = ld.ep_g->array(mfi);
                Array4<Real const> const& cp_g   = ld.cp_g->array(mfi);
                Array4<Real const> const& dhdt_o = conv_s_old[lev]->const_array(mfi);
                Array4<Real const> const& dhdt   = conv_s[lev]->const_array(mfi);
                Array4<Real const> const& lapT_o = lapT_old[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                  int conv_comp = 1;

                  const Real num =         (rho_o(i,j,k) * epg(i,j,k));
                  const Real denom = 1.0 / (rho_n(i,j,k) * epg(i,j,k));

                  // Crank-Nicolson so we only add half of the diffusive term here
                  Real h_g = num * h_g_o(i,j,k) + 0.5 * l_dt * (dhdt_o(i,j,k,conv_comp) + dhdt(i,j,k,conv_comp))
                                                + 0.5 * l_dt * lapT_o(i,j,k);

                  h_g_n(i,j,k) = h_g * denom;

                  T_g_n(i,j,k) = h_g_n(i,j,k) / cp_g(i,j,k);

                });
            } // mfi
        } // lev
    } // advect_enthalpy

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
                Array4<Real const> const& laps_o = laps_old[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                  const Real epg_loc = epg(i,j,k);

                  for (int n = 0; n < l_ntrac; ++n)
                  {
                    int conv_comp = 2+n;

                    Real tra = (rho_o(i,j,k)*epg_loc)*tra_o(i,j,k,n)
                             + 0.5 * l_dt * (dtdt_o(i,j,k,conv_comp) + dtdt(i,j,k,conv_comp));
                    tra /= (rho_n(i,j,k)*epg_loc);

                    // Crank-Nicolson so we add the explicit half here
                    tra += 0.5 * laps_o(i,j,k,n);

                    tra_n(i,j,k,n) = tra;
                  }
                });
            } // mfi
        } // lev
    } // advect_tracer

    // *************************************************************************
    // Update species mass fraction
    // *************************************************************************
    if (advect_fluid_species)
    {
      const int nspecies_g = FLUID::nspecies_g;

      for (int lev = 0; lev <= finest_level; lev++)
      {
        auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*ld.vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.tilebox();
          Array4<Real const> const& X_g_o  = ld.X_go->const_array(mfi);
          Array4<Real      > const& X_g_n  = ld.X_g->array(mfi);
          Array4<Real const> const& rho_o  = ld.ro_go->const_array(mfi);
          Array4<Real const> const& rho_n  = ld.ro_g->const_array(mfi);
          Array4<Real      > const& epg    = ld.ep_g->array(mfi);
          Array4<Real const> const& dXdt_o = conv_X_old[lev]->const_array(mfi);
          Array4<Real const> const& dXdt   = conv_X[lev]->const_array(mfi);
          Array4<Real const> const& lapX_o = lapX_old[lev]->const_array(mfi);

          ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const Real epg_loc = epg(i,j,k);

            for (int n = 0; n < nspecies_g; ++n)
            {
              int conv_comp = n;

              const Real num =         (rho_o(i,j,k) * epg_loc);
              const Real denom = 1.0 / (rho_n(i,j,k) * epg_loc);

              // Crank-Nicolson so we only add half of the diffusive term here
              Real X_g = num * X_g_o(i,j,k,n)
                  + 0.5 * l_dt * (dXdt_o(i,j,k,conv_comp) + dXdt(i,j,k,conv_comp))
                  + 0.5 * l_dt * lapX_o(i,j,k,n);

              X_g_n(i,j,k,n) = X_g * denom;
            }
          });
        } // mfi
      } // lev
    } // advect_fluid_species

    // *************************************************************************
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
         Array4<Real const> const& dudt_o   = conv_u_old[lev]->const_array(mfi);
         Array4<Real const> const& dudt     = conv_u[lev]->const_array(mfi);
         Array4<Real const> const& gp       = ld.gp->const_array(mfi);
         Array4<Real const> const& rho_nph  = density_nph[lev].const_array(mfi);
         Array4<Real const> const& epg      = ld.ep_g->const_array(mfi);
         Array4<Real const> const& divtau_o = divtau_old[lev]->const_array(mfi);

         // We need this until we remove static attribute from mfix::gravity
         const RealVect gp0_dev(gp0);
         const RealVect gravity_dev(gravity);

         amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
           const Real epg_loc = epg(i,j,k);

           Real vel_nx = epg_loc*vel_o(i,j,k,0) + .5*l_dt*(dudt_o(i,j,k,0)+dudt(i,j,k,0));
           Real vel_ny = epg_loc*vel_o(i,j,k,1) + .5*l_dt*(dudt_o(i,j,k,1)+dudt(i,j,k,1));
           Real vel_nz = epg_loc*vel_o(i,j,k,2) + .5*l_dt*(dudt_o(i,j,k,2)+dudt(i,j,k,2));

           vel_nx /= epg_loc;
           vel_ny /= epg_loc;
           vel_nz /= epg_loc;

           // Crank-Nicolson so we should only add half of the explicit term here, but
           //     we go ahead and add all of it now before doing the implicit drag solve,
           //     then we will subtract half of it after the drag solve
           vel_nx += l_dt * divtau_o(i,j,k,0);
           vel_ny += l_dt * divtau_o(i,j,k,1);
           vel_nz += l_dt * divtau_o(i,j,k,2);

           Real inv_dens = 1.0 / rho_nph(i,j,k);
           vel_nx += l_dt * (gravity_dev[0]-(gp(i,j,k,0)+gp0_dev[0])*inv_dens);
           vel_ny += l_dt * (gravity_dev[1]-(gp(i,j,k,1)+gp0_dev[1])*inv_dens);
           vel_nz += l_dt * (gravity_dev[2]-(gp(i,j,k,2)+gp0_dev[2])*inv_dens);

           vel_n(i,j,k,0) = vel_nx;
           vel_n(i,j,k,1) = vel_ny;
           vel_n(i,j,k,2) = vel_nz;
         });
       }
    }

    // *************************************************************************************
    // Add the drag term implicitly
    // *************************************************************************************
    if (DEM::solve or PIC::solve)
        mfix_add_txfr_implicit(l_dt);

    // *************************************************************************************
    // Subtract off half of the explicit diffusion term (see comment above)
    // *************************************************************************************
    for (int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab::Saxpy(*m_leveldata[lev]->vel_g, -l_dt/2.0, *divtau_old[lev], 0, 0,     3, 0);
    }

    // *************************************************************************************
    // Solve for u^star s.t. u^star = u_go + dt/2 (R_u^* + R_u^n) + dt/2 (Lu)^n + dt/2 (Lu)^star
    // Note we multiply ep_g by ro_g so that we pass in a single array holding (ro_g * ep_g)
    // *************************************************************************************

    // NOTE: we do this call before multiplying ep_g by ro_g
    if (advect_enthalpy) {
      diffusion_op->diffuse_temperature(get_T_g(), get_ep_g(), get_ro_g(),
          get_h_g(), get_cp_g(), get_k_g(), get_T_g_on_eb(), get_k_g_on_eb(),
          0.5*l_dt);
    }

    // Convert "ep_g" into (rho * ep_g)
    for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Multiply(*m_leveldata[lev]->ep_g,
                           *m_leveldata[lev]->ro_g, 0, 0, 1,
                            m_leveldata[lev]->ep_g->nGrow());

    diffusion_op->diffuse_velocity(get_vel_g(), get_ep_g(), get_mu_g(), 0.5*l_dt);

    // mfix_set_tracer_bcs (new_time, get_trac(), 0);
    if (advect_tracer)
        diffusion_op->diffuse_scalar(get_trac(), get_ep_g(), mu_s, 0.5*l_dt);

    if (advect_fluid_species) {
      diffusion_op->diffuse_species(get_X_g(), get_ep_g(), get_D_g(), 0.5*l_dt);
    }

    // Convert (rho * ep_g) back into ep_g
    for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Divide(*m_leveldata[lev]->ep_g,
                         *m_leveldata[lev]->ro_g, 0, 0, 1,
                          m_leveldata[lev]->ep_g->nGrow());

    // *************************************************************************************
    // Apply projection
    // *************************************************************************************
    Vector< MultiFab* > depdt(finest_level+1);
    for (int lev(0); lev <= finest_level; ++lev)
        depdt[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0, 1).release();

    mfix_apply_nodal_projection(depdt, new_time, l_dt, l_prev_dt, proj_2);

    for (int lev(0); lev <= finest_level; ++lev)
      delete depdt[lev];

    // *************************************************************************************
    // Correct small cells
    // *************************************************************************************
    mfix_correct_small_cells(get_vel_g());

    for (int lev = 0; lev <= finest_level; lev++)
    {
       delete conv_u[lev];
       delete conv_s[lev];

       if (advect_fluid_species)
         delete conv_X[lev];
    }
}
