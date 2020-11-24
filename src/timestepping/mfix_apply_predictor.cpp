#include <mfix.H>

#include <AMReX_VisMF.H>
#include <mfix_mf_helpers.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_calc_fluid_coeffs.H>

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
                            Vector< MultiFab* >& conv_X_old,
                            Vector< MultiFab* >& ro_RHS_old,
                            Vector< MultiFab* >& divtau_old,
                            Vector< MultiFab* >& lap_trac_old,
                            Vector< MultiFab* >& enthalpy_RHS_old,
                            Vector< MultiFab* >& lap_T_old,
                            Vector< MultiFab* >& lap_T_star,
                            Vector< MultiFab* >& species_RHS_old,
                            Vector< MultiFab* >& lap_X_old,
                            Vector< MultiFab* >& lap_X_star,
                            Real time,
                            Real l_dt,
                            Real l_prev_dt,
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
    // Note that "conv_s_old" returns update to (ep_g rho), (ep_g rho h_g) and (ep_g rho tracer)
    // *************************************************************************************
    bool update_laplacians = open_system_constraint;

    mfix_compute_convective_term(update_laplacians, conv_u_old, conv_s_old,
        conv_X_old, lap_T_old, lap_X_old, get_vel_g_old(), get_ep_g(),
        get_ro_g_old(), get_MW_g(), get_T_g_old(), get_cp_g(), get_k_g(),
        get_h_g_old(), get_T_g_on_eb(), get_k_g_on_eb(), get_trac_old(),
        get_X_gk_old(), get_D_gk(), get_h_gk(), get_txfr(), get_ro_gk_txfr(),
        time);

    // *************************************************************************************
    // Compute right hand side terms on the old status
    // *************************************************************************************
    bool explicit_diffusion_pred = true;

    {
      for (int lev = 0; lev <= finest_level; lev++)
          ro_RHS_old[lev]->setVal(0.);

      mfix_density_rhs(ro_RHS_old, get_ro_gk_txfr());
    }

    if (explicit_diffusion_pred)
    {
      diffusion_op->ComputeDivTau(divtau_old, get_vel_g_old(), get_ro_g_old(),
          get_ep_g(), get_mu_g());

      for (int lev = 0; lev <= finest_level; lev++)
        EB_set_covered(*divtau_old[lev], 0, divtau_old[lev]->nComp(),
            divtau_old[lev]->nGrow(), 0.);
    }
    else {
      for (int lev = 0; lev <= finest_level; lev++)
          divtau_old[lev]->setVal(0.);
    }

    if (advect_enthalpy) {
      for (int lev = 0; lev <= finest_level; lev++)
          enthalpy_RHS_old[lev]->setVal(0.);

      bool update_laplacian = explicit_diffusion_pred and
                              (not open_system_constraint);

      mfix_enthalpy_rhs(update_laplacian, enthalpy_RHS_old, lap_T_old,
          get_T_g_old(), get_ep_g(), get_ro_g_old(), get_k_g(), get_T_g_on_eb(),
          get_k_g_on_eb(), get_X_gk_old(), get_D_gk(), get_h_gk());
    }

    if (advect_tracer) {
      for (int lev = 0; lev <= finest_level; lev++)
          lap_trac_old[lev]->setVal(0.);

      mfix_scalar_rhs(explicit_diffusion_pred, lap_trac_old, get_trac_old(),
          get_ep_g(), get_ro_g_old(), mu_s);
    }

    // Species
    if (advect_fluid_species) {
      for (int lev = 0; lev <= finest_level; lev++)
        species_RHS_old[lev]->setVal(0.);

      bool update_laplacian = explicit_diffusion_pred and
                              (not open_system_constraint);

      mfix_species_X_rhs(update_laplacian, species_RHS_old, lap_X_old,
          get_X_gk_old(), get_ep_g(), get_ro_g_old(), get_D_gk(),
          get_ro_gk_txfr());
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
                Array4<Real const> const& rho_RHS = ro_RHS_old[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [rho_new,rho_nph,rho_o,epg,drdt_o,l_dt,rho_RHS]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                  int conv_comp = 0;

                  const Real epg_loc = epg(i,j,k);
                  Real rho = epg_loc*rho_o(i,j,k) + l_dt * drdt_o(i,j,k,conv_comp)
                                                  + l_dt * rho_RHS(i,j,k);

                  rho /= epg_loc;

                  rho_new(i,j,k) = rho;
                  rho_nph(i,j,k) = 0.5 * (rho_o(i,j,k) + rho);
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
                Array4<Real const> const& h_g_o   = ld.h_go->const_array(mfi);
                Array4<Real      > const& h_g_n   = ld.h_g->array(mfi);
                Array4<Real      > const& T_g_n   = ld.T_g->array(mfi);
                Array4<Real const> const& rho_o   = ld.ro_go->const_array(mfi);
                Array4<Real const> const& rho_n   = ld.ro_g->const_array(mfi);
                Array4<Real const> const& lap_T_o = lap_T_old[lev]->const_array(mfi);
                Array4<Real const> const& h_RHS_o = enthalpy_RHS_old[lev]->const_array(mfi);
                Array4<Real const> const& epg     = ld.ep_g->const_array(mfi);
                Array4<Real const> const& cp_g    = ld.cp_g->const_array(mfi);
                Array4<Real const> const& dhdt_o  = conv_s_old[lev]->const_array(mfi);

                // explicit_diffusion_pred is handled inside RHS computation
                // no need to separate computation in here anymore
                amrex::ParallelFor(bx, [h_g_o,h_g_n,T_g_n,rho_o,rho_n,h_RHS_o,
                    epg,cp_g,dhdt_o,l_dt,lap_T_o,explicit_diffusion_pred]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    int conv_comp = 1;

                    const Real num =         (rho_o(i,j,k) * epg(i,j,k));
                    const Real denom = 1.0 / (rho_n(i,j,k) * epg(i,j,k));

                    Real h_g = num * h_g_o(i,j,k) + l_dt * dhdt_o(i,j,k,conv_comp)
                                                  + l_dt * h_RHS_o(i,j,k);

                    if (explicit_diffusion_pred)
                      h_g += l_dt * lap_T_o(i,j,k);

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
                Array4<Real      > const& tra_n     = ld.trac->array(mfi);
                Array4<Real const> const& tra_o     = ld.trac_o->const_array(mfi);
                Array4<Real const> const& rho_o     = ld.ro_go->const_array(mfi);
                Array4<Real const> const& rho_n     = ld.ro_g->const_array(mfi);
                Array4<Real const> const& lap_tra_o = lap_trac_old[lev]->const_array(mfi);
                Array4<Real const> const& epg       = ld.ep_g->const_array(mfi);
                Array4<Real const> const& dtdt_o    = conv_s_old[lev]->const_array(mfi);

                // explicit_diffusion_pred is handled inside RHS computation
                // no need to separate computation in here anymore
                amrex::ParallelFor(bx, [tra_n,tra_o,rho_o,rho_n,lap_tra_o,epg,
                    dtdt_o,l_ntrac,l_dt,explicit_diffusion_pred]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                  for (int n = 0; n < l_ntrac; ++n)
                  {
                    int conv_comp = 2+n;

                    const Real epg_loc = epg(i,j,k);

                    Real tra = (rho_o(i,j,k)*epg_loc)*tra_o(i,j,k,n);
                    tra += l_dt * dtdt_o(i,j,k,conv_comp);

                    if (explicit_diffusion_pred)
                      tra += l_dt * lap_tra_o(i,j,k,n);

                    tra /= (rho_n(i,j,k)*epg_loc);

                    tra_n(i,j,k,n) = tra;
                  }
                });
            } // mfi
        } // lev
    } // advect_tracer

    // *************************************************************************
    // Update species mass fractions
    // *************************************************************************
    if (advect_fluid_species)
    {
      const int nspecies_g = FLUID::nspecies;

      for (int lev = 0; lev <= finest_level; lev++)
      {
        auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*ld.vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.tilebox();
          Array4<Real      > const& X_gk_n  = ld.X_gk->array(mfi);
          Array4<Real const> const& X_gk_o  = ld.X_gko->const_array(mfi);
          Array4<Real const> const& rho_n   = ld.ro_g->const_array(mfi);
          Array4<Real const> const& rho_o   = ld.ro_go->const_array(mfi);
          Array4<Real const> const& epg     = ld.ep_g->const_array(mfi);
          Array4<Real const> const& dXdt_o  = conv_X_old[lev]->const_array(mfi);
          Array4<Real const> const& lap_X_o = lap_X_old[lev]->const_array(mfi);
          Array4<Real const> const& X_RHS_o = species_RHS_old[lev]->const_array(mfi);
          
          // explicit_diffusion_pred is handled inside RHS computation
          // no need to separate computation in here anymore
          ParallelFor(bx, [nspecies_g,epg,rho_o,rho_n,X_gk_o,dXdt_o,lap_X_o,
              l_dt,X_gk_n,X_RHS_o,explicit_diffusion_pred]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const Real num =         (rho_o(i,j,k) * epg(i,j,k));
            const Real denom = 1.0 / (rho_n(i,j,k) * epg(i,j,k));

            for (int n = 0; n < nspecies_g; ++n)
            {
              Real X_gk = num * X_gk_o(i,j,k,n) + l_dt * dXdt_o(i,j,k,n)
                                                + l_dt * X_RHS_o(i,j,k,n);

              if (explicit_diffusion_pred)
                X_gk += l_dt * lap_X_o(i,j,k,n);

              X_gk_n(i,j,k,n) = X_gk * denom;
            }
          });
        } // mfi
      } // lev
    } // advect_fluid_species

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

         amrex::ParallelFor(bx, [epg,vel_o,dudt_o,divtau_o,gp0_dev,gravity_dev,
             gp,l_dt,vel_n,rho_nph,explicit_diffusion_pred]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
           const Real epg_loc = epg(i,j,k);

           Real vel_nx = epg_loc*vel_o(i,j,k,0) + l_dt*dudt_o(i,j,k,0);
           Real vel_ny = epg_loc*vel_o(i,j,k,1) + l_dt*dudt_o(i,j,k,1);
           Real vel_nz = epg_loc*vel_o(i,j,k,2) + l_dt*dudt_o(i,j,k,2);

           vel_nx /= epg_loc;
           vel_ny /= epg_loc;
           vel_nz /= epg_loc;

           if (explicit_diffusion_pred)
           {
             vel_nx += l_dt * divtau_o(i,j,k,0);
             vel_ny += l_dt * divtau_o(i,j,k,1);
             vel_nz += l_dt * divtau_o(i,j,k,2);
           }

           Real inv_dens = 1.0 / rho_nph(i,j,k);
           vel_nx += l_dt * (gravity_dev[0]-(gp(i,j,k,0)+gp0_dev[0])*inv_dens);
           vel_ny += l_dt * (gravity_dev[1]-(gp(i,j,k,1)+gp0_dev[1])*inv_dens);
           vel_nz += l_dt * (gravity_dev[2]-(gp(i,j,k,2)+gp0_dev[2])*inv_dens);

           vel_n(i,j,k,0) = vel_nx;
           vel_n(i,j,k,1) = vel_ny;
           vel_n(i,j,k,2) = vel_nz;
         });

       } // mfi
    } // lev

    // *************************************************************************************
    // Add the drag term implicitly
    // *************************************************************************************
    if (DEM::solve or PIC::solve)
        mfix_add_txfr_implicit(l_dt);

    // *************************************************************************************
    // If doing implicit diffusion, solve here for u^*
    // Note we multiply ep_g by ro_g so that we pass in a single array holding (ro_g * ep_g)
    // *************************************************************************************
    if (not explicit_diffusion_pred)
    {
      mfix_set_density_bcs(time, get_ro_g());
      mfix_set_tracer_bcs(time, get_trac());

      if (advect_enthalpy)
        mfix_set_temperature_bcs(time, get_T_g());

      mfix_set_scalar_bcs(time, get_mu_g(), get_cp_g(), get_k_g(), get_MW_g());

      if (advect_enthalpy)
        mfix_set_enthalpy_bcs(time, get_h_g());

      if (advect_fluid_species)
        mfix_set_species_bcs(time, get_X_gk(), get_D_gk(), get_cp_gk(), get_h_gk());

      // NOTE: we do this call before multiplying ep_g by ro_g
      // Diffuse enthalpy
      if (advect_enthalpy) {
        diffusion_op->diffuse_temperature(get_T_g(), get_ep_g(), get_ro_g(), get_h_g(),
            get_cp_g(), get_k_g(), get_T_g_on_eb(), get_k_g_on_eb(), l_dt);
      }

      // Convert "ep_g" into (rho * ep_g)
      for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Multiply(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                           0, 0, 1, m_leveldata[lev]->ep_g->nGrow());

      // Set velocity boundary conditions
      mfix_set_velocity_bcs(new_time, get_vel_g(), 0);
      // Diffuse velocity
      diffusion_op->diffuse_velocity(get_vel_g(), get_ep_g(), get_mu_g(), l_dt);

      // Diffuse tracer
      if (advect_tracer) {
        diffusion_op->diffuse_scalar(get_trac(), get_ep_g(), mu_s, l_dt);
      }

      // Diffuse species mass fractions
      if (advect_fluid_species) {
        diffusion_op->diffuse_species(get_X_gk(), get_ep_g(), get_D_gk(), l_dt);
      }

      // Convert (rho * ep_g) back into ep_g
      for (int lev = 0; lev <= finest_level; lev++)
          MultiFab::Divide(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                           0, 0, 1, m_leveldata[lev]->ep_g->nGrow());
    }

    // *************************************************************************************
    // Rescale species in order to respect sum = 1
    // *************************************************************************************
    if (advect_fluid_species) {
      mfix_normalize_fluid_species(get_X_gk());
    }

    // *************************************************************************************
    // Update fluid (if fluid is a mixture) and fluid species specific heat,
    // enthalpy and temperature
    // *************************************************************************************
    if (advect_fluid_species) {
      mfix_update_fluid_and_species(get_cp_gk(), get_h_gk(), get_MW_g(),
          get_cp_g(), get_h_g(), get_T_g(), get_X_gk());
    }

    // *************************************************************************************
    // Project velocity field -- depdt=0 for now
    // *************************************************************************************
    Vector< MultiFab* > depdt(finest_level+1);
    Vector< MultiFab* > constraint_RHS(finest_level+1);
    Vector< MultiFab* > S_cc(finest_level+1);

    for (int lev(0); lev <= finest_level; ++lev) {
      depdt[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0, 1).release();
      constraint_RHS[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0, 1).release();
      S_cc[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0, 1).release();
    }

    if (open_system_constraint) {
      bool update_laplacians = true;

      mfix_open_system_rhs(constraint_RHS, update_laplacians, lap_T_star,
          lap_X_star, get_ep_g(), get_ro_g(), get_MW_g(), get_T_g(), get_cp_g(),
          get_k_g(), get_T_g_on_eb(), get_k_g_on_eb(), get_X_gk(), get_D_gk(),
          get_h_gk(), get_txfr(), get_ro_gk_txfr());
    }

    for (int lev(0); lev <= finest_level; ++lev) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*S_cc[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        Array4< Real > const& depdt_array = depdt[lev]->array(mfi);
        Array4< Real > const& constraint_RHS_array = constraint_RHS[lev]->array(mfi);
        Array4< Real > const& S_cc_array = S_cc[lev]->array(mfi);

        amrex::ParallelFor(bx, [S_cc_array,depdt_array,constraint_RHS_array]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          S_cc_array(i,j,k) = constraint_RHS_array(i,j,k) - depdt_array(i,j,k);
        });
      }
    }

    mfix_apply_nodal_projection(S_cc, new_time, l_dt, l_prev_dt, proj_2);

    for (int lev(0); lev <= finest_level; ++lev) {
      delete depdt[lev];
      delete constraint_RHS[lev];
      delete S_cc[lev];
    }

    // *************************************************************************************
    // Correct small cells
    // *************************************************************************************
    mfix_correct_small_cells (get_vel_g());
}
