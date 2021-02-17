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
                            Vector< MultiFab* >& ro_RHS_old,
                            Vector< MultiFab* >& /*divtau_old*/,
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
    BL_PROFILE("mfix::mfix_apply_corrector");

    // We use the new-time value for things computed on the "*" state
    Real new_time = time + l_dt;

    fillpatch_all(get_vel_g(), get_ro_g(), get_h_g(), get_trac(),
                  get_X_gk(), new_time);

    mfix_set_density_bcs(time, get_ro_g());

    // *************************************************************************************
    // Allocate space for the forcing terms
    // *************************************************************************************
    // Forcing terms
    Vector<MultiFab> vel_forces, tra_forces;
    Vector<MultiFab> vel_eta, tra_eta;

    for (int lev = 0; lev <= finest_level; ++lev) {

      vel_forces.emplace_back(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost_force(),
                              MFInfo(), EBFactory(lev));

      vel_eta.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(), EBFactory(lev));

      if (advect_tracer) {
        tra_forces.emplace_back(grids[lev], dmap[lev], ntrac, nghost_force(),
                                MFInfo(), EBFactory(lev));

        tra_eta.emplace_back(grids[lev], dmap[lev], ntrac, 1, MFInfo(), EBFactory(lev));
      }
    }

    // *************************************************************************************
    // Allocate space for half-time density and convective terms
    // *************************************************************************************
    Vector<MultiFab> density_nph;
    Vector<MultiFab*> conv_u(finest_level+1);
    Vector<MultiFab*> conv_s(finest_level+1);
    Vector<MultiFab*> conv_X(finest_level+1);
    Vector<MultiFab*> ro_RHS(finest_level+1);
    Vector<MultiFab*> lap_trac(finest_level+1);
    Vector<MultiFab*> enthalpy_RHS(finest_level+1);
    Vector<MultiFab*> species_RHS(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        density_nph.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(),  *ebfactory[lev]);

        conv_u[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
        // 3 components: one for density, one for tracer and one for enthalpy
        conv_s[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);

        ro_RHS[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
        lap_trac[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), *ebfactory[lev]);
        enthalpy_RHS[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);

        conv_u[lev]->setVal(0.0);
        conv_s[lev]->setVal(0.0);

        ro_RHS[lev]->setVal(0.0);
        lap_trac[lev]->setVal(0.0);
        enthalpy_RHS[lev]->setVal(0.0);

      if (advect_fluid_species) {
        conv_X[lev] = new MultiFab(grids[lev], dmap[lev], FLUID::nspecies, 0,
            MFInfo(), *ebfactory[lev]);
        species_RHS[lev] = new MultiFab(grids[lev], dmap[lev], FLUID::nspecies,
            0, MFInfo(), *ebfactory[lev]);

        conv_X[lev]->setVal(0.0);
        species_RHS[lev]->setVal(0.0);
      }
    }


    // *************************************************************************************
    // Allocate space for the MAC velocities
    // *************************************************************************************
    Vector<MultiFab> ep_u_mac(finest_level+1), ep_v_mac(finest_level+1), ep_w_mac(finest_level+1);
    Vector<MultiFab> rhs_mac(finest_level+1);

    int ngmac = nghost_mac();

    for (int lev = 0; lev <= finest_level; ++lev) {
      ep_u_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(0)),
                           dmap[lev], 1, ngmac, MFInfo(), *ebfactory[lev]);
      ep_v_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(1)),
                           dmap[lev], 1, ngmac, MFInfo(), *ebfactory[lev]);
      ep_w_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(2)),
                           dmap[lev], 1, ngmac, MFInfo(), *ebfactory[lev]);

      rhs_mac[lev].define(grids[lev], dmap[lev], 1, ngmac, MFInfo(), EBFactory(lev));
      rhs_mac[lev].setVal(0.);

      if (ngmac > 0) {
        ep_u_mac[lev].setBndry(0.0);
        ep_v_mac[lev].setBndry(0.0);
        ep_w_mac[lev].setBndry(0.0);
      }
    }


    // *************************************************************************************
    // Compute explicit diffusive terms
    // *************************************************************************************

    const bool explicit_diffusive_enthalpy = false;
    const bool explicit_diffusive_trac     = false;
    const bool explicit_diffusive_species  = false;

    {
      // We do not need to calculate the Laplacians if the open system constraint is used
      // because it they were previously computed at the end of the predictor step,
      // prior to the nodal projection.

      const bool update_lapT = (advect_enthalpy      && (explicit_diffusive_enthalpy && !open_system_constraint));
      const bool update_lapS = (advect_tracer        &&  explicit_diffusive_trac);
      const bool update_lapX = (advect_fluid_species && (explicit_diffusive_species  && !open_system_constraint));

      compute_laps(update_lapT, update_lapS, update_lapX, lap_T_old, lap_trac_old, lap_X_old,
                   get_T_g_old(), get_trac_old(), get_X_gk_old(),
                   get_ep_g_const(), get_ro_g_old_const());
    }

    // *************************************************************************************
    // Compute RHS for the MAC projection
    // *************************************************************************************

    if (open_system_constraint) {
      mfix_open_system_rhs(GetVecOfPtrs(rhs_mac), lap_T_star, lap_X_star,
         get_ep_g_const(), get_ro_g_const(), get_MW_g_const(),
         get_T_g_const(), get_cp_g_const(), get_X_gk_const(), get_D_gk_const(),
         get_h_gk_const(), get_txfr_const(), get_ro_gk_txfr_const());
    }

    // *************************************************************************************
    // Compute the explicit advective term R_u^*
    // *************************************************************************************

    mfix_compute_convective_term(conv_u, conv_s, conv_X,
        GetVecOfPtrs(vel_forces), GetVecOfPtrs(tra_forces),
        get_vel_g_const(), get_ep_g_const(), get_ro_g_const(),
        get_h_g_const(), get_trac_const(), get_X_gk_const(),
        GetVecOfPtrs(ep_u_mac), GetVecOfPtrs(ep_v_mac), GetVecOfPtrs(ep_w_mac),
        GetVecOfConstPtrs(rhs_mac), get_divtau(), l_dt, new_time);

    // *************************************************************************************
    // Compute right hand side terms on the intermediate status
    // *************************************************************************************

    {
      mfix_density_rhs(ro_RHS, get_ro_gk_txfr_const());
    }

    if (advect_enthalpy) {

      mfix_enthalpy_rhs(enthalpy_RHS, get_ep_g_const(), get_ro_g_const(),
           get_X_gk_const(), get_D_gk_const(), get_h_gk_const());
    }

    if (advect_tracer) {
      for (int lev = 0; lev <= finest_level; lev++)
          lap_trac[lev]->setVal(0.);

      mfix_scalar_rhs(explicit_diffusive_trac, lap_trac, get_trac(), get_ep_g(),
          get_ro_g(), mu_s);
    }

    if (advect_fluid_species) {

      mfix_species_X_rhs(species_RHS, get_ro_gk_txfr_const());
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
                Array4<Real const> const& rho_o     = ld.ro_go->const_array(mfi);
                Array4<Real      > const& rho_new   = ld.ro_g->array(mfi);
                Array4<Real      > const& rho_nph   = density_nph[lev].array(mfi);
                Array4<Real      > const& epg       = ld.ep_g->array(mfi);
                Array4<Real const> const& drdt_o    = conv_s_old[lev]->const_array(mfi);
                Array4<Real const> const& drdt      = conv_s[lev]->const_array(mfi);
                Array4<Real const> const& rho_rhs_o = ro_RHS_old[lev]->const_array(mfi);
                Array4<Real const> const& rho_rhs   = ro_RHS[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [epg,rho_o,l_dt,drdt,drdt_o,rho_new,
                    rho_nph,rho_rhs_o,rho_rhs]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                  int conv_comp = 0;

                  const Real epg_loc = epg(i,j,k);
                  const Real rho_o_loc = rho_o(i,j,k);

                  Real rho = epg_loc*rho_o_loc;
                  rho += .5*l_dt*(drdt_o(i,j,k,conv_comp)+drdt(i,j,k,conv_comp));
                  rho += .5*l_dt*(rho_rhs_o(i,j,k)+rho_rhs(i,j,k));

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
                Array4<Real const> const& h_g_o   = ld.h_go->const_array(mfi);
                Array4<Real      > const& h_g_n   = ld.h_g->array(mfi);
                Array4<Real      > const& T_g_n   = ld.T_g->array(mfi);
                Array4<Real const> const& rho_o   = ld.ro_go->const_array(mfi);
                Array4<Real const> const& rho_n   = ld.ro_g->const_array(mfi);
                Array4<Real const> const& epg     = ld.ep_g->array(mfi);
                Array4<Real const> const& cp_g    = ld.cp_g->array(mfi);
                Array4<Real const> const& dhdt_o  = conv_s_old[lev]->const_array(mfi);
                Array4<Real const> const& dhdt    = conv_s[lev]->const_array(mfi);
                Array4<Real const> const& h_rhs_o = enthalpy_RHS_old[lev]->const_array(mfi);
                Array4<Real const> const& h_rhs   = enthalpy_RHS[lev]->const_array(mfi);
                Array4<Real const> const& lap_T_o = lap_T_old[lev]->const_array(mfi);
                Array4<Real const> const& lap_T   = lap_T_star[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [h_g_o,h_g_n,T_g_n,rho_o,rho_n,epg,cp_g,
                    dhdt_o,dhdt,h_rhs_o,h_rhs,l_dt,lap_T_o,lap_T,
                    explicit_diffusive_enthalpy]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                  int conv_comp = 1;

                  const Real num =         (rho_o(i,j,k) * epg(i,j,k));
                  const Real denom = 1.0 / (rho_n(i,j,k) * epg(i,j,k));

                  Real h_g = num*h_g_o(i,j,k);
                  h_g += .5*l_dt*(dhdt_o(i,j,k,conv_comp)+dhdt(i,j,k,conv_comp));
                  h_g += .5*l_dt*(h_rhs_o(i,j,k)+h_rhs(i,j,k));

                  if (explicit_diffusive_enthalpy) {
                    h_g += .5*l_dt*(lap_T_o(i,j,k)+lap_T(i,j,k));
                  }
                  else {
                    // Crank-Nicolson so we only add half of the diffusive term
                    // here, but we go ahead and add all of it now before doing
                    // the implicit enthalpy solve, then we will subtract half
                    // of it after the enthalpy solve
                    h_g += l_dt * lap_T_o(i,j,k);
                  }

                  h_g *= denom;
                  h_g_n(i,j,k) = h_g;
                  T_g_n(i,j,k) = h_g / cp_g(i,j,k);
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
                Array4<Real const> const& tra_o     = ld.trac_o->const_array(mfi);
                Array4<Real      > const& tra_n     = ld.trac->array(mfi);
                Array4<Real const> const& rho_o     = ld.ro_go->const_array(mfi);
                Array4<Real const> const& rho_n     = ld.ro_g->const_array(mfi);
                Array4<Real> const& epg             = ld.ep_g->array(mfi);
                Array4<Real const> const& dtdt_o    = conv_s_old[lev]->const_array(mfi);
                Array4<Real const> const& dtdt      = conv_s[lev]->const_array(mfi);
                Array4<Real const> const& lap_tra_o = lap_trac_old[lev]->const_array(mfi);
                Array4<Real const> const& lap_tra   = lap_trac[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [tra_o,tra_n,rho_o,rho_n,epg,dtdt_o,dtdt,
                    lap_tra_o,lap_tra,l_dt,l_ntrac,explicit_diffusive_trac]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                  const Real epg_loc = epg(i,j,k);

                  const Real num =         (rho_o(i,j,k) * epg_loc);
                  const Real denom = 1.0 / (rho_n(i,j,k) * epg_loc);

                  for (int n = 0; n < l_ntrac; ++n)
                  {
                    int conv_comp = 2+n;

                    Real tra = num*tra_o(i,j,k,n);
                    tra += .5*l_dt*(dtdt_o(i,j,k,conv_comp)+dtdt(i,j,k,conv_comp));

                    // Crank-Nicolson so we add the explicit half here
                    tra += .5*l_dt*lap_tra_o(i,j,k,n);

                    if (explicit_diffusive_trac)
                      tra += 0.5 * l_dt * lap_tra(i,j,k,n);

                    tra_n(i,j,k,n) = tra * denom;
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
          Array4<Real const> const& X_gk_o  = ld.X_gko->const_array(mfi);
          Array4<Real      > const& X_gk_n  = ld.X_gk->array(mfi);
          Array4<Real const> const& rho_o   = ld.ro_go->const_array(mfi);
          Array4<Real const> const& rho_n   = ld.ro_g->const_array(mfi);
          Array4<Real      > const& epg     = ld.ep_g->array(mfi);
          Array4<Real const> const& dXdt_o  = conv_X_old[lev]->const_array(mfi);
          Array4<Real const> const& dXdt    = conv_X[lev]->const_array(mfi);
          Array4<Real const> const& lap_X_o = lap_X_old[lev]->const_array(mfi);
          Array4<Real const> const& lap_X   = lap_X_star[lev]->const_array(mfi);
          Array4<Real const> const& X_rhs_o = species_RHS_old[lev]->const_array(mfi);
          Array4<Real const> const& X_rhs   = species_RHS[lev]->const_array(mfi);

          ParallelFor(bx, [X_gk_o,X_gk_n,rho_o,rho_n,epg,dXdt_o,dXdt,X_rhs_o,
              X_rhs,l_dt,nspecies_g,lap_X_o,lap_X,explicit_diffusive_species]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const Real epg_loc = epg(i,j,k);

            const Real num =         (rho_o(i,j,k) * epg_loc);
            const Real denom = 1.0 / (rho_n(i,j,k) * epg_loc);

            for (int n = 0; n < nspecies_g; ++n)
            {
              Real X_gk = num*X_gk_o(i,j,k,n);
              X_gk += .5*l_dt*(dXdt_o(i,j,k,n)+dXdt(i,j,k,n));
              X_gk += .5*l_dt*(X_rhs_o(i,j,k,n)+X_rhs(i,j,k,n));

              // Crank-Nicolson so we only add half of the diffusive term here
              X_gk += .5*l_dt*lap_X_o(i,j,k,n);

              if (explicit_diffusive_species)
                X_gk += .5*l_dt*lap_X(i,j,k,n);

              X_gk_n(i,j,k,n) = X_gk * denom;
            }
          });
        } // mfi
      } // lev
    } // advect_fluid_species



    // *************************************************************************************
    // Define (or if advection_type != "MOL", re-define) the forcing terms, without the
    //    viscous terms and using the half-time density
    // *************************************************************************************
    compute_vel_forces(GetVecOfPtrs(vel_forces), get_vel_g_const(),
                       GetVecOfConstPtrs(density_nph));


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
         Array4<Real const> const& epg      = ld.ep_g->const_array(mfi);
         Array4<Real const> const& divtau_o = ld.divtau_o->const_array(mfi);
         Array4<Real const> const& vel_f    = vel_forces[lev].const_array(mfi);

         // We need this until we remove static attribute from mfix::gravity
         const RealVect gp0_dev(gp0);
         const RealVect gravity_dev(gravity);

         amrex::ParallelFor(bx, [vel_n,vel_o,dudt_o,dudt,gp,vel_f,epg,divtau_o,
             gp0_dev,gravity_dev,l_dt]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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

           vel_nx += l_dt * vel_f(i,j,k,0);
           vel_ny += l_dt * vel_f(i,j,k,1);
           vel_nz += l_dt * vel_f(i,j,k,2);

           vel_n(i,j,k,0) = vel_nx;
           vel_n(i,j,k,1) = vel_ny;
           vel_n(i,j,k,2) = vel_nz;
         });
       }
    }

    // *************************************************************************************
    // Add the drag and enthalpy terms implicitly
    // *************************************************************************************
    if (DEM::solve or PIC::solve)
        mfix_add_txfr_implicit(l_dt);

    // *************************************************************************************
    // Subtract off half of the explicit diffusion terms (see comment above)
    // *************************************************************************************
    if (advect_enthalpy and (not explicit_diffusive_enthalpy))
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
              Array4<Real      > const& h_g_n   = ld.h_g->array(mfi);
              Array4<Real      > const& T_g_n   = ld.T_g->array(mfi);
              Array4<Real const> const& cp_g    = ld.cp_g->array(mfi);
              Array4<Real const> const& lap_T_o = lap_T_old[lev]->const_array(mfi);

              amrex::ParallelFor(bx, [h_g_n,T_g_n,cp_g,lap_T_o,l_dt]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                Real h_g = h_g_n(i,j,k);
                h_g -= .5 * l_dt * lap_T_o(i,j,k);

                h_g_n(i,j,k) = h_g;
                T_g_n(i,j,k) = h_g / cp_g(i,j,k);
              });
          } // mfi
      } // lev
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
      MultiFab::Saxpy(*m_leveldata[lev]->vel_g, -l_dt/2.0, *m_leveldata[lev]->divtau_o, 0, 0, 3, 0);
    }

    // *************************************************************************************
    // Solve for u^star s.t. u^star = u_go + dt/2 (R_u^* + R_u^n) + dt/2 (Lu)^n + dt/2 (Lu)^star
    // Note we multiply ep_g by ro_g so that we pass in a single array holding (ro_g * ep_g)
    // *************************************************************************************

    // NOTE: we do this call before multiplying ep_g by ro_g
    if (advect_enthalpy and (not explicit_diffusive_enthalpy)) {
      diffusion_op->diffuse_temperature(get_T_g(), get_ep_g(), get_ro_g(), get_h_g(),
          get_cp_g(), get_k_g(), get_T_g_on_eb(), get_k_g_on_eb(), 0.5*l_dt);
    }

    // Convert "ep_g" into (rho * ep_g)
    for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Multiply(*m_leveldata[lev]->ep_g,
                           *m_leveldata[lev]->ro_g, 0, 0, 1,
                            m_leveldata[lev]->ep_g->nGrow());

    diffusion_op->diffuse_velocity(get_vel_g(), get_ep_g(), get_mu_g(), 0.5*l_dt);

    // mfix_set_tracer_bcs (new_time, get_trac(), 0);
    if (advect_tracer and (not explicit_diffusive_trac))
        diffusion_op->diffuse_scalar(get_trac(), get_ep_g(), mu_s, 0.5*l_dt);

    if (advect_fluid_species and (not explicit_diffusive_species)) {
      diffusion_op->diffuse_species(get_X_gk(), get_ep_g(), get_D_gk(), 0.5*l_dt);
    }

    // Convert (rho * ep_g) back into ep_g
    for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Divide(*m_leveldata[lev]->ep_g,
                         *m_leveldata[lev]->ro_g, 0, 0, 1,
                          m_leveldata[lev]->ep_g->nGrow());

    // *************************************************************************************
    // Normalize species mass fractions in order to respect sum = 1
    // *************************************************************************************
    if (advect_fluid_species) {
      mfix_normalize_fluid_species(get_X_gk());
    }

    // *************************************************************************************
    // Update fluid and fluid species specific heat, enthalpy and temperature
    // *************************************************************************************
    if (advect_fluid_species) {
      mfix_update_fluid_and_species(get_cp_gk(), get_h_gk(), get_MW_g(),
          get_cp_g(), get_h_g(), get_T_g(), get_X_gk());
    }

    // *************************************************************************************
    // Apply projection
    // *************************************************************************************
    Vector< MultiFab* > depdt(finest_level+1);
    Vector< MultiFab* > S_cc(finest_level+1);

    for (int lev(0); lev <= finest_level; ++lev) {
      depdt[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0, 1).release();
      S_cc[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0, 1).release();
    }

    if (open_system_constraint) {

      const bool update_lapT = advect_enthalpy;
      const bool update_lapS = false;
      const bool update_lapX = advect_fluid_species;

      compute_laps(update_lapT, update_lapS, update_lapX, lap_T_star, lap_trac, lap_X_star,
         get_T_g(), get_trac_old(), get_X_gk_old(), get_ep_g_const(), get_ro_g_old_const());

      mfix_open_system_rhs(S_cc, lap_T_star, lap_X_star, get_ep_g_const(), get_ro_g_const(),
          get_MW_g_const(), get_T_g_const(), get_cp_g_const(), get_X_gk_const(),
          get_D_gk_const(), get_h_gk_const(), get_txfr_const(), get_ro_gk_txfr_const());

    }


    if (m_use_depdt_constraint) {
      for (int lev(0); lev <= finest_level; ++lev) {
        MultiFab::Subtract(*S_cc[lev], *depdt[lev],0,0,1,0);
      }
    }

    mfix_apply_nodal_projection(S_cc, new_time, l_dt, l_prev_dt, proj_2,
                                GetVecOfConstPtrs(density_nph));

    for (int lev(0); lev <= finest_level; ++lev) {
      delete depdt[lev];
      delete S_cc[lev];
    }

    // *************************************************************************************
    // Correct small cells
    // *************************************************************************************
    mfix_correct_small_cells (get_vel_g(), GetVecOfConstPtrs(ep_u_mac),
         GetVecOfConstPtrs(ep_v_mac), GetVecOfConstPtrs(ep_w_mac));

    // *************************************************************************************
    // Free up memory from temporary data
    // *************************************************************************************
    for (int lev = 0; lev <= finest_level; lev++)
    {
       delete conv_u[lev];
       delete conv_s[lev];

       delete ro_RHS[lev];
       delete lap_trac[lev];
       delete enthalpy_RHS[lev];

       if (advect_fluid_species) {
         delete conv_X[lev];
         delete species_RHS[lev];
       }
    }
}
