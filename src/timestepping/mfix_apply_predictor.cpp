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
                            Vector< MultiFab* >& ro_RHS,
                            Vector< MultiFab* >& lap_trac_old,
                            Vector< MultiFab* >& lap_trac,
                            Vector< MultiFab* >& enthalpy_RHS_old,
                            Vector< MultiFab* >& enthalpy_RHS,
                            Vector< MultiFab* >& lap_T_old,
                            Vector< MultiFab* >& lap_T,
                            Vector< MultiFab* >& species_RHS_old,
                            Vector< MultiFab* >& species_RHS,
                            Vector< MultiFab* >& lap_X_old,
                            Vector< MultiFab* >& lap_X,
                            Vector< Real >& rhs_pressure_g_old,
                            Vector< Real >& rhs_pressure_g,
                            Real time,
                            Real l_dt,
                            Real l_prev_dt,
                            bool proj_2,
                            Real& coupling_timing)
{
    BL_PROFILE("mfix::mfix_apply_predictor");

    // We use the new-time value for things computed on the "*" state
    Real new_time = time + l_dt;

    // Averaged quantities for closed system constraint computations
    Vector<Real> avgSigma(finest_level+1, 0);
    Vector<Real> avgTheta(finest_level+1, 0);

    // Local flag for explicit diffusion
    bool l_explicit_diff = (predictor_diff_type() == DiffusionType::Explicit);

    mfix_set_density_bcs(time, get_ro_g_old());

    // *************************************************************************************
    // Allocate space for the forcing terms
    // *************************************************************************************
    // Forcing terms
    Vector<MultiFab> vel_forces, tra_forces;
    Vector<MultiFab> vel_eta, tra_eta;

    for (int lev = 0; lev <= finest_level; ++lev) {
      vel_forces.emplace_back(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost_force(),
                              MFInfo(), EBFactory(lev));
      vel_forces[lev].setVal(0.);

      if (advect_tracer) {
        tra_forces.emplace_back(grids[lev], dmap[lev], ntrac, nghost_force(),
                                MFInfo(), EBFactory(lev));
        tra_forces[lev].setVal(0.);
      }

      vel_eta.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(), EBFactory(lev));
      vel_eta[lev].setVal(0.);

      if (advect_tracer) {
        tra_eta.emplace_back(grids[lev], dmap[lev], ntrac, 1, MFInfo(), EBFactory(lev));
        tra_eta[lev].setVal(0.);
      }
    }

    // *************************************************************************************
    // Allocate space for half-time density
    // *************************************************************************************
    Vector<MultiFab> density_nph;
    for (int lev = 0; lev <= finest_level; ++lev) {
        density_nph.emplace_back(grids[lev], dmap[lev], 1, nghost_state(), MFInfo(),  *ebfactory[lev]);
        density_nph[lev].setVal(0.0);
    }



    // *************************************************************************************
    // Allocate space for the MAC velocities and RHS
    // *************************************************************************************
    Vector<MultiFab> ep_u_mac(finest_level+1), ep_v_mac(finest_level+1), ep_w_mac(finest_level+1);
    Vector<MultiFab> rhs_mac(finest_level+1),  depdt(finest_level+1);

    int ngmac = nghost_mac();

    for (int lev = 0; lev <= finest_level; ++lev) {
      ep_u_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(0)),
                           dmap[lev], 1, ngmac, MFInfo(), EBFactory(lev));
      ep_u_mac[lev].setVal(0.);

      ep_v_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(1)),
                           dmap[lev], 1, ngmac, MFInfo(), EBFactory(lev));
      ep_v_mac[lev].setVal(0.);

      ep_w_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(2)),
                           dmap[lev], 1, ngmac, MFInfo(), EBFactory(lev));
      ep_w_mac[lev].setVal(0.);

      if (ngmac > 0) {
        ep_u_mac[lev].setBndry(0.0);
        ep_v_mac[lev].setBndry(0.0);
        ep_w_mac[lev].setBndry(0.0);
      }

      rhs_mac[lev].define(grids[lev], dmap[lev], 1, ngmac, MFInfo(), EBFactory(lev));
      rhs_mac[lev].setVal(0.);

      if (m_use_depdt_constraint){
        depdt[lev].define(grids[lev], dmap[lev], 1, ngmac, MFInfo(), EBFactory(lev));
        depdt[lev].setVal(0.);
      }
    }


    // *************************************************************************************
    // Compute explicit diffusive terms
    // *************************************************************************************

    if (need_divtau()) {
      diffusion_op->ComputeDivTau(get_divtau(), get_vel_g_old(), get_ep_g(), get_mu_g());
    }

    {
      const bool constraint = !(m_idealgas_constraint == IdealGasConstraint::None);

      const bool update_lapT = (advect_enthalpy      && (l_explicit_diff || constraint));
      const bool update_lapS = (advect_tracer        &&  l_explicit_diff);
      const bool update_lapX = (advect_fluid_species && (l_explicit_diff || constraint));

      compute_laps(update_lapT, update_lapS, update_lapX, lap_T_old, lap_trac_old, lap_X_old,
                   get_T_g_old(), get_trac_old(), get_X_gk_old(), get_ep_g_const(),
                   get_ro_g_old_const());
    }


    // *************************************************************************************
    // Compute right hand side terms on the old status
    // *************************************************************************************

    if (advect_density) {
      mfix_density_rhs(ro_RHS_old, get_chem_txfr_const());
    }


    if (advect_enthalpy) {
      mfix_enthalpy_rhs(enthalpy_RHS_old, get_ep_g_const(), get_ro_g_old_const(),
          get_X_gk_old(), get_D_gk_const(), get_h_gk_const());
    }

    if (advect_tracer) {
      mfix_scalar_rhs(/*trac_RHS_old,*/ get_trac_old_const(), get_ep_g_const(),
          get_ro_g_old_const(), mu_s);
    }

    // Species
    if (advect_fluid_species) {
      mfix_species_X_rhs(species_RHS_old, get_chem_txfr_const());
    }


    // *************************************************************************************
    // Compute RHS for the MAC projection
    // *************************************************************************************

    if (m_idealgas_constraint == IdealGasConstraint::OpenSystem) {

      mfix_open_system_rhs(GetVecOfPtrs(rhs_mac), GetVecOfConstPtrs(lap_T_old),
          GetVecOfConstPtrs(enthalpy_RHS_old), GetVecOfConstPtrs(lap_X_old),
          GetVecOfConstPtrs(species_RHS_old), get_ro_g_old_const(),
          get_MW_g_const(), get_T_g_old_const(), get_cp_g_const(),
          get_X_gk_old(), get_D_gk_const(), get_h_gk_const(),
          get_txfr_const(), get_chem_txfr_const());

    } else if (m_idealgas_constraint == IdealGasConstraint::ClosedSystem) {

      mfix_closed_system_rhs(GetVecOfPtrs(rhs_mac), GetVecOfConstPtrs(lap_T_old),
          GetVecOfConstPtrs(enthalpy_RHS_old), GetVecOfConstPtrs(lap_X_old),
          GetVecOfConstPtrs(species_RHS_old), get_ep_g_const(),
          get_ro_g_old_const(), get_MW_g_const(), get_T_g_old_const(),
          get_cp_g_const(), get_X_gk_old(), get_D_gk_const(),
          get_h_gk_const(), get_txfr_const(), get_chem_txfr_const(),
          get_pressure_g_old_const(), avgSigma, avgTheta);
    }

    if (m_use_depdt_constraint) {
      for (int lev(0); lev <= finest_level; ++lev) {
        MultiFab::Subtract(rhs_mac[lev], depdt[lev], 0, 0, 1, 0);
      }
    }


    // *************************************************************************************
    // Compute the explicit advective terms
    // Note that "conv_u_old" returns update to (ep_g u)
    // Note that "conv_s_old" returns update to (ep_g rho), (ep_g rho h_g) and (ep_g rho tracer)
    // *************************************************************************************

    mfix_compute_convective_term(conv_u_old, conv_s_old, conv_X_old,
        GetVecOfPtrs(vel_forces), GetVecOfPtrs(tra_forces),
        get_vel_g_old_const(), get_ep_g_const(), get_ro_g_old_const(),
        get_h_g_old_const(), get_trac_old_const(), get_X_gk_old_const(), get_txfr_const(),
        GetVecOfPtrs(ep_u_mac), GetVecOfPtrs(ep_v_mac), GetVecOfPtrs(ep_w_mac),
        GetVecOfConstPtrs(rhs_mac), get_divtau(), l_dt, time);



    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    if (!advect_density)
    {
        for (int lev = 0; lev <= finest_level; lev++)
            MultiFab::Copy(density_nph[lev], *(m_leveldata[lev]->ro_go), 0, 0, 1, nghost_state());

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

        Real half_time = time + 0.5*l_dt;
        mfix_set_density_bcs(half_time, GetVecOfPtrs(density_nph));

    } // not constant density


    // **************************************************************************
    // Update thermodynamic pressure
    // **************************************************************************
    if (advect_enthalpy && (m_idealgas_constraint == IdealGasConstraint::ClosedSystem))
    {
      for (int lev = 0; lev <= finest_level; ++lev) {
        rhs_pressure_g_old[lev] = avgSigma[lev] / avgTheta[lev];

        auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*(ld.pressure_g),TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.tilebox();

          Array4<Real      > const& p_g     = ld.pressure_g->array(mfi);
          Array4<Real const> const& p_g_old = ld.pressure_go->const_array(mfi);
          const Real Dpressure_Dt           = rhs_pressure_g_old[lev];

          amrex::ParallelFor(bx, [p_g,p_g_old,Dpressure_Dt,l_dt]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            p_g(i,j,k) = p_g_old(i,j,k) + l_dt*Dpressure_Dt;
          });
        } // mfi
      }
    }


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

                const Real Dpressure_Dt           = rhs_pressure_g_old[lev];
                const int closed_system = (m_idealgas_constraint == IdealGasConstraint::ClosedSystem);

                amrex::ParallelFor(bx, [h_g_o,h_g_n,T_g_n,rho_o,rho_n,h_RHS_o,
                    epg,cp_g,dhdt_o,l_dt,lap_T_o,l_explicit_diff,Dpressure_Dt,
                    closed_system]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    int conv_comp = 1;

                    const Real num =         (rho_o(i,j,k) * epg(i,j,k));
                    const Real denom = 1.0 / (rho_n(i,j,k) * epg(i,j,k));

                    Real h_g = num * h_g_o(i,j,k) + l_dt * dhdt_o(i,j,k,conv_comp)
                                                  + l_dt * h_RHS_o(i,j,k);

                    if (closed_system)
                      h_g += l_dt * epg(i,j,k) * Dpressure_Dt;

                    if (l_explicit_diff)
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

                amrex::ParallelFor(bx, [tra_n,tra_o,rho_o,rho_n,lap_tra_o,epg,
                    dtdt_o,l_ntrac,l_dt,l_explicit_diff]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                  for (int n = 0; n < l_ntrac; ++n)
                  {
                    int conv_comp = 2+n;

                    const Real epg_loc = epg(i,j,k);

                    Real tra = (rho_o(i,j,k)*epg_loc)*tra_o(i,j,k,n);
                    tra += l_dt * dtdt_o(i,j,k,conv_comp);

                    if (l_explicit_diff)
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

          ParallelFor(bx, [nspecies_g,epg,rho_o,rho_n,X_gk_o,dXdt_o,lap_X_o,
              l_dt,X_gk_n,X_RHS_o,l_explicit_diff]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const Real num =         (rho_o(i,j,k) * epg(i,j,k));
            const Real denom = 1.0 / (rho_n(i,j,k) * epg(i,j,k));

            for (int n = 0; n < nspecies_g; ++n)
            {
              Real X_gk = num * X_gk_o(i,j,k,n) + l_dt * dXdt_o(i,j,k,n)
                                                + l_dt * X_RHS_o(i,j,k,n);

              if (l_explicit_diff)
                X_gk += l_dt * lap_X_o(i,j,k,n);

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
    compute_vel_forces(GetVecOfPtrs(vel_forces), get_vel_g_old_const(),
                       GetVecOfConstPtrs(density_nph), get_txfr_const());

    // *************************************************************************************
    // Update velocity with convective update, diffusive update, gp and gravity source terms
    // *************************************************************************************
    for (int lev = 0; lev <= finest_level; lev++)
    {

       auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(*ld.vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
         // Tilebox
         Box bx = mfi.tilebox ();

         Array4<Real      > const& vel_n    = ld.vel_g->array(mfi);
         Array4<Real const> const& vel_o    = ld.vel_go->const_array(mfi);
         Array4<Real const> const& divtau_o = ld.divtau_o->const_array(mfi);
         Array4<Real const> const& dudt_o   = conv_u_old[lev]->const_array(mfi);
         Array4<Real const> const& gp       = ld.gp->const_array(mfi);
         Array4<Real const> const& ro_g_o   = ld.ro_go->const_array(mfi);
         Array4<Real const> const& epg      = ld.ep_g->const_array(mfi);
         Array4<Real const> const& vel_f    = vel_forces[lev].const_array(mfi);

         // We need this until we remove static attribute from mfix::gravity
         const RealVect gp0_dev(gp0);
         const RealVect gravity_dev(gravity);

         amrex::ParallelFor(bx, [epg,ro_g_o,vel_o,dudt_o,divtau_o,gp0_dev,
             gravity_dev,gp,l_dt,vel_n,vel_f,l_explicit_diff]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
           const Real epg_loc = epg(i,j,k);
           const Real rog_loc = ro_g_o(i,j,k);
           const Real denom = 1.0 / (epg_loc*rog_loc);

           Real vel_nx = epg_loc*vel_o(i,j,k,0) + l_dt*dudt_o(i,j,k,0);
           Real vel_ny = epg_loc*vel_o(i,j,k,1) + l_dt*dudt_o(i,j,k,1);
           Real vel_nz = epg_loc*vel_o(i,j,k,2) + l_dt*dudt_o(i,j,k,2);

           vel_nx /= epg_loc;
           vel_ny /= epg_loc;
           vel_nz /= epg_loc;

           if (l_explicit_diff)
           {
             vel_nx += l_dt * (divtau_o(i,j,k,0) * denom);
             vel_ny += l_dt * (divtau_o(i,j,k,1) * denom);
             vel_nz += l_dt * (divtau_o(i,j,k,2) * denom);
           }

           vel_nx += l_dt * vel_f(i,j,k,0);
           vel_ny += l_dt * vel_f(i,j,k,1);
           vel_nz += l_dt * vel_f(i,j,k,2);

           vel_n(i,j,k,0) = vel_nx;
           vel_n(i,j,k,1) = vel_ny;
           vel_n(i,j,k,2) = vel_nz;
         });

       } // mfi
    } // lev

    // *************************************************************************************
    // Add the drag and convective heat transfer terms implicitly to vel_g and h_g
    // *************************************************************************************
    if (DEM::solve || PIC::solve)
      mfix_add_txfr_implicit(l_dt, get_vel_g(), get_h_g(), get_T_g(), get_txfr_const(),
                GetVecOfConstPtrs(density_nph), get_ep_g_const(), get_cp_g_const());


    // *************************************************************************************
    // If doing implicit diffusion, solve here for u^*
    // Note we multiply ep_g by ro_g so that we pass in a single array holding (ro_g * ep_g)
    // *************************************************************************************
    if (!l_explicit_diff) {

      mfix_set_density_bcs(time, get_ro_g());

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
        mfix_set_tracer_bcs(time, get_trac());
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
    Vector< MultiFab* > S_cc(finest_level+1);

    for (int lev(0); lev <= finest_level; ++lev) {
      S_cc[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0, 1).release();
    }

    if (!(m_idealgas_constraint == IdealGasConstraint::None)) {

      // Calculate drag coefficient
      if (DEM::solve || PIC::solve) {

        Real start_drag = ParallelDescriptor::second();
        amrex::Print() << "\nRecalculating drag ..." << std::endl;
        mfix_calc_txfr_fluid(get_txfr(), get_ep_g(), get_ro_g(), get_vel_g(),
                             get_mu_g(), get_cp_g(), get_k_g(), new_time);

        if (REACTIONS::solve) {
          mfix_calc_chem_txfr(get_chem_txfr(), get_ep_g(), get_ro_g(),
                              get_X_gk(), get_D_gk(), get_cp_gk(), get_h_gk(),
                              new_time);
        }

        coupling_timing += ParallelDescriptor::second() - start_drag;
      }

      const bool update_lapT = advect_enthalpy;
      const bool update_lapS = advect_tracer;
      const bool update_lapX = advect_fluid_species;

      compute_laps(update_lapT, update_lapS, update_lapX, lap_T, lap_trac, lap_X,
                   get_T_g(), get_trac(), get_X_gk(), get_ep_g_const(), get_ro_g_const());

      if (advect_density) {
        mfix_density_rhs(ro_RHS, get_chem_txfr_const());
      }

      if (advect_enthalpy) {
        mfix_enthalpy_rhs(enthalpy_RHS, get_ep_g_const(), get_ro_g_const(),
            get_X_gk(), get_D_gk_const(), get_h_gk_const());
      }

      if (advect_tracer) {
        mfix_scalar_rhs(/*trac_RHS,*/ get_trac_const(), get_ep_g_const(),
                        get_ro_g_const(), mu_s);
      }

      // Species
      if (advect_fluid_species) {
        mfix_species_X_rhs(species_RHS, get_chem_txfr_const());
      }
  
      if (m_idealgas_constraint == IdealGasConstraint::OpenSystem) {
  
        mfix_open_system_rhs(S_cc, GetVecOfConstPtrs(lap_T),
            GetVecOfConstPtrs(enthalpy_RHS), GetVecOfConstPtrs(lap_X),
            GetVecOfConstPtrs(species_RHS), get_ro_g_const(), get_MW_g_const(),
            get_T_g_const(), get_cp_g_const(), get_X_gk(),
            get_D_gk_const(), get_h_gk_const(), get_txfr_const(),
            get_chem_txfr_const());
  
      } else if (m_idealgas_constraint == IdealGasConstraint::ClosedSystem) {
  
        mfix_closed_system_rhs(S_cc, GetVecOfConstPtrs(lap_T),
            GetVecOfConstPtrs(enthalpy_RHS), GetVecOfConstPtrs(lap_X),
            GetVecOfConstPtrs(species_RHS), get_ep_g_const(), get_ro_g_const(),
            get_MW_g_const(), get_T_g_const(), get_cp_g_const(),
            get_X_gk(), get_D_gk_const(), get_h_gk_const(),
            get_txfr_const(), get_chem_txfr_const(), get_pressure_g_const(),
            avgSigma, avgTheta);
  
        // Update the thermodynamic pressure rhs in here so we do not have to call
        // the closed_system_rhs again in the corrector
        for (int lev = 0; lev <= finest_level; ++lev) {
          rhs_pressure_g[lev] = avgSigma[lev] / avgTheta[lev];
        }
      }
    }
  
    if (m_use_depdt_constraint) {
      for (int lev(0); lev <= finest_level; ++lev) {
        MultiFab::Subtract(*S_cc[lev], depdt[lev], 0, 0, 1, 0);
      }
    }
  
    mfix_apply_nodal_projection(S_cc, new_time, l_dt, l_prev_dt, proj_2,
                                get_vel_g_old(), get_vel_g(), get_p_g(), get_gp(),
                                get_ep_g(), get_txfr(), GetVecOfConstPtrs(density_nph));

    for (int lev(0); lev <= finest_level; ++lev) {
      delete S_cc[lev];
    }

    // *************************************************************************************
    // Correct small cells
    // *************************************************************************************
    mfix_correct_small_cells (get_vel_g(), GetVecOfConstPtrs(ep_u_mac),
         GetVecOfConstPtrs(ep_v_mac), GetVecOfConstPtrs(ep_w_mac));
}
