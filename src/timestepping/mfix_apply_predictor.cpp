#include <mfix.H>

#include <AMReX_VisMF.H>
#include <mfix_mf_helpers.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_fluid.H>
#include <mfix_species.H>
#include <mfix_solvers.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif


using namespace Solvers;


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
                            Vector< MultiFab* >& lap_trac_old,
                            Vector< MultiFab* >& enthalpy_RHS_old,
                            Vector< MultiFab* >& lap_T_old,
                            Vector< MultiFab* >& species_RHS_old,
                            Vector< MultiFab* >& vel_RHS_old,
                            Vector< MultiFab* >& div_J_old,
                            Vector< MultiFab* >& div_hJ_old,
                            Vector< Real >& rhs_pressure_g_old,
                            Vector< Real >& rhs_pressure_g,
                            Vector< MultiFab* > eb_flow_vel,
                            Vector< MultiFab* > eb_flow_scalars,
                            Vector< MultiFab* > eb_flow_species,
                            Real time,
                            Real l_dt,
                            Real l_prev_dt,
                            bool proj_2,
                            Real& /*coupling_timing*/)
{
    BL_PROFILE("mfix::mfix_apply_predictor");

    // We use the new-time value for things computed on the "*" state
    Real new_time = time + l_dt;

    Vector< MultiFab* > enthalpy_RHS(finest_level+1);
    Vector< MultiFab* > species_RHS(finest_level+1);

    for (int lev = 0; lev <= finest_level; lev++)
    {
       if (fluid.solve_enthalpy()) {
         enthalpy_RHS[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
         enthalpy_RHS[lev]->setVal(0.0);
       }

       if (fluid.solve_species()) {
         species_RHS[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies(), 0, MFInfo(), *ebfactory[lev]);
         species_RHS[lev]->setVal(0.0);
       }
    }

    // Averaged quantities for closed system constraint computations
    Vector<Real> avgSigma(finest_level+1, 0);
    Vector<Real> avgTheta(finest_level+1, 0);

    // Local flag for explicit diffusion
    bool l_explicit_diff = (predictor_diff_type() == DiffusionType::Explicit);

    m_boundary_conditions.set_density_bcs(time, get_ro_g_old());

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

      if (fluid.solve_tracer()) {
        tra_forces.emplace_back(grids[lev], dmap[lev], ntrac, nghost_force(),
                                MFInfo(), EBFactory(lev));
        tra_forces[lev].setVal(0.);
      }

      vel_eta.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(), EBFactory(lev));
      vel_eta[lev].setVal(0.);

      if (fluid.solve_tracer()) {
        tra_eta.emplace_back(grids[lev], dmap[lev], ntrac, 1, MFInfo(), EBFactory(lev));
        tra_eta[lev].setVal(0.);
      }
    }

    // *************************************************************************************
    // Allocate space for half-time density
    // *************************************************************************************
    Vector<MultiFab> density_nph;
    for (int lev = 0; lev <= finest_level; ++lev) {
        density_nph.emplace_back(grids[lev], dmap[lev], 1, nghost_state(), MFInfo(), *ebfactory[lev]);
        density_nph[lev].setVal(0.0);
    }


    // *************************************************************************************
    // Allocate space for the MAC velocities and RHS
    // *************************************************************************************
    Vector<MultiFab> ep_u_mac(finest_level+1), ep_v_mac(finest_level+1), ep_w_mac(finest_level+1);
    Vector<MultiFab> rhs_mac(finest_level+1),  depdt(finest_level+1);

    int ngmac = 1; //nghost_mac();

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
    Vector<Array<MultiFab*, AMREX_SPACEDIM>> J_gk(finest_level+1);

    if (fluid.solve_species()) {
      for (int lev = 0; lev <= finest_level; ++lev) {
        const int nspecies_g = fluid.nspecies();

        for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
          J_gk[lev][dir] = new MultiFab(amrex::convert(grids[lev], IntVect::TheDimensionVector(dir)),
                                        dmap[lev], nspecies_g, 1, MFInfo(), EBFactory(lev));
          J_gk[lev][dir]->setVal(0.);
        }
      }
    }

    if (need_divtau()) {
      diffusion_op->ComputeDivTau(get_divtau(), get_vel_g_old(), get_ep_g(), get_T_g_old(), GetVecOfConstPtrs(eb_flow_vel));
    }

    {
      const bool constraint = !(fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IncompressibleFluid);

      const bool update_lapT = (fluid.solve_enthalpy() && (l_explicit_diff || constraint));
      const bool update_lapS = (fluid.solve_tracer() &&  l_explicit_diff);
      const bool update_flux = (fluid.solve_species() && (l_explicit_diff || constraint));

      if (fluid.solve_enthalpy()) {
        m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g_old());
        m_boundary_conditions.set_enthalpy_bcs(time, fluid,get_h_g_old());
      }

      if (fluid.solve_species())
        m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk_old());

      compute_laps(update_lapT, update_lapS, update_flux, lap_T_old, lap_trac_old,
                   J_gk, get_T_g_old(), get_trac_old(), get_X_gk_old(), get_ep_g_const(),
                   get_ro_g_old_const());

      // We call the bc routines again to enforce the ext_dir condition
      // on the faces (the diffusion operator may move those to ghost cell centers)
      if (fluid.solve_enthalpy()) {
        m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g_old());
        m_boundary_conditions.set_enthalpy_bcs(time, fluid,get_h_g_old());
      }

      if (fluid.solve_species())
        m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk_old());
    }


    // *************************************************************************************
    //
    // *************************************************************************************
    if (fluid.solve_species())
      diffusion_op->ComputeDivJ(div_J_old, J_gk);


    // *************************************************************************************
    //
    // *************************************************************************************
    Vector< Array<MultiFab, AMREX_SPACEDIM> > h_gk_fc(finest_level+1);

    if (fluid.solve_enthalpy() && fluid.solve_species()) {
      for (int lev = 0; lev <= finest_level; ++lev) {
        for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
          h_gk_fc[lev][dir].define(amrex::convert(grids[lev], IntVect::TheDimensionVector(dir)),
                                   dmap[lev], fluid.nspecies(), 1, MFInfo(), EBFactory(lev));
          h_gk_fc[lev][dir].setVal(0.);
        }
      }

      const int update_enthalpies = 1;
      diffusion_op->ComputeDivhJ(div_hJ_old, h_gk_fc, J_gk, get_T_g_old_const(), update_enthalpies);
    }


    // *************************************************************************************
    // Compute right hand side terms on the old status
    // *************************************************************************************
    if (fluid.solve_density()) {
      mfix_density_rhs(ro_RHS_old, get_txfr_const());
    }

    if (fluid.solve_enthalpy()) {

      if (fluid.solve_species())
        m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk_old());

      m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g_old());
      m_boundary_conditions.set_enthalpy_bcs(time, fluid, get_h_g_old());

      mfix_enthalpy_rhs(enthalpy_RHS_old, get_ep_g_const(), get_ro_g_old_const(),
          get_X_gk_old(), get_T_g_old_const(), get_txfr_const());

      if (fluid.solve_species())
        m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk_old());

      m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g_old());
      m_boundary_conditions.set_enthalpy_bcs(time, fluid, get_h_g_old());
    }

    // Species
    if (fluid.solve_species()) {
      m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk_old());

      mfix_species_X_rhs(species_RHS_old, get_txfr_const());

      m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk_old());
    }

    // Linear momentum
    if (reactions.solve()) {
      mfix_momentum_rhs(vel_RHS_old, get_ep_g_const(), get_txfr_const());
    }


    // *************************************************************************************
    // Compute RHS for the MAC projection
    // *************************************************************************************
    if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IncompressibleFluid) {

      mfix_incompressible_fluid_rhs(GetVecOfPtrs(rhs_mac));

    } else {

      const int nspecies_g = fluid.nspecies();

      for (int lev(0); lev < nlev; ++lev) {
        if (fluid.solve_enthalpy()) {
          MultiFab::Add(*enthalpy_RHS_old[lev], *lap_T_old[lev], 0, 0, 1, 0);

          if (fluid.solve_species()) {
            MultiFab::Subtract(*enthalpy_RHS_old[lev], *div_hJ_old[lev], 0, 0, 1, 0);
          }
        }

        if (fluid.solve_species()) {
          MultiFab::Subtract(*species_RHS_old[lev], *div_J_old[lev], 0, 0, nspecies_g, 0);
        }
      }

      if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasOpenSystem) {

        mfix_idealgas_opensystem_rhs(GetVecOfPtrs(rhs_mac),
            GetVecOfConstPtrs(enthalpy_RHS_old), GetVecOfConstPtrs(species_RHS_old),
            get_ro_g_old_const(), get_T_g_old_const(), get_X_gk_old_const());

      } else if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem) {

        mfix_idealgas_closedsystem_rhs(GetVecOfPtrs(rhs_mac),
            GetVecOfConstPtrs(enthalpy_RHS_old), GetVecOfConstPtrs(species_RHS_old),
            get_ep_g_const(), get_ro_g_old_const(), get_T_g_old_const(),
            get_X_gk_old_const(), get_thermodynamic_p_g_old_const(), avgSigma, avgTheta);

      }

      for (int lev(0); lev < nlev; ++lev) {
        if (fluid.solve_enthalpy()) {
          MultiFab::Subtract(*enthalpy_RHS_old[lev], *lap_T_old[lev], 0, 0, 1, 0);

          if (fluid.solve_species()) {
            MultiFab::Add(*enthalpy_RHS_old[lev], *div_hJ_old[lev], 0, 0, 1, 0);
          }
        }

        if (fluid.solve_species()) {
          MultiFab::Add(*species_RHS_old[lev], *div_J_old[lev], 0, 0, nspecies_g, 0);
        }
      }
    }

    if (m_use_depdt_constraint) {
      for (int lev(0); lev <= finest_level; ++lev) {
        MultiFab::Subtract(rhs_mac[lev], depdt[lev], 0, 0, 1, 0);
      }
    }

    m_boundary_conditions.set_density_bcs(time, get_ro_g_old());

    if (fluid.solve_enthalpy()) {
      m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g_old());
      m_boundary_conditions.set_enthalpy_bcs(time, fluid, get_h_g_old());
    }

    if (fluid.solve_species())
      m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk_old());


    // *************************************************************************************
    // Compute the explicit advective terms
    // Note that "conv_u_old" returns update to (ep_g u)
    // Note that "conv_s_old" returns update to (ep_g rho), (ep_g rho h_g) and (ep_g rho tracer)
    // *************************************************************************************
    compute_MAC_projected_velocities(time, l_dt, get_vel_g_old_const(),
        GetVecOfPtrs(ep_u_mac), GetVecOfPtrs(ep_v_mac), GetVecOfPtrs(ep_w_mac),
        get_ep_g_const(), get_ro_g_old_const(), get_txfr_const(), GetVecOfConstPtrs(eb_flow_vel),
        GetVecOfPtrs(vel_forces), GetVecOfConstPtrs(rhs_mac));

    if (fluid.solve_species())
      m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk_old());

    mfix_compute_convective_term(conv_u_old, conv_s_old, conv_X_old,
        GetVecOfPtrs(vel_forces), GetVecOfPtrs(tra_forces), get_vel_g_old_const(),
        get_ep_g(), get_ro_g_old_const(), get_h_g_old_const(),
        get_trac_old_const(), get_X_gk_old_const(), get_txfr_const(),
        GetVecOfConstPtrs(eb_flow_vel), GetVecOfConstPtrs(eb_flow_scalars),
        GetVecOfConstPtrs(eb_flow_species), GetVecOfPtrs(ep_u_mac),
        GetVecOfPtrs(ep_v_mac), GetVecOfPtrs(ep_w_mac), l_dt, time);


    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    if (!fluid.solve_density())
    {
        for (int lev = 0; lev <= finest_level; lev++)
            MultiFab::Copy(density_nph[lev], *(m_leveldata[lev]->ro_go), 0, 0, 1, nghost_state());

    } else {

        const int nspecies_g = fluid.nspecies();
        const int use_species_advection = fluid.isMixture() && fluid.solve_species();

        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*ld.vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();

                Array4<Real const> dummy_arr;

                Array4<Real      > const& rho_new = ld.ro_g->array(mfi);
                Array4<Real      > const& rho_nph = density_nph[lev].array(mfi);
                Array4<Real const> const& rho_o   = ld.ro_go->const_array(mfi);
                Array4<Real const> const& epg     = ld.ep_g->const_array(mfi);
                Array4<Real const> const& drdt_o  = conv_s_old[lev]->const_array(mfi);
                Array4<Real const> const& dXdt_o  = use_species_advection ?
                                                    conv_X_old[lev]->const_array(mfi) : dummy_arr;
                Array4<Real const> const& rho_rhs = ro_RHS_old[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [rho_new,rho_nph,rho_o,epg,drdt_o,l_dt,
                    dXdt_o,use_species_advection, nspecies_g,rho_rhs]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                  const Real epg_loc = epg(i,j,k);

                  Real rho = epg_loc*rho_o(i,j,k);
                  if (use_species_advection) {
                    for (int n(0); n < nspecies_g; ++n)
                      rho += l_dt * dXdt_o(i,j,k,n);
                  } else {
                    constexpr int conv_comp = 0;
                    rho += l_dt * drdt_o(i,j,k,conv_comp);
                  }

                  rho += l_dt * rho_rhs(i,j,k);

                  rho /= epg_loc;

                  rho_new(i,j,k) = rho;
                  rho_nph(i,j,k) = 0.5 * (rho_o(i,j,k) + rho);
                });
            } // mfi
        } // lev

        Real half_time = time + 0.5*l_dt;
        m_boundary_conditions.set_density_bcs(half_time, GetVecOfPtrs(density_nph));

    } // not constant density


    // *************************************************************************
    // Update species mass fractions
    // *************************************************************************
    if (fluid.solve_species())
    {
      const int nspecies_g = fluid.nspecies();

      for (int lev = 0; lev <= finest_level; lev++)
      {
        auto& ld = *m_leveldata[lev];

        const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ld.X_gk->Factory());
        const auto& flags = factory.getMultiEBCellFlagFab();

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
          Array4<Real const> const& div_J_o = div_J_old[lev]->const_array(mfi);
          Array4<Real const> const& X_rhs_o = species_RHS_old[lev]->const_array(mfi);

          auto const& flags_arr = flags.const_array(mfi);

          ParallelFor(bx, [nspecies_g,epg,rho_o,rho_n,X_gk_o,dXdt_o,l_dt,X_gk_n,
              X_rhs_o,l_explicit_diff,div_J_o,flags_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const Real epg_loc = epg(i,j,k);

            const Real num =         (rho_o(i,j,k) * epg_loc);
            const Real denom = 1.0 / (rho_n(i,j,k) * epg_loc);

            Real X_gk_sum(0.);

            for (int n = 0; n < nspecies_g; ++n)
            {
              Real X_gk = num * X_gk_o(i,j,k,n);
              X_gk += l_dt * dXdt_o(i,j,k,n);
              X_gk += l_dt * X_rhs_o(i,j,k,n);

              if (l_explicit_diff) {
                X_gk -= l_dt * div_J_o(i,j,k,n);
              }

              X_gk *= denom;
              X_gk = amrex::max(amrex::min(X_gk, 1.), 0.);

              X_gk_sum += X_gk;

              X_gk_n(i,j,k,n) = X_gk;
            }

            if (!flags_arr(i,j,k).isCovered()) {
              for (int n = 0; n < nspecies_g; ++n)
              {
                X_gk_n(i,j,k,n) /= X_gk_sum;
              }
            }
          });
        } // mfi
      } // lev

      if (!l_explicit_diff) {
        // When using implicit diffusion for species, we "Add" (subtract) the
        // correction term computed at time t^{star,star} to the RHS before
        // doing the implicit diffusion
        m_boundary_conditions.set_epg_bcs(time, get_ep_g(), 0);
        m_boundary_conditions.set_density_bcs(time, get_ro_g());
        m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk());

        // Convert "ep_g" into (rho * ep_g)
        for (int lev = 0; lev <= finest_level; lev++)
          MultiFab::Multiply(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                             0, 0, 1, m_leveldata[lev]->ep_g->nGrow());

        // Diffuse species mass fractions
        diffusion_op->diffuse_species(get_X_gk(), J_gk, get_ep_g(), get_T_g_old(),
                                      get_species_bcrec(), l_dt);

        // Convert (rho * ep_g) back into ep_g
        for (int lev = 0; lev <= finest_level; lev++)
            MultiFab::Divide(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                             0, 0, 1, m_leveldata[lev]->ep_g->nGrow());

        m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk());

        // *********************************************************************
        // Correction
        // *********************************************************************
        diffusion_op->ComputeDivJ(div_J_old, J_gk);

        for (int lev = 0; lev <= finest_level; lev++) {

          auto& ld = *m_leveldata[lev];

          const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ld.X_gk->Factory());
          const auto& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(*ld.vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            Box const& bx = mfi.tilebox();

            Array4<Real      > const& X_gk_n  = ld.X_gk->array(mfi);
            Array4<Real const> const& X_gk_o  = ld.X_gko->const_array(mfi);
            Array4<Real const> const& rho_n   = ld.ro_g->const_array(mfi);
            Array4<Real const> const& rho_o   = ld.ro_go->const_array(mfi);
            Array4<Real const> const& epg     = ld.ep_g->const_array(mfi);
            Array4<Real const> const& dXdt_o  = conv_X_old[lev]->const_array(mfi);
            Array4<Real const> const& div_J_n = div_J_old[lev]->const_array(mfi);
            Array4<Real const> const& X_RHS_o = species_RHS_old[lev]->const_array(mfi);

            auto const& flags_arr = flags.const_array(mfi);

            ParallelFor(bx, [nspecies_g,epg,rho_o,rho_n,X_gk_o,dXdt_o,l_dt,X_gk_n,
                X_RHS_o,l_explicit_diff,div_J_n,flags_arr]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              const Real epg_loc = epg(i,j,k);

              const Real num =         (rho_o(i,j,k) * epg_loc);
              const Real denom = 1.0 / (rho_n(i,j,k) * epg_loc);

              Real X_gk_sum(0.);

              for (int n = 0; n < nspecies_g; ++n)
              {
                Real X_gk = num * X_gk_o(i,j,k,n);
                X_gk += l_dt * dXdt_o(i,j,k,n);
                X_gk += l_dt * X_RHS_o(i,j,k,n);

                X_gk -= l_dt * div_J_n(i,j,k,n);

                X_gk *= denom;
                X_gk = amrex::max(amrex::min(X_gk, 1.), 0.);

                X_gk_sum += X_gk;

                X_gk_n(i,j,k,n) = X_gk;
              }

              if (!flags_arr(i,j,k).isCovered()) {
                for (int n = 0; n < nspecies_g; ++n)
                {
                  X_gk_n(i,j,k,n) /= X_gk_sum;
                }
              }
            });
          }
        }

        // Note we need to call the bc routines again to enforce the ext_dir
        // condition on the faces (the diffusion operator moved those to ghost
        // cell centers)
        m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk());
      }

      // Update ghost cells
      for (int lev = 0; lev <= finest_level; lev++)
        m_leveldata[lev]->X_gk->FillBoundary(geom[lev].periodicity());

    } // solve_species


    // *************************************************************************
    // Update thermodynamic pressure
    // *************************************************************************
    if (fluid.solve_enthalpy() &&
        (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem))
    {
      for (int lev = 0; lev <= finest_level; ++lev) {
        rhs_pressure_g_old[lev] = avgSigma[lev] / avgTheta[lev];

        auto& ld = *m_leveldata[lev];

        Real* p_g     = ld.thermodynamic_p_g;
        Real* p_g_old = ld.thermodynamic_p_go;

        const Real Dpressure_Dt = rhs_pressure_g_old[lev];

        *p_g = *p_g_old + l_dt*Dpressure_Dt;
      }
    }


    // *************************************************************************************
    // Update enthalpy and temperature
    // *************************************************************************************
    if (fluid.solve_enthalpy()) {

      const auto& fluid_parms = fluid.parameters();
      const int fluid_is_a_mixture = fluid.isMixture();
      const int nspecies_g = fluid.nspecies();

      const int closed_system = (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem);

      if (!l_explicit_diff && fluid.solve_species()) {
        const int update_enthalpies = 0;
        diffusion_op->ComputeDivhJ(div_hJ_old, h_gk_fc, J_gk, get_T_g_old_const(), update_enthalpies);
      }

      for (int lev = 0; lev <= finest_level; lev++) {

        auto& ld = *m_leveldata[lev];

        const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ld.T_g->Factory());
        const auto& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*ld.vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

          Box const& bx = mfi.tilebox();

          Array4<Real const> dummy_arr;

          Array4<Real      > const& h_g_n    = ld.h_g->array(mfi);
          Array4<Real const> const& h_g_o    = ld.h_go->const_array(mfi);
          Array4<Real      > const& T_g_n    = ld.T_g->array(mfi);
          Array4<Real const> const& T_g_o    = ld.T_go->array(mfi);
          Array4<Real const> const& X_gk_n   = fluid_is_a_mixture ? ld.X_gk->array(mfi) : dummy_arr;
          Array4<Real const> const& rho_o    = ld.ro_go->const_array(mfi);
          Array4<Real const> const& rho_n    = ld.ro_g->const_array(mfi);
          Array4<Real const> const& epg      = ld.ep_g->const_array(mfi);
          Array4<Real const> const& dhdt_o   = conv_s_old[lev]->const_array(mfi);
          Array4<Real const> const& lap_T_o  = lap_T_old[lev]->const_array(mfi);
          Array4<Real const> const& div_hJ_o = fluid.solve_species()?
            div_hJ_old[lev]->const_array(mfi) : dummy_arr;

          const Real Dpressure_Dt = rhs_pressure_g_old[lev];

          auto const& flags_arr = flags.const_array(mfi);

          const int solve_species = fluid.solve_species();

          amrex::ParallelFor(bx, [h_g_o,h_g_n,T_g_o,T_g_n,rho_o,rho_n,epg,dhdt_o,
              l_dt,lap_T_o,l_explicit_diff,Dpressure_Dt,X_gk_n,closed_system,
              fluid_parms,fluid_is_a_mixture,nspecies_g,div_hJ_o,flags_arr,
              solve_species,abstol=newton_abstol,reltol=newton_reltol,
              maxiter=newton_maxiter]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

            if (!cell_is_covered) {
              int conv_comp = 1;

              const Real epg_loc = epg(i,j,k);

              const Real num =         (rho_o(i,j,k) * epg_loc);
              const Real denom = 1.0 / (rho_n(i,j,k) * epg_loc);

              Real h_g = num * h_g_o(i,j,k);

              h_g += l_dt * dhdt_o(i,j,k,conv_comp);

              if (solve_species)
                h_g -= l_dt * div_hJ_o(i,j,k);

              if (closed_system)
                h_g += l_dt * epg(i,j,k) * Dpressure_Dt;

              if (l_explicit_diff)
                h_g += l_dt * lap_T_o(i,j,k);

              h_g *= denom;

              h_g_n(i,j,k) = h_g;

              // Newton-Raphson solver for solving implicit equation for
              // temperature
              Newton::FluidEnthalpy::Residue residue(i, j, k, fluid_is_a_mixture,
                  cell_is_covered, nspecies_g, fluid_parms, X_gk_n, h_g);

              Newton::FluidEnthalpy::Gradient gradient(i, j, k, fluid_is_a_mixture,
                  nspecies_g, fluid_parms, X_gk_n);

              Real Tg(T_g_o(i,j,k));

              auto output = Newton::solve(Tg, residue, gradient, abstol, reltol, maxiter);

              if (output.iterations == -1)
                amrex::Abort("Newton solver did not converge");

              T_g_n(i,j,k) = Tg;
            }
          });
        } // mfi
      } // lev

      if (!l_explicit_diff) {

        m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g());
        m_boundary_conditions.set_enthalpy_bcs(time, fluid, get_h_g());

        // Diffuse temperature
        diffusion_op->diffuse_temperature(get_T_g(), get_ep_g(), get_ro_g(),
                                          get_h_g(), get_X_gk(), get_T_g_on_eb(), l_dt,
                                          newton_abstol, newton_reltol, newton_maxiter);

        // Note we need to call the bc routines again to enforce the ext_dir condition
        // on the faces (the diffusion operator moved those to ghost cell centers)
        m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g());
        m_boundary_conditions.set_enthalpy_bcs(time, fluid, get_h_g());
      }

      // ***********************************************************************
      // Add the dconvective heat transfer terms implicitly to h_g
      // ***********************************************************************
      if (m_dem.solve() || m_pic.solve())
        mfix_add_enthalpy_txfr_implicit(l_dt, get_h_g(), get_T_g(), get_X_gk_const(),
            get_txfr_const(), get_ro_g_const(), get_ep_g_const());

    } // fluid.solve_enthalpy


    // *************************************************************************************
    // Update tracer(s)
    // *************************************************************************************
    if (fluid.solve_tracer())
    {
        int l_ntrac = ntrac;

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

        if (!l_explicit_diff) {

          // Convert "ep_g" into (rho * ep_g)
          for (int lev = 0; lev <= finest_level; lev++)
            MultiFab::Multiply(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                               0, 0, 1, m_leveldata[lev]->ep_g->nGrow());

          // Diffuse tracer
          m_boundary_conditions.set_tracer_bcs(time, fluid, get_trac());

          diffusion_op->diffuse_scalar(get_trac(), get_ep_g(), mu_s, get_tracer_bcrec(), l_dt);

          // Note we need to call the bc routines again to enforce the ext_dir condition
          // on the faces (the diffusion operator moved those to ghost cell centers)
          m_boundary_conditions.set_tracer_bcs(time, fluid, get_trac());

          // Convert (rho * ep_g) back into ep_g
          for (int lev = 0; lev <= finest_level; lev++)
              MultiFab::Divide(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                               0, 0, 1, m_leveldata[lev]->ep_g->nGrow());
        }
    } // fluid.solve_tracer


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

         Array4<Real const> empty_array;

         Array4<Real      > const& vel_n     = ld.vel_g->array(mfi);
         Array4<Real const> const& vel_o     = ld.vel_go->const_array(mfi);
         Array4<Real const> const& rog_nph   = density_nph[lev].const_array(mfi);
         Array4<Real const> const& epg       = ld.ep_g->const_array(mfi);
         Array4<Real const> const& divtau_o  = ld.divtau_o->const_array(mfi);
         Array4<Real const> const& dudt_o    = conv_u_old[lev]->const_array(mfi);
         Array4<Real const> const& gp        = ld.gp->const_array(mfi);
         Array4<Real const> const& vel_f     = vel_forces[lev].const_array(mfi);
         Array4<Real const> const& vel_rhs_o = reactions.solve() ?
                                               vel_RHS_old[lev]->const_array(mfi) : empty_array;

         const int l_solve_reactions = reactions.solve();

         amrex::ParallelFor(bx, [epg,rog_nph,vel_o,dudt_o,divtau_o,gp,l_dt,vel_n,
             vel_f,l_explicit_diff,vel_rhs_o,l_solve_reactions]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
           const Real epg_loc = epg(i,j,k);
           const Real rog_loc = rog_nph(i,j,k);
           const Real denom = 1.0 / (epg_loc*rog_loc);

           Real vel_nx = epg_loc*vel_o(i,j,k,0) + l_dt*dudt_o(i,j,k,0);
           Real vel_ny = epg_loc*vel_o(i,j,k,1) + l_dt*dudt_o(i,j,k,1);
           Real vel_nz = epg_loc*vel_o(i,j,k,2) + l_dt*dudt_o(i,j,k,2);

           vel_nx /= epg_loc;
           vel_ny /= epg_loc;
           vel_nz /= epg_loc;

           if (l_explicit_diff) {
             vel_nx += l_dt * (divtau_o(i,j,k,0) * denom);
             vel_ny += l_dt * (divtau_o(i,j,k,1) * denom);
             vel_nz += l_dt * (divtau_o(i,j,k,2) * denom);
           }

           if (l_solve_reactions) {
             vel_nx += l_dt * (vel_rhs_o(i,j,k,0) * denom);
             vel_ny += l_dt * (vel_rhs_o(i,j,k,1) * denom);
             vel_nz += l_dt * (vel_rhs_o(i,j,k,2) * denom);
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
    // Add the drag transfer terms implicitly to vel_g
    // *************************************************************************************
    if (m_dem.solve() || m_pic.solve())
      mfix_add_vel_txfr_implicit(l_dt, get_vel_g(), get_txfr_const(),
          GetVecOfConstPtrs(density_nph), get_ep_g_const());


    // *************************************************************************************
    // If doing implicit diffusion, solve here for u^*
    // Note we multiply ep_g by ro_g so that we pass in a single array holding (ro_g * ep_g)
    // *************************************************************************************
    if (!l_explicit_diff) {

      m_boundary_conditions.set_density_bcs(time, get_ro_g());

      // Convert "ep_g" into (rho * ep_g)
      for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Multiply(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                           0, 0, 1, m_leveldata[lev]->ep_g->nGrow());

      // Set velocity boundary conditions
      m_boundary_conditions.set_velocity_bcs(new_time, get_vel_g(), 0);

      // Diffuse velocity
      diffusion_op->diffuse_velocity(get_vel_g(), get_ep_g(), get_T_g(), l_dt, GetVecOfConstPtrs(eb_flow_vel));

      // Convert (rho * ep_g) back into ep_g
      for (int lev = 0; lev <= finest_level; lev++)
          MultiFab::Divide(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                           0, 0, 1, m_leveldata[lev]->ep_g->nGrow());
    }


    // *************************************************************************************
    // Project velocity field -- depdt=0 for now
    // *************************************************************************************
    Vector< MultiFab* > S_cc(finest_level+1);

    for (int lev(0); lev <= finest_level; ++lev) {
      S_cc[lev] = new MultiFab(grids[lev], dmap[lev], 1, ngmac, MFInfo(), EBFactory(lev));
      S_cc[lev]->setVal(0.);
    }

    if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IncompressibleFluid) {

      mfix_incompressible_fluid_rhs(S_cc);

    } else {

      const int closed_system = int(fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem);

      const int nspecies_g = fluid.nspecies();

      for (int lev = 0; lev <= finest_level; lev++)
      {
          auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(*S_cc[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
              const int solve_species = fluid.solve_species();
              const int solve_enthalpy = fluid.solve_enthalpy();

              Box const& bx = mfi.tilebox();

              Array4<Real      > dummy_arr;
              Array4<Real const> const_dummy_arr;

              Array4<Real> const& species_RHS_arr  = solve_species ? species_RHS[lev]->array(mfi) : dummy_arr;
              Array4<Real> const& enthalpy_RHS_arr = solve_enthalpy ? enthalpy_RHS[lev]->array(mfi) : dummy_arr;

              Array4<Real const> const& X_gk_n  = solve_species ? ld.X_gk->const_array(mfi) : const_dummy_arr;
              Array4<Real const> const& X_gk_o  = solve_species ? ld.X_gko->const_array(mfi) : const_dummy_arr;
              Array4<Real const> const& h_g_n   = solve_enthalpy ? ld.h_g->const_array(mfi) : const_dummy_arr;
              Array4<Real const> const& h_g_o   = solve_enthalpy ? ld.h_go->const_array(mfi) : const_dummy_arr;
              Array4<Real const> const& rho_n   = ld.ro_g->const_array(mfi);
              Array4<Real const> const& rho_o   = ld.ro_go->const_array(mfi);
              Array4<Real const> const& epg     = ld.ep_g->const_array(mfi);

              Array4<Real const> const& dhdt_o  = solve_enthalpy ? conv_s_old[lev]->const_array(mfi) : const_dummy_arr;
              Array4<Real const> const& dXdt_o  = solve_species ? conv_X_old[lev]->const_array(mfi) : const_dummy_arr;

              const Real Dpressure_Dt = rhs_pressure_g_old[lev];

              amrex::ParallelFor(bx, [species_RHS_arr,enthalpy_RHS_arr,X_gk_n,
                  X_gk_o,h_g_n,h_g_o,rho_n,rho_o,epg,dhdt_o,dXdt_o,l_dt,nspecies_g,
                  closed_system,Dpressure_Dt,solve_species,solve_enthalpy]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                const Real epg_loc = epg(i,j,k);
                const Real ro_g_n  = rho_n(i,j,k);
                const Real ro_g_o  = rho_o(i,j,k);

                if (solve_enthalpy) {
                  enthalpy_RHS_arr(i,j,k) =
                    epg_loc*(ro_g_n*h_g_n(i,j,k) - ro_g_o*h_g_o(i,j,k)) / l_dt - dhdt_o(i,j,k,1);

                  if (closed_system)
                    enthalpy_RHS_arr(i,j,k) -= epg_loc*Dpressure_Dt;
                }

                if (solve_species) {
                  for (int n_g(0); n_g < nspecies_g; ++n_g) {
                    species_RHS_arr(i,j,k,n_g) =
                      epg_loc*(ro_g_n*X_gk_n(i,j,k,n_g) - ro_g_o*X_gk_o(i,j,k,n_g)) / l_dt - dXdt_o(i,j,k,n_g);
                  }
                }
              });
          } // mfi
      } // lev

      if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasOpenSystem) {

        mfix_idealgas_opensystem_rhs(S_cc, GetVecOfConstPtrs(enthalpy_RHS),
            GetVecOfConstPtrs(species_RHS), get_ro_g_const(),
            get_T_g_const(), get_X_gk_const());

      } else if (fluid.constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem) {

        mfix_idealgas_closedsystem_rhs(S_cc, GetVecOfConstPtrs(enthalpy_RHS),
            GetVecOfConstPtrs(species_RHS), get_ep_g_const(), get_ro_g_const(),
            get_T_g_const(), get_X_gk_const(), get_thermodynamic_p_g_const(),
            avgSigma, avgTheta);

        // Update the thermodynamic pressure rhs in here so we do not have to call
        // the closed_system_rhs again in the corrector
        for (int lev = 0; lev <= finest_level; ++lev) {
          rhs_pressure_g[lev] = fluid.solve_enthalpy()? avgSigma[lev] / avgTheta[lev] : 0.;
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
                                get_ep_g(), get_txfr(), GetVecOfConstPtrs(density_nph),
                                GetVecOfConstPtrs(eb_flow_vel));

    for (int lev(0); lev <= finest_level; ++lev) {
      delete S_cc[lev];
    }


    // *************************************************************************************
    // Correct small cells
    // *************************************************************************************
    mfix_correct_small_cells (get_vel_g(), GetVecOfConstPtrs(ep_u_mac),
         GetVecOfConstPtrs(ep_v_mac), GetVecOfConstPtrs(ep_w_mac),
         GetVecOfConstPtrs(eb_flow_vel));


    // *************************************************************************
    // Free up memory from temporary data
    // *************************************************************************
    for (int lev = 0; lev <= finest_level; lev++)
    {
      if (fluid.solve_enthalpy())
        delete enthalpy_RHS[lev];

      if (fluid.solve_species()) {
         delete species_RHS[lev];
      }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
      for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
        delete J_gk[lev][dir];
      }
    }
}
