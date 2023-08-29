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
    BL_PROFILE("mfix::mfix_apply_corrector");

    Vector< MultiFab* > lap_T(finest_level+1);
    Vector< MultiFab* > lap_trac(finest_level+1);
    Vector< MultiFab* > ro_RHS(finest_level+1);
    Vector< MultiFab* > enthalpy_RHS(finest_level+1);
    Vector< MultiFab* > species_RHS(finest_level+1);
    Vector< MultiFab* > vel_RHS(finest_level+1);

    for (int lev = 0; lev <= finest_level; lev++)
    {
       if (fluid.solve_enthalpy()) {
         lap_T[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
         lap_T[lev]->setVal(0.0);

         enthalpy_RHS[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
         enthalpy_RHS[lev]->setVal(0.0);
       }

       if (fluid.solve_tracer()) {
         lap_trac[lev] = new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), *ebfactory[lev]);
         lap_trac[lev]->setVal(0.0);
       }

       if (fluid.solve_density()) {
         ro_RHS[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
         ro_RHS[lev]->setVal(0.0);
       }

       if (fluid.solve_species()) {
         species_RHS[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies(), 0, MFInfo(), *ebfactory[lev]);
         species_RHS[lev]->setVal(0.0);
       }

       if (reactions.solve()) {
         vel_RHS[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
         vel_RHS[lev]->setVal(0.0);
       }
    }

    Vector< Real > avgSigma(finest_level+1, 0.);
    Vector< Real > avgTheta(finest_level+1, 0.);

    // We use the new-time value for things computed on the "*" state
    Real new_time = time + l_dt;

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

      if (fluid.solve_tracer()) {
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

    for (int lev = 0; lev <= finest_level; ++lev)
    {
      density_nph.emplace_back(grids[lev], dmap[lev], 1, nghost_state(), MFInfo(), *ebfactory[lev]);

      conv_u[lev] = new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
      // 2+ntrac components: one for density, ntrac for tracers and one for enthalpy
      conv_s[lev] = new MultiFab(grids[lev], dmap[lev], 2+ntrac, 0, MFInfo(), *ebfactory[lev]);

      density_nph[lev].setVal(0.0);

      conv_u[lev]->setVal(0.0);
      conv_s[lev]->setVal(0.0);

      if (fluid.solve_species()) {
        conv_X[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies(), 0, MFInfo(), *ebfactory[lev]);
        conv_X[lev]->setVal(0.0);
      }
    }


    // *************************************************************************************
    // Allocate space for the MAC velocities
    // *************************************************************************************
    Vector<MultiFab> ep_u_mac(finest_level+1), ep_v_mac(finest_level+1), ep_w_mac(finest_level+1);
    Vector<MultiFab> rhs_mac(finest_level+1);

    int ngmac = 1; //nghost_mac();

    for (int lev = 0; lev <= finest_level; ++lev) {
      ep_u_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(0)),
                           dmap[lev], 1, ngmac, MFInfo(), *ebfactory[lev]);
      ep_v_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(1)),
                           dmap[lev], 1, ngmac, MFInfo(), *ebfactory[lev]);
      ep_w_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(2)),
                           dmap[lev], 1, ngmac, MFInfo(), *ebfactory[lev]);

      if (ngmac > 0) {
        ep_u_mac[lev].setBndry(0.0);
        ep_v_mac[lev].setBndry(0.0);
        ep_w_mac[lev].setBndry(0.0);
      }

      rhs_mac[lev].define(grids[lev], dmap[lev], 1, ngmac, MFInfo(), EBFactory(lev));
      rhs_mac[lev].setVal(0.);
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

    const bool explicit_diffusive_enthalpy = false;
    const bool explicit_diffusive_trac     = false;
//    const bool explicit_diffusive_species  = false;
    const bool explicit_diffusive_species  = true;

    m_boundary_conditions.set_density_bcs(time, get_ro_g());

    {
      const bool constraint = !(fluid.constraint_type()== MFIXFluidPhase::ConstraintType::IncompressibleFluid);

      const bool update_lapT = (fluid.solve_enthalpy() && (explicit_diffusive_enthalpy || constraint));
      const bool update_lapS = (fluid.solve_tracer()   &&  explicit_diffusive_trac);
      const bool update_flux = (fluid.solve_species() && (explicit_diffusive_species || constraint));

      if (fluid.solve_enthalpy()) {
        m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g());
        m_boundary_conditions.set_enthalpy_bcs(time, fluid,get_h_g());
      }

      if (fluid.solve_species())
        m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk());

      compute_laps(update_lapT, update_lapS, update_flux, lap_T, lap_trac, J_gk,
                   get_T_g(), get_trac(), get_X_gk(), get_ep_g_const(),
                   get_ro_g_const());

      // We call the bc routines again to enforce the ext_dir condition
      // on the faces (the diffusion operator can move those to ghost cell centers)
      if (fluid.solve_enthalpy()) {
        m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g());
        m_boundary_conditions.set_enthalpy_bcs(time, fluid,get_h_g());
      }

      if (fluid.solve_species())
        m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk());
    }


    // *************************************************************************
    //
    // *************************************************************************
    Vector<MultiFab*> div_J(finest_level+1);

    if (fluid.solve_species()) {
      for (int lev = 0; lev <= finest_level; ++lev) {
        div_J[lev] = new MultiFab(grids[lev], dmap[lev], fluid.nspecies(), 0, MFInfo(), EBFactory(lev));
        div_J[lev]->setVal(0.);
      }

      diffusion_op->ComputeDivJ(div_J, J_gk);
    }


    // ************************************************************************
    //
    // *************************************************************************
    Vector<MultiFab*> div_hJ(finest_level+1);
    Vector< Array<MultiFab, AMREX_SPACEDIM> > h_gk_fc(finest_level+1);

    if (fluid.solve_enthalpy() && fluid.solve_species()) {
      for (int lev = 0; lev <= finest_level; ++lev) {
        div_hJ[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), EBFactory(lev));
        div_hJ[lev]->setVal(0.);

        for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
          h_gk_fc[lev][dir].define(amrex::convert(grids[lev], IntVect::TheDimensionVector(dir)),
                                   dmap[lev], fluid.nspecies(), 1, MFInfo(), EBFactory(lev));
          h_gk_fc[lev][dir].setVal(0.);
        }
      }

      const int update_enthalpies = 1;
      diffusion_op->ComputeDivhJ(div_hJ, h_gk_fc, J_gk, get_T_g_const(), update_enthalpies);
    }


    // ************************************************************************
    // Compute right hand side terms on the intermediate status
    // ************************************************************************
    if (fluid.solve_density()) {
      mfix_density_rhs(ro_RHS, get_txfr_const());
    }

    if (fluid.solve_enthalpy()) {
      mfix_enthalpy_rhs(enthalpy_RHS, get_ep_g_const(), get_ro_g_const(),
           get_X_gk(), get_T_g_const(), get_txfr_const());
    }

    if (fluid.solve_species()) {
      mfix_species_X_rhs(species_RHS, get_txfr_const());
    }

    // Linear momentum RHS
    if (reactions.solve()) {
      mfix_momentum_rhs(vel_RHS, get_ep_g_const(), get_txfr_const());
    }


    // *************************************************************************************
    // Compute RHS for the MAC projection
    // *************************************************************************************
    if (fluid.constraint_type()== MFIXFluidPhase::ConstraintType::IncompressibleFluid) {

      mfix_incompressible_fluid_rhs(GetVecOfPtrs(rhs_mac));

    } else {

      const int nspecies_g = fluid.nspecies();

      for (int lev(0); lev < nlev; ++lev) {
        if (fluid.solve_enthalpy()) {
          MultiFab::Add(*enthalpy_RHS[lev], *lap_T[lev], 0, 0, 1, 0);

          if (fluid.solve_species()) {
            MultiFab::Subtract(*enthalpy_RHS[lev], *div_hJ[lev], 0, 0, 1, 0);
          }
        }

        if (fluid.solve_species()) {
          MultiFab::Subtract(*species_RHS[lev], *div_J[lev], 0, 0, nspecies_g, 0);
        }
      }

      if (fluid.constraint_type()== MFIXFluidPhase::ConstraintType::IdealGasOpenSystem) {

        mfix_idealgas_opensystem_rhs(GetVecOfPtrs(rhs_mac),
            GetVecOfConstPtrs(enthalpy_RHS), GetVecOfConstPtrs(species_RHS),
            get_ro_g_const(), get_T_g_const(), get_X_gk_const());

      } else if (fluid.constraint_type()== MFIXFluidPhase::ConstraintType::IdealGasClosedSystem) {

        mfix_idealgas_closedsystem_rhs(GetVecOfPtrs(rhs_mac),
            GetVecOfConstPtrs(enthalpy_RHS), GetVecOfConstPtrs(species_RHS),
            get_ep_g_const(), get_ro_g_const(), get_T_g_const(), get_X_gk_const(),
            get_thermodynamic_p_g_const(), avgSigma, avgTheta);
      }

      for (int lev(0); lev < nlev; ++lev) {
        if (fluid.solve_enthalpy()) {
          MultiFab::Subtract(*enthalpy_RHS[lev], *lap_T[lev], 0, 0, 1, 0);

          if (fluid.solve_species()) {
            MultiFab::Add(*enthalpy_RHS[lev], *div_hJ[lev], 0, 0, 1, 0);
          }
        }

        if (fluid.solve_species()) {
          MultiFab::Add(*species_RHS[lev], *div_J[lev], 0, 0, nspecies_g, 0);
        }
      }
    }


    // *************************************************************************************
    // Compute the explicit advective term R_u^*
    // *************************************************************************************
    compute_MAC_projected_velocities(time, l_dt, get_vel_g_const(),
        GetVecOfPtrs(ep_u_mac), GetVecOfPtrs(ep_v_mac), GetVecOfPtrs(ep_w_mac),
        get_ep_g_const(), get_ro_g_const(), get_txfr_const(), GetVecOfConstPtrs(eb_flow_vel),
        GetVecOfPtrs(vel_forces), GetVecOfConstPtrs(rhs_mac));

    if (fluid.solve_species())
      m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk());

    mfix_compute_convective_term(conv_u, conv_s, conv_X, GetVecOfPtrs(vel_forces),
        GetVecOfPtrs(tra_forces), get_vel_g_const(), get_ep_g(), get_ro_g_const(),
        get_h_g_const(), get_trac_const(), get_X_gk_const(), get_txfr_const(),
        GetVecOfConstPtrs(eb_flow_vel), GetVecOfConstPtrs(eb_flow_scalars), GetVecOfConstPtrs(eb_flow_species),
        GetVecOfPtrs(ep_u_mac), GetVecOfPtrs(ep_v_mac), GetVecOfPtrs(ep_w_mac),
        l_dt, new_time);


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

                Array4<Real      > const& rho_new   = ld.ro_g->array(mfi);
                Array4<Real      > const& rho_nph   = density_nph[lev].array(mfi);
                Array4<Real const> const& rho_o     = ld.ro_go->const_array(mfi);
                Array4<Real const> const& epg       = ld.ep_g->array(mfi);
                Array4<Real const> const& drdt_o    = conv_s_old[lev]->const_array(mfi);
                Array4<Real const> const& drdt      = conv_s[lev]->const_array(mfi);
                Array4<Real const> const& dXdt_o    = use_species_advection ?
                                                      conv_X_old[lev]->const_array(mfi) : dummy_arr;
                Array4<Real const> const& dXdt      = use_species_advection ?
                                                      conv_X[lev]->const_array(mfi) : dummy_arr;
                Array4<Real const> const& rho_rhs_o = ro_RHS_old[lev]->const_array(mfi);
                Array4<Real const> const& rho_rhs   = ro_RHS[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [epg,rho_o,l_dt,drdt,drdt_o,rho_new,
                    use_species_advection,rho_nph,rho_rhs_o,rho_rhs,dXdt,dXdt_o,
                    nspecies_g]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                  const Real epg_loc = epg(i,j,k);
                  const Real rho_o_loc = rho_o(i,j,k);

                  Real rho = epg_loc*rho_o_loc;

                  if( use_species_advection) {
                    for (int n(0); n < nspecies_g; ++n)
                      rho += .5*l_dt*(dXdt_o(i,j,k,n)+dXdt(i,j,k,n));
                  } else {
                    constexpr int conv_comp = 0;
                    rho += .5*l_dt*(drdt_o(i,j,k,conv_comp)+drdt(i,j,k,conv_comp));
                  }

                  rho += .5*l_dt*(rho_rhs_o(i,j,k)+rho_rhs(i,j,k));

                  rho /= epg_loc;

                  rho_new(i,j,k) = rho;
                  rho_nph(i,j,k) = 0.5 * (rho_o_loc + rho);
                });
            } // mfi
        } // lev

        Real half_time = time + 0.5*l_dt;
        m_boundary_conditions.set_density_bcs(half_time, GetVecOfPtrs(density_nph));

    } // not constant density


    // *************************************************************************
    // Update species mass fraction
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
          Array4<Real const> const& dXdt    = conv_X[lev]->const_array(mfi);
          Array4<Real const> const& div_J_o = div_J_old[lev]->const_array(mfi);
          Array4<Real const> const& div_J_n = div_J[lev]->const_array(mfi);
          Array4<Real const> const& X_rhs_o = species_RHS_old[lev]->const_array(mfi);
          Array4<Real const> const& X_rhs   = species_RHS[lev]->const_array(mfi);

          auto const& flags_arr = flags.const_array(mfi);

          ParallelFor(bx, [X_gk_o,X_gk_n,rho_o,rho_n,epg,dXdt_o,dXdt,X_rhs_o,
              X_rhs,l_dt,nspecies_g,explicit_diffusive_species,div_J_o,div_J_n,
              flags_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const Real epg_loc = epg(i,j,k);

            const Real num =         (rho_o(i,j,k) * epg_loc);
            const Real denom = 1.0 / (rho_n(i,j,k) * epg_loc);

            Real X_gk_sum(0.);

            for (int n = 0; n < nspecies_g; ++n)
            {
              Real X_gk = num*X_gk_o(i,j,k,n);
              X_gk += .5*l_dt*(dXdt_o(i,j,k,n)+dXdt(i,j,k,n));
              X_gk += .5*l_dt*(X_rhs_o(i,j,k,n)+X_rhs(i,j,k,n));

              if (explicit_diffusive_species) {
                X_gk -= .5*l_dt*(div_J_o(i,j,k,n)+div_J_n(i,j,k,n));
              } else {
                // Crank-Nicolson so we only add half of the diffusive term
                // here, but we go ahead and add all of it now, then we will
                // subtract half of it
                X_gk -= .5*l_dt * div_J_o(i,j,k,n);
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

      if (!explicit_diffusive_species) {
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
        diffusion_op->diffuse_species(get_X_gk(), J_gk, get_ep_g(), get_T_g(),
                                      get_species_bcrec(), 0.5*l_dt);

        // Convert (rho * ep_g) back into ep_g
        for (int lev = 0; lev <= finest_level; lev++)
            MultiFab::Divide(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                             0, 0, 1, m_leveldata[lev]->ep_g->nGrow());

        m_boundary_conditions.set_species_bcs(time, fluid,get_X_gk());

        // *********************************************************************
        // Correction
        // *********************************************************************
        diffusion_op->ComputeDivJ(div_J, J_gk);

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
            Array4<Real const> const& dXdt    = conv_X[lev]->const_array(mfi);
            Array4<Real const> const& div_J_o = div_J_old[lev]->const_array(mfi);
            Array4<Real const> const& div_J_n = div_J[lev]->const_array(mfi);
            Array4<Real const> const& X_rhs_o = species_RHS_old[lev]->const_array(mfi);
            Array4<Real const> const& X_rhs   = species_RHS[lev]->const_array(mfi);

            auto const& flags_arr = flags.const_array(mfi);

            ParallelFor(bx, [X_gk_o,X_gk_n,rho_o,rho_n,epg,dXdt_o,dXdt,X_rhs_o,
                X_rhs,l_dt,nspecies_g,explicit_diffusive_species,div_J_o,div_J_n,
                flags_arr]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              const Real epg_loc = epg(i,j,k);

              const Real num =         (rho_o(i,j,k) * epg_loc);
              const Real denom = 1.0 / (rho_n(i,j,k) * epg_loc);

              Real X_gk_sum(0.);

              for (int n = 0; n < nspecies_g; ++n)
              {
                Real X_gk = num*X_gk_o(i,j,k,n);
                X_gk += .5*l_dt*(dXdt_o(i,j,k,n)+dXdt(i,j,k,n));
                X_gk += .5*l_dt*(X_rhs_o(i,j,k,n)+X_rhs(i,j,k,n));

                X_gk -= .5*l_dt*(div_J_o(i,j,k,n)+div_J_n(i,j,k,n));

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

    } // fluid.solve_species


    // **************************************************************************
    // Update thermodynamic pressure
    // **************************************************************************
    if (fluid.solve_enthalpy() &&
        (fluid.constraint_type()== MFIXFluidPhase::ConstraintType::IdealGasClosedSystem))
    {
      for (int lev = 0; lev <= finest_level; ++lev) {

        auto& ld = *m_leveldata[lev];

        Real* p_g     = ld.thermodynamic_p_g;
        Real* p_g_old = ld.thermodynamic_p_go;

        const Real Dpressure_Dt           = rhs_pressure_g[lev];
        const Real Dpressure_Dt_old       = rhs_pressure_g_old[lev];

        *p_g = *p_g_old + .5*l_dt*(Dpressure_Dt_old+Dpressure_Dt);
      }
    }


    // *************************************************************************************
    // Update enthalpy and temperature
    // *************************************************************************************
    if (fluid.solve_enthalpy()) {

      const auto& fluid_parms = fluid.parameters();
      const int fluid_is_a_mixture = fluid.isMixture();
      const int nspecies_g = fluid.nspecies();

      const int closed_system = (fluid.constraint_type()== MFIXFluidPhase::ConstraintType::IdealGasClosedSystem);

      if (!explicit_diffusive_species && fluid.solve_species()) {
        const int update_enthalpies = 0;
        diffusion_op->ComputeDivhJ(div_hJ, h_gk_fc, J_gk, get_T_g_const(), update_enthalpies);
      }

      for (int lev(0); lev < nlev; ++lev) {

        auto& ld = *m_leveldata[lev];

        const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ld.T_g->Factory());
        const auto& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*ld.vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

          const int solve_species = fluid.solve_species();

          Box const& bx = mfi.tilebox();

          Array4<Real const> dummy_arr;

          Array4<Real      > const& h_g_n    = ld.h_g->array(mfi);
          Array4<Real const> const& h_g_o    = ld.h_go->const_array(mfi);
          Array4<Real      > const& T_g_n    = ld.T_g->array(mfi);
          Array4<Real const> const& T_g_o    = ld.T_go->array(mfi);
          Array4<Real const> const& X_gk_n   = fluid_is_a_mixture ? ld.X_gk->array(mfi) : dummy_arr;
          Array4<Real const> const& rho_o    = ld.ro_go->const_array(mfi);
          Array4<Real const> const& rho_n    = ld.ro_g->const_array(mfi);
          Array4<Real const> const& epg      = ld.ep_g->array(mfi);
          Array4<Real const> const& dhdt_o   = conv_s_old[lev]->const_array(mfi);
          Array4<Real const> const& dhdt_n   = conv_s[lev]->const_array(mfi);
          Array4<Real const> const& lap_T_o  = lap_T_old[lev]->const_array(mfi);
          Array4<Real const> const& lap_T_n  = lap_T[lev]->const_array(mfi);
          Array4<Real const> const& h_rhs_o  = enthalpy_RHS_old[lev]->const_array(mfi);
          Array4<Real const> const& h_rhs_n  = enthalpy_RHS[lev]->const_array(mfi);
          Array4<Real const> const& div_hJ_o = solve_species ? div_hJ_old[lev]->const_array(mfi) : dummy_arr;
          Array4<Real const> const& div_hJ_n = solve_species ? div_hJ[lev]->const_array(mfi) : dummy_arr;

          const Real Dpressure_Dt           = rhs_pressure_g[lev];
          const Real Dpressure_Dt_old       = rhs_pressure_g_old[lev];

          auto const& flags_arr = flags.const_array(mfi);

          amrex::ParallelFor(bx, [h_g_o,h_g_n,T_g_o,T_g_n,rho_o,rho_n,epg,
              dhdt_o,dhdt_n,h_rhs_o,h_rhs_n,l_dt,lap_T_o,lap_T_n,Dpressure_Dt,
              Dpressure_Dt_old,closed_system,explicit_diffusive_enthalpy,
              fluid_parms,X_gk_n,nspecies_g,fluid_is_a_mixture,flags_arr,
              div_hJ_o,div_hJ_n,solve_species,abstol=newton_abstol,
              reltol=newton_reltol,maxiter=newton_maxiter]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

            if (!cell_is_covered) {
              int conv_comp = 1;
              const Real epg_loc = epg(i,j,k);

              const Real num =         (rho_o(i,j,k) * epg_loc);
              const Real denom = 1.0 / (rho_n(i,j,k) * epg_loc);

              Real h_g = num*h_g_o(i,j,k);
              h_g += .5*l_dt*(dhdt_o(i,j,k,conv_comp)+dhdt_n(i,j,k,conv_comp));
              h_g += .5*l_dt*(h_rhs_o(i,j,k)+h_rhs_n(i,j,k));

              if (solve_species)
                h_g -= .5*l_dt*(div_hJ_o(i,j,k)+div_hJ_n(i,j,k));

              if (explicit_diffusive_enthalpy) {
                h_g += .5*l_dt*(lap_T_o(i,j,k)+lap_T_n(i,j,k));
              }
              else {
                // Crank-Nicolson so we only add half of the diffusive term
                // here, but we go ahead and add all of it now before doing
                // the implicit enthalpy solve, then we will subtract half
                // of it after the enthalpy solve
                h_g += l_dt * lap_T_o(i,j,k);
              }

              if (closed_system) {
                h_g += .5*l_dt*epg_loc*(Dpressure_Dt+Dpressure_Dt_old);
              }

              h_g *= denom;

              h_g_n(i,j,k) = h_g;

              // ************************************************************
              // Newton-Raphson solver for solving implicit equation for
              // temperature
              // ************************************************************
              // Residual computation
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

      // *************************************************************************************
      // Subtract off half of the explicit diffusion terms (see comment above)
      // *************************************************************************************
      if (!explicit_diffusive_enthalpy) {

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

            Array4<Real      > const& h_g_n   = ld.h_g->array(mfi);
            Array4<Real      > const& T_g_n   = ld.T_g->array(mfi);
            Array4<Real const> const& T_g_o   = ld.T_go->array(mfi);
            Array4<Real const> const& ep_g    = ld.ep_g->const_array(mfi);
            Array4<Real const> const& ro_g_n  = ld.ro_g->const_array(mfi);
            Array4<Real const> const& X_gk_n  = fluid.solve_species()? ld.X_gk->array(mfi) : dummy_arr;
            Array4<Real const> const& lap_T_o = lap_T_old[lev]->const_array(mfi);

            auto const& flags_arr = flags.const_array(mfi);

            amrex::ParallelFor(bx, [ep_g,ro_g_n,h_g_n,T_g_o,T_g_n,lap_T_o,l_dt,
                fluid_parms,fluid_is_a_mixture,nspecies_g,X_gk_n,flags_arr,
                abstol=newton_abstol,reltol=newton_reltol,maxiter=newton_maxiter]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

              if (!cell_is_covered) {

                const Real epg_loc = ep_g(i,j,k);

                const Real denom = 1.0 / (ro_g_n(i,j,k) * epg_loc);

                Real h_g = h_g_n(i,j,k);
                h_g -= .5 * l_dt * (lap_T_o(i,j,k) * denom);

                h_g_n(i,j,k) = h_g;

                // ************************************************************
                // Newton-Raphson solver for solving implicit equation for
                // temperature
                // ************************************************************
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

        m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g());
        m_boundary_conditions.set_enthalpy_bcs(time, fluid, get_h_g());

        // NOTE: we do this call before multiplying ep_g by ro_g
        diffusion_op->diffuse_temperature(get_T_g(), get_ep_g(), get_ro_g(),
                                          get_h_g(), get_X_gk(), get_T_g_on_eb(), 0.5*l_dt,
                                          newton_abstol, newton_reltol, newton_maxiter);

        // We call the bc routines again to enforce the ext_dir condition
        // on the faces (the diffusion operator can move those to ghost cell centers)
        m_boundary_conditions.set_temperature_bcs(time, fluid, get_T_g());
        m_boundary_conditions.set_enthalpy_bcs(time, fluid, get_h_g());
      }

      // ***********************************************************************
      // Add the drag and enthalpy terms implicitly
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
                Array4<Real const> const& epg       = ld.ep_g->const_array(mfi);
                Array4<Real const> const& dtdt_o    = conv_s_old[lev]->const_array(mfi);
                Array4<Real const> const& dtdt_n    = conv_s[lev]->const_array(mfi);
                Array4<Real const> const& lap_tra_o = lap_trac_old[lev]->const_array(mfi);
                Array4<Real const> const& lap_tra_n = lap_trac[lev]->const_array(mfi);

                amrex::ParallelFor(bx, [tra_o,tra_n,rho_o,rho_n,epg,dtdt_o,dtdt_n,
                    lap_tra_o,lap_tra_n,l_dt,l_ntrac,explicit_diffusive_trac]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                  const Real epg_loc = epg(i,j,k);

                  const Real num =         (rho_o(i,j,k) * epg_loc);
                  const Real denom = 1.0 / (rho_n(i,j,k) * epg_loc);

                  for (int n = 0; n < l_ntrac; ++n)
                  {
                    int conv_comp = 2+n;

                    Real tra = num*tra_o(i,j,k,n);
                    tra += .5*l_dt*(dtdt_o(i,j,k,conv_comp)+dtdt_n(i,j,k,conv_comp));

                    if (explicit_diffusive_trac) {
                      tra += .5*l_dt*(lap_tra_o(i,j,k,n)+lap_tra_n(i,j,k,n));
                    }
                    else {
                      // Crank-Nicolson so we add the explicit half here
                      tra += .5*l_dt*lap_tra_o(i,j,k,n);
                    }

                    tra_n(i,j,k,n) = tra * denom;
                  }
                });
            } // mfi
        } // lev

        if (!explicit_diffusive_trac) {

            // Convert "ep_g" into (rho * ep_g)
            for (int lev = 0; lev <= finest_level; lev++)
                MultiFab::Multiply(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                                   0, 0, 1, m_leveldata[lev]->ep_g->nGrow());

            m_boundary_conditions.set_tracer_bcs(time, fluid, get_trac());

            // mfix_set_tracer_bcs (new_time, fluid, get_trac(), 0);
            diffusion_op->diffuse_scalar(get_trac(), get_ep_g(), mu_s, get_tracer_bcrec(), 0.5*l_dt);

            // We call the bc routines again to enforce the ext_dir condition
            // on the faces (the diffusion operator can move those to ghost cell centers)
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
    compute_vel_forces(GetVecOfPtrs(vel_forces), get_vel_g_const(),
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
       for (MFIter mfi(*m_leveldata[lev]->vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
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
         Array4<Real const> const& dudt_n    = conv_u[lev]->const_array(mfi);
         Array4<Real const> const& gp        = ld.gp->const_array(mfi);
         Array4<Real const> const& vel_f     = vel_forces[lev].const_array(mfi);
         Array4<Real const> const& vel_rhs_o = reactions.solve() ?
                                               vel_RHS_old[lev]->const_array(mfi) : empty_array;
         Array4<Real const> const& vel_rhs_n = reactions.solve() ?
                                               vel_RHS[lev]->const_array(mfi) : empty_array;

         const int l_solve_reactions = reactions.solve();

         amrex::ParallelFor(bx, [vel_n,vel_o,dudt_o,dudt_n,gp,vel_f,epg,rog_nph,
             divtau_o,l_dt,vel_rhs_o,vel_rhs_n,l_solve_reactions]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
           const Real epg_loc = epg(i,j,k);
           const Real rog_loc = rog_nph(i,j,k);
           const Real denom = 1.0 / (epg_loc*rog_loc);

           Real vel_nx = epg_loc*vel_o(i,j,k,0) + .5*l_dt*(dudt_o(i,j,k,0)+dudt_n(i,j,k,0));
           Real vel_ny = epg_loc*vel_o(i,j,k,1) + .5*l_dt*(dudt_o(i,j,k,1)+dudt_n(i,j,k,1));
           Real vel_nz = epg_loc*vel_o(i,j,k,2) + .5*l_dt*(dudt_o(i,j,k,2)+dudt_n(i,j,k,2));

           vel_nx /= epg_loc;
           vel_ny /= epg_loc;
           vel_nz /= epg_loc;

           // Crank-Nicolson so we should only add half of the explicit term here, but
           //     we go ahead and add all of it now before doing the implicit drag solve,
           //     then we will subtract half of it after the drag solve
           vel_nx += l_dt * (divtau_o(i,j,k,0) * denom);
           vel_ny += l_dt * (divtau_o(i,j,k,1) * denom);
           vel_nz += l_dt * (divtau_o(i,j,k,2) * denom);

           if (l_solve_reactions) {
             vel_nx += .5*l_dt * ((vel_rhs_o(i,j,k,0)+vel_rhs_n(i,j,k,0)) * denom);
             vel_ny += .5*l_dt * ((vel_rhs_o(i,j,k,1)+vel_rhs_n(i,j,k,1)) * denom);
             vel_nz += .5*l_dt * ((vel_rhs_o(i,j,k,2)+vel_rhs_n(i,j,k,2)) * denom);
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
    // Add the drag terms implicitly
    // *************************************************************************************
    if (m_dem.solve() || m_pic.solve())
      mfix_add_vel_txfr_implicit(l_dt, get_vel_g(), get_txfr_const(),
          GetVecOfConstPtrs(density_nph), get_ep_g_const());


    // *************************************************************************
    // Subtract half of the old divtau and then do implicit diffusion
    // *************************************************************************
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
         Array4<Real const> const& epg      = ld.ep_g->const_array(mfi);
         Array4<Real const> const& rog_nph  = density_nph[lev].const_array(mfi);
         Array4<Real const> const& divtau_o = ld.divtau_o->const_array(mfi);

         amrex::ParallelFor(bx, [vel_n,epg,rog_nph,divtau_o,l_dt]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
           const Real epg_loc = epg(i,j,k);
           const Real rog_loc = rog_nph(i,j,k);
           const Real denom = 1.0 / (epg_loc*rog_loc);

           vel_n(i,j,k,0) -= (.5 * l_dt) * (divtau_o(i,j,k,0) * denom);
           vel_n(i,j,k,1) -= (.5 * l_dt) * (divtau_o(i,j,k,1) * denom);
           vel_n(i,j,k,2) -= (.5 * l_dt) * (divtau_o(i,j,k,2) * denom);
         });
       }
    }

    // *************************************************************************************
    // Solve for u^star s.t. u^star = u_go + dt/2 (R_u^* + R_u^n) + dt/2 (Lu)^n + dt/2 (Lu)^star
    // Note we multiply ep_g by ro_g so that we pass in a single array holding (ro_g * ep_g)
    // *************************************************************************************

    // Convert "ep_g" into (rho * ep_g)
    for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Multiply(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                           0, 0, 1, m_leveldata[lev]->ep_g->nGrow());

    diffusion_op->diffuse_velocity(get_vel_g(), get_ep_g(), get_T_g(),
                                   0.5*l_dt, GetVecOfConstPtrs(eb_flow_vel));

    // Convert (rho * ep_g) back into ep_g
    for (int lev = 0; lev <= finest_level; lev++)
        MultiFab::Divide(*m_leveldata[lev]->ep_g, *m_leveldata[lev]->ro_g,
                         0, 0, 1, m_leveldata[lev]->ep_g->nGrow());


    // *************************************************************************************
    // Apply projection
    // *************************************************************************************
    Vector< MultiFab* > depdt(finest_level+1);
    Vector< MultiFab* > S_cc(finest_level+1);

    for (int lev(0); lev <= finest_level; ++lev) {
      depdt[lev] = new MultiFab(grids[lev], dmap[lev], 1, ngmac, MFInfo(), EBFactory(lev));
      S_cc[lev] = new MultiFab(grids[lev], dmap[lev], 1, ngmac, MFInfo(), EBFactory(lev));

      depdt[lev]->setVal(0.);
      S_cc[lev]->setVal(0.);
    }

    if (fluid.constraint_type()== MFIXFluidPhase::ConstraintType::IncompressibleFluid) {

      mfix_incompressible_fluid_rhs(S_cc);

    } else {

      const int closed_system = int(fluid.constraint_type()== MFIXFluidPhase::ConstraintType::IdealGasClosedSystem);

      const int nspecies_g = fluid.nspecies();

      for (int lev = 0; lev <= finest_level; lev++)
      {
          auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(*S_cc[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

              const int solve_species  = fluid.solve_species();
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

              Array4<Real const> const& dhdt_n  = solve_enthalpy ? conv_s[lev]->const_array(mfi) : const_dummy_arr;
              Array4<Real const> const& dXdt_n  = solve_species ? conv_X[lev]->const_array(mfi) : const_dummy_arr;
              Array4<Real const> const& dhdt_o  = solve_enthalpy ? conv_s_old[lev]->const_array(mfi) : const_dummy_arr;
              Array4<Real const> const& dXdt_o  = solve_species ? conv_X_old[lev]->const_array(mfi) : const_dummy_arr;

              const Real Dpressure_Dt           = rhs_pressure_g[lev];
              const Real Dpressure_Dt_old       = rhs_pressure_g_old[lev];

              amrex::ParallelFor(bx, [species_RHS_arr,enthalpy_RHS_arr,X_gk_n,
                  X_gk_o,h_g_n,h_g_o,rho_n,rho_o,epg,dhdt_o,dXdt_o,l_dt,
                  dhdt_n,dXdt_n,nspecies_g,closed_system,Dpressure_Dt_old,
                  Dpressure_Dt,solve_species,solve_enthalpy]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                const Real epg_loc = epg(i,j,k);
                const Real ro_g_n  = rho_n(i,j,k);
                const Real ro_g_o  = rho_o(i,j,k);

                if (solve_enthalpy) {
                  enthalpy_RHS_arr(i,j,k) =
                    epg_loc*(ro_g_n*h_g_n(i,j,k) - ro_g_o*h_g_o(i,j,k)) / l_dt - dhdt_n(i,j,k,1);

                  if (closed_system)
                    enthalpy_RHS_arr(i,j,k) -= epg_loc * Dpressure_Dt;
                }

                if (solve_species) {
                  for (int n_g(0); n_g < nspecies_g; ++n_g) {
                    species_RHS_arr(i,j,k,n_g) =
                      epg_loc*(ro_g_n*X_gk_n(i,j,k,n_g) - ro_g_o*X_gk_o(i,j,k,n_g)) / l_dt - dXdt_n(i,j,k,n_g);
                  }
                }
              });
          } // mfi
      } // lev

      if (fluid.constraint_type()== MFIXFluidPhase::ConstraintType::IdealGasOpenSystem) {

        mfix_idealgas_opensystem_rhs(S_cc, GetVecOfConstPtrs(enthalpy_RHS),
            GetVecOfConstPtrs(species_RHS), get_ro_g_const(),
            get_T_g_const(), get_X_gk_const());

      } else if (fluid.constraint_type()== MFIXFluidPhase::ConstraintType::IdealGasClosedSystem) {

        mfix_idealgas_closedsystem_rhs(S_cc, GetVecOfConstPtrs(enthalpy_RHS),
            GetVecOfConstPtrs(species_RHS), get_ep_g_const(), get_ro_g_const(),
            get_T_g_const(), get_X_gk_const(), get_thermodynamic_p_g_const(), avgSigma,
            avgTheta);
      }
    }

    if (m_use_depdt_constraint) {
      for (int lev(0); lev <= finest_level; ++lev) {
        MultiFab::Subtract(*S_cc[lev], *depdt[lev], 0, 0, 1, 0);
      }
    }

    mfix_apply_nodal_projection(S_cc, new_time, l_dt, l_prev_dt, proj_2,
                                get_vel_g_old(), get_vel_g(), get_p_g(), get_gp(),
                                get_ep_g(), get_txfr(), GetVecOfConstPtrs(density_nph),
                                GetVecOfConstPtrs(eb_flow_vel));

    for (int lev(0); lev <= finest_level; ++lev) {
      delete depdt[lev];
      delete S_cc[lev];
    }


    // *************************************************************************************
    // Correct small cells
    // *************************************************************************************
    mfix_correct_small_cells (get_vel_g(), GetVecOfConstPtrs(ep_u_mac),
         GetVecOfConstPtrs(ep_v_mac), GetVecOfConstPtrs(ep_w_mac),
         GetVecOfConstPtrs(eb_flow_vel));


    // *************************************************************************************
    // Free up memory from temporary data
    // *************************************************************************************
    for (int lev = 0; lev <= finest_level; lev++)
    {
       delete conv_u[lev];
       delete conv_s[lev];

       if (fluid.solve_species()) {
         delete conv_X[lev];
       }
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
       delete ro_RHS[lev];
       delete lap_trac[lev];
       delete enthalpy_RHS[lev];
       delete lap_T[lev];

       if (fluid.solve_species()) {
         delete species_RHS[lev];
       }

       if (reactions.solve())
         delete vel_RHS[lev];
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
      delete div_J[lev];
      delete div_hJ[lev];

      for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
        delete J_gk[lev][dir];
      }
    }
}
