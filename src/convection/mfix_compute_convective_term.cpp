#include <mfix.H>
#include <MOL.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_EB_utils.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>

#include <redistribution.H>
//
// Compute the three components of the convection term
//
void
mfix::mfix_compute_convective_term (const bool update_laplacians,
                                    Vector< MultiFab* >& conv_u_in,
                                    Vector< MultiFab* >& conv_s_in,
                                    Vector< MultiFab* >& conv_X_in,
                                    Vector< MultiFab* >& lap_T_out,
                                    Vector< MultiFab* >& lap_X_out,
                                    Vector< MultiFab* > const& vel_in,
                                    Vector< MultiFab* > const& ep_g_in,
                                    Vector< MultiFab* > const& ep_u_mac,
                                    Vector< MultiFab* > const& ep_v_mac,
                                    Vector< MultiFab* > const& ep_w_mac,
                                    Vector< MultiFab* > const& ro_g_in,
                                    Vector< MultiFab* > const& MW_g_in,
                                    Vector< MultiFab* > const& T_g_in,
                                    Vector< MultiFab* > const& cp_g_in,
                                    Vector< MultiFab* > const& k_g_in,
                                    Vector< MultiFab* > const& h_g_in,
                                    Vector< MultiFab* > const& T_g_on_eb_in,
                                    Vector< MultiFab* > const& k_g_on_eb_in,
                                    Vector< MultiFab* > const& trac_in,
                                    Vector< MultiFab* > const& X_gk_in,
                                    Vector< MultiFab* > const& D_gk_in,
                                    Vector< MultiFab* > const& h_gk_in,
                                    Vector< MultiFab* > const& txfr_in,
                                    Vector< MultiFab* > const& ro_gk_txfr_in,
                                    Real l_dt, Real time)
{
    BL_PROFILE("mfix::mfix_compute_convective_term");

    int ngmac = nghost_mac();

    const int l_nspecies = FLUID::nspecies;

    // First do FillPatch of {velocity, density, tracer, enthalpy} so we know
    // the ghost cells of these arrays are all filled
    for (int lev = 0; lev < nlev; lev++) {

      int state_comp, num_comp;

      // State with ghost cells
      MultiFab Sborder_u(grids[lev], dmap[lev], vel_in[lev]->nComp(), nghost,
                         MFInfo(), *ebfactory[lev]);
      FillPatchVel(lev, time, Sborder_u, 0, Sborder_u.nComp(), bcs_u);

      // Copy each FAB back from Sborder_u into the vel array, complete with filled ghost cells
      MultiFab::Copy(*vel_in[lev], Sborder_u, 0, 0, vel_in[lev]->nComp(), vel_in[lev]->nGrow());

      MultiFab Sborder_s(grids[lev], dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]);

      // We FillPatch density even if not advecting it because we need it in the projections
      state_comp =  0; // comp = 0 --> density
      num_comp = 1;
      FillPatchScalar(lev, time, Sborder_s, state_comp, num_comp, bcs_s);
      MultiFab::Copy(*ro_g_in[lev], Sborder_s, 0, 0, num_comp, ro_g_in[lev]->nGrow());

      if (advect_tracer) {
        state_comp =  1; // comp = 1 --> tracer
        num_comp = 1;
        FillPatchScalar(lev, time, Sborder_s, state_comp, num_comp, bcs_s);
        MultiFab::Copy(*trac_in[lev], Sborder_s, 0, 0, num_comp, trac_in[lev]->nGrow());

      }

      if (advect_enthalpy) {
        state_comp =  5; // comp = 1 --> enthalpy
        num_comp = 1;
        FillPatchScalar(lev, time, Sborder_s, state_comp, num_comp, bcs_s);
        MultiFab::Copy(*h_g_in[lev], Sborder_s, 0, 0, num_comp, h_g_in[lev]->nGrow());

      }

      if (advect_fluid_species) {
        MultiFab Sborder_X(grids[lev], dmap[lev], FLUID::nspecies, nghost,
                           MFInfo(), *ebfactory[lev]);

        Sborder_X.setVal(0);

        state_comp = 0;
        num_comp = FLUID::nspecies;

        FillPatchSpecies(lev, time, Sborder_X, state_comp, num_comp, bcs_X);

        MultiFab::Copy(*X_gk_in[lev], Sborder_X, 0, 0, num_comp,
                       X_gk_in[lev]->nGrow());

      }
    }



    // We first compute the velocity forcing terms to be used in predicting
    //    to faces before the MAC projection
    if (m_advection_type != "MOL") {
#if 0
      bool include_pressure_gradient = !(m_use_mac_phi_in_godunov);
      compute_vel_forces(vel_forces, vel, density, tracer, tracer, include_pressure_gradient);

      if (m_godunov_include_diff_in_forcing)
        for (int lev = 0; lev <= finest_level; ++lev)
          MultiFab::Add(*vel_forces[lev], m_leveldata[lev]->divtau_o, 0, 0, AMREX_SPACEDIM, 0);

      if (nghost_force() > 0)
        fillpatch_force(m_cur_time, vel_forces, nghost_force());
#endif
    }





    // Do projection on all AMR levels in one shot -- note that the {u_mac, v_mac, w_mac}
    //    arrays returned from this call are in fact {ep * u_mac, ep * v_mac, ep * w_mac}
    //    on face CENTROIDS
    compute_MAC_projected_velocities(time, vel_in, ep_u_mac, ep_v_mac, ep_w_mac,
       ep_g_in, ro_g_in, MW_g_in, T_g_in, cp_g_in, k_g_in, T_g_on_eb_in,
       k_g_on_eb_in, X_gk_in, D_gk_in, h_gk_in, txfr_in, ro_gk_txfr_in,
       update_laplacians, lap_T_out, lap_X_out);





    // We now re-compute the velocity forcing terms including the pressure gradient,
    //    and compute the tracer forcing terms for the first time
    if (m_advection_type != "MOL") {
#if 0
      compute_vel_forces(vel_forces, vel, density, tracer, tracer);

      if (m_godunov_include_diff_in_forcing)
        for (int lev = 0; lev <= finest_level; ++lev)
          MultiFab::Add(*vel_forces[lev], m_leveldata[lev]->divtau_o, 0, 0, AMREX_SPACEDIM, 0);

      if (nghost_force() > 0)
        fillpatch_force(m_cur_time, vel_forces, nghost_force());

      // Note this is forcing for (rho s), not for s
      if (advect_tracer) {
        compute_tra_forces(tra_forces, get_density_old_const());
        if (m_godunov_include_diff_in_forcing)
          for (int lev = 0; lev <= finest_level; ++lev)
            MultiFab::Add(*tra_forces[lev], m_leveldata[lev]->laps_o, 0, 0, ntrac, 0);
        if (nghost_force() > 0)
          fillpatch_force(m_cur_time, tra_forces, nghost_force());
      }
#endif
    }


    for (int lev=0; lev < nlev; ++lev) {
      if (ngmac > 0) {
        ep_u_mac[lev]->FillBoundary(geom[lev].periodicity());
        ep_v_mac[lev]->FillBoundary(geom[lev].periodicity());
        ep_w_mac[lev]->FillBoundary(geom[lev].periodicity());
      }

      MultiFab divu(vel_in[lev]->boxArray(),vel_in[lev]->DistributionMap(),1,2);
      divu.setVal(0.);
      Array<MultiFab const*, AMREX_SPACEDIM> u;
      u[0] = ep_u_mac[lev];
      u[1] = ep_v_mac[lev];
      u[2] = ep_w_mac[lev];

      auto const& fact = EBFactory(lev);
      if (fact.isAllRegular())
        computeDivergence(divu,u,geom[lev]);
      else
        EB_computeDivergence(divu,u,geom[lev],true);

      divu.FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      // Turn off tiling -- HACK HACK HACK
      for (MFIter mfi(*ep_g_in[lev],false); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();
        mfix::compute_convective_term(bx, lev, mfi,
                                      conv_u_in[lev]->array(mfi),
                                      conv_s_in[lev]->array(mfi),
                                      (l_nspecies>0) ? conv_X_in[lev]->array(mfi) : Array4<Real>{},
                                      advect_density, advect_enthalpy, advect_tracer, advect_fluid_species,
                                      vel_in[lev]->const_array(mfi),
                                      ep_g_in[lev]->const_array(mfi),
                                      ro_g_in[lev]->const_array(mfi),
                                      (advect_enthalpy) ?  h_g_in[lev]->const_array(mfi) : Array4<Real const>{},
                                      (advect_tracer) ? trac_in[lev]->const_array(mfi) : Array4<Real const>{},
                                      (l_nspecies > 0)  ? X_gk_in[lev]->const_array(mfi) : Array4<Real const>{},
                                      l_nspecies,
                                      divu.const_array(mfi),
                                      ep_u_mac[lev]->const_array(mfi),
                                      ep_v_mac[lev]->const_array(mfi),
                                      ep_w_mac[lev]->const_array(mfi),
                                      l_dt);
      } // mfi

    } // lev

}




















void
mfix::compute_convective_term (Box const& bx, int lev, MFIter const& mfi,
                               Array4<Real> const& dvdt, // velocity
                               Array4<Real> const& dsdt, // density, enthalpy, tracer
                               Array4<Real> const& dXdt, // tracer
                               const bool l_advect_density,
                               const bool l_advect_enthalpy,
                               const bool l_advect_tracer,
                               const bool l_advect_species,
                               Array4<Real const> const& vel,
                               Array4<Real const> const& ep_g,
                               Array4<Real const> const& rho,
                               Array4<Real const> const& hg,
                               Array4<Real const> const& tra,
                               Array4<Real const> const& Xgk,
                               const int l_nspecies,
                               Array4<Real const> const& divu,
                               Array4<Real const> const& ep_umac,
                               Array4<Real const> const& ep_vmac,
                               Array4<Real const> const& ep_wmac,
                               const Real l_dt)
{

  auto const& fact = EBFactory(lev);
  EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
  Array4<EBCellFlag const> const& flag = flagfab.const_array();
  if (flagfab.getType(bx) == FabType::covered) {
    amrex::ParallelFor(bx, [=]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      dvdt(i,j,k,0) = 0.0;
      dvdt(i,j,k,1) = 0.0;
      dvdt(i,j,k,2) = 0.0;

      dsdt(i,j,k,0) = 0.0; // density
      dsdt(i,j,k,1) = 0.0; // enthalpy
      dsdt(i,j,k,2) = 0.0; // tracer

    });

    if (l_advect_species) {
      amrex::ParallelFor(bx, l_nspecies, [=]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        dXdt(i,j,k,n) = 0.0;
      });
    }

    return;
  }

    bool regular = (flagfab.getType(amrex::grow(bx,2)) == FabType::regular);

    Array4<Real const> fcx, fcy, fcz, ccc, vfrac, apx, apy, apz;
    if (!regular) {
        fcx = fact.getFaceCent()[0]->const_array(mfi);
        fcy = fact.getFaceCent()[1]->const_array(mfi);
        fcz = fact.getFaceCent()[2]->const_array(mfi);
        ccc   = fact.getCentroid().const_array(mfi);
        apx = fact.getAreaFrac()[0]->const_array(mfi);
        apy = fact.getAreaFrac()[1]->const_array(mfi);
        apz = fact.getAreaFrac()[2]->const_array(mfi);
        vfrac = fact.getVolFrac().const_array(mfi);
    }

    Box rho_box = amrex::grow(bx,2);
    if (m_advection_type != "MOL")  rho_box.grow(1);
    if (!regular) rho_box.grow(2);

    FArrayBox rhohgfab, rhotracfab, rhoXgkfab;
    Elixir eli_rh, eli_rt, eli_rXgk;
    Array4<Real> rhohg, rhotrac, rhoXgk;

    if (l_advect_enthalpy) {
      rhohgfab.resize(rho_box, 1);
      eli_rh = rhohgfab.elixir();
      rhohg = rhohgfab.array();
      amrex::ParallelFor(rho_box, [rho, hg, rhohg]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        rhohg(i,j,k) = rho(i,j,k) * hg(i,j,k);
      });
    }

    if (l_advect_tracer) {
      rhotracfab.resize(rho_box, ntrac);
      eli_rt = rhotracfab.elixir();
      rhotrac = rhotracfab.array();
      amrex::ParallelFor(rho_box, ntrac, [rho, tra, rhotrac]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        rhotrac(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
      });
    }

    if (l_advect_species) {
      rhoXgkfab.resize(rho_box, l_nspecies);
      eli_rXgk = rhoXgkfab.elixir();
      rhoXgk = rhoXgkfab.array();
      amrex::ParallelFor(rho_box, l_nspecies, [rho, Xgk, rhoXgk]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        rhoXgk(i,j,k,n) = rho(i,j,k) * Xgk(i,j,k,n);
      });
    }

    // Figure out how big the temp array needs to be to hold
    // all of the variables we plan to advect. Since we reuse
    // the temp array, it only needs to be as big as the
    // largest scalar array.
    int nmaxcomp = AMREX_SPACEDIM;
    if (l_advect_tracer)  nmaxcomp = std::max(nmaxcomp,ntrac);
    if (l_advect_species) nmaxcomp = std::max(nmaxcomp,l_nspecies);


    const GpuArray<int, 3> bc_types =
      {bc_list.get_minf(), bc_list.get_pout(), bc_list.get_pinf()};

    if (m_advection_type == "Godunov") {
      amrex::Abort("Godunov not implemented!");

    } else if (m_advection_type == "MOL") {
      Box tmpbox = amrex::surroundingNodes(bx);
      int tmpcomp = nmaxcomp*3; // fx, fy, fz
      Box gbx = bx;
      if (!regular) {
        gbx.grow(2);
        tmpbox.grow(3);
        tmpcomp += nmaxcomp; // dUdt_tmp (eb only)
      }

      FArrayBox tmpfab(tmpbox, tmpcomp);
      Elixir eli = tmpfab.elixir();

      Array4<Real> fx = tmpfab.array(0);
      Array4<Real> fy = tmpfab.array(nmaxcomp);
      Array4<Real> fz = tmpfab.array(nmaxcomp*2);

      if (!regular) {

        Array4<Real> dUdt_tmp = tmpfab.array(nmaxcomp*3);

        // ***************************************************************************
        // Compute div (ep_u_mac * (u)) -- the update for velocity
        // ***************************************************************************
        {// Velocity
          const int scomp = 0; // starting component
          const int ncomp = 3; // number of components
          const int ccomp = 0; // convection (dsdt) component

          mol::compute_convective_fluxes_eb(gbx, ncomp, scomp, vel,
                                            fx, fy, fz,
                                            ep_umac, ep_vmac, ep_wmac,
                                            flag, fcx, fcy, fcz, ccc,
                                            bc_types, m_vel_g_bc_types,
                                            bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                            bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                            bc_klo[lev]->array(), bc_khi[lev]->array(),
                                            geom[lev]);

          mol::compute_convective_rate_eb(gbx, ncomp, scomp, dUdt_tmp, fx, fy, fz,
                                          flag, vfrac, apx, apy, apz, geom[lev]);

          redistribution::apply_eb_redistribution(bx, dvdt, dUdt_tmp, ep_g, mfi,
                                                  ccomp, ncomp, flag, vfrac, geom[lev]);

        } // end velocity

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho)) -- the update for density
        // ***************************************************************************
        if (l_advect_density) {

          const int scomp = 0; // starting component
          const int ncomp = 1; // number of components
          const int ccomp = 0; // convection (dsdt) component

          mol::compute_convective_fluxes_eb(gbx, ncomp, scomp, rho,
                                            fx, fy, fz,
                                            ep_umac, ep_vmac, ep_wmac,
                                            flag, fcx, fcy, fcz, ccc,
                                            bc_types, m_ro_g_bc_types,
                                            bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                            bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                            bc_klo[lev]->array(), bc_khi[lev]->array(),
                                            geom[lev]);

          mol::compute_convective_rate_eb(gbx, ncomp, scomp, dUdt_tmp, fx, fy, fz,
                                          flag, vfrac, apx, apy, apz, geom[lev]);

          redistribution::apply_eb_redistribution(bx, dsdt, dUdt_tmp, ep_g, mfi,
                                                  ccomp, ncomp, flag, vfrac, geom[lev]);

        } // end density

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho*h_g)) -- the update for enthalpy
        // ***************************************************************************
        if (l_advect_enthalpy) {

          const int scomp = 0; // starting component
          const int ncomp = 1; // number of components
          const int ccomp = 1; // convection (dsdt) component

          mol::compute_convective_fluxes_eb(gbx, ncomp, scomp, rhohg,
                                            fx, fy, fz,
                                            ep_umac, ep_vmac, ep_wmac,
                                            flag, fcx, fcy, fcz, ccc,
                                            bc_types, m_T_g_bc_types,
                                            bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                            bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                            bc_klo[lev]->array(), bc_khi[lev]->array(),
                                            geom[lev]);

          mol::compute_convective_rate_eb(gbx, ncomp, scomp, dUdt_tmp, fx, fy, fz,
                                          flag, vfrac, apx, apy, apz, geom[lev]);

          redistribution::apply_eb_redistribution(bx, dsdt, dUdt_tmp, ep_g, mfi,
                                                  ccomp, ncomp, flag, vfrac, geom[lev]);
        } // end enthalpy

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho*trac)) -- the update for tracer
        // ***************************************************************************
        if (l_advect_tracer) {

          const int scomp = 0; // starting component
          const int ncomp = 1; // number of components
          const int ccomp = 2; // convection (dsdt) component

          mol::compute_convective_fluxes_eb(gbx, ncomp, scomp, rhotrac,
                                            fx, fy, fz,
                                            ep_umac, ep_vmac, ep_wmac,
                                            flag, fcx, fcy, fcz, ccc,
                                            bc_types, m_trac_g_bc_types,
                                            bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                            bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                            bc_klo[lev]->array(), bc_khi[lev]->array(),
                                            geom[lev]);

          mol::compute_convective_rate_eb(gbx, ncomp, scomp, dUdt_tmp, fx, fy, fz,
                                          flag, vfrac, apx, apy, apz, geom[lev]);

          redistribution::apply_eb_redistribution(bx, dsdt, dUdt_tmp, ep_g, mfi,
                                                  ccomp, ncomp, flag, vfrac, geom[lev]);
        } // end tracer

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho*X_g)) -- the update for species mass fraction
        // ***************************************************************************
        if (l_advect_species) {

          const int scomp = 0; // starting component
          const int ncomp = l_nspecies; // number of components
          const int ccomp = 0; // convection (dXdt) component

          mol::compute_convective_fluxes_eb(gbx, ncomp, scomp, rhoXgk,
                                            fx, fy, fz,
                                            ep_umac, ep_vmac, ep_wmac,
                                            flag, fcx, fcy, fcz, ccc,
                                            bc_types, m_X_gk_bc_types,
                                            bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                            bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                            bc_klo[lev]->array(), bc_khi[lev]->array(),
                                            geom[lev]);

          mol::compute_convective_rate_eb(gbx, ncomp, scomp, dUdt_tmp, fx, fy, fz,
                                          flag, vfrac, apx, apy, apz, geom[lev]);

          redistribution::apply_eb_redistribution(bx, dXdt, dUdt_tmp, ep_g, mfi,
                                                  ccomp, ncomp, flag, vfrac, geom[lev]);

        } // end species


      } else { // regular fab

        {// Velocity
          const int scomp = 0; // starting component
          const int ncomp = 3; // number of components
          const int ccomp = 0; // convection (dvdt) component

          mol::compute_convective_fluxes(bx, ncomp, scomp, vel,
                                         fx, fy, fz,
                                         ep_umac, ep_vmac, ep_wmac,
                                         bc_types, m_vel_g_bc_types,
                                         bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                         bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                         bc_klo[lev]->array(), bc_khi[lev]->array(),
                                         geom[lev]);

          mol::compute_convective_rate(bx, ncomp, ccomp, dvdt, fx, fy, fz, geom[lev]);

        } // end velocity

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho)) -- the update for density
        // ***************************************************************************
        if (l_advect_density) {

          const int scomp = 0; // starting component
          const int ncomp = 1; // number of components
          const int ccomp = 0; // convection (dsdt) component

          mol::compute_convective_fluxes(bx, ncomp, scomp, rho,
                                         fx, fy, fz,
                                         ep_umac, ep_vmac, ep_wmac,
                                         bc_types, m_ro_g_bc_types,
                                         bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                         bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                         bc_klo[lev]->array(), bc_khi[lev]->array(),
                                         geom[lev]);

          mol::compute_convective_rate(bx, ncomp, ccomp, dsdt, fx, fy, fz, geom[lev]);

        } // end density

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho*h_g)) -- the update for enthalpy
        // ***************************************************************************
        if (l_advect_enthalpy) {

          const int scomp = 0; // starting component
          const int ncomp = 1; // number of components
          const int ccomp = 1; // convection (dsdt) component

          mol::compute_convective_fluxes(bx, ncomp, scomp, rhohg,
                                         fx, fy, fz,
                                         ep_umac, ep_vmac, ep_wmac,
                                         bc_types, m_T_g_bc_types,
                                         bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                         bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                         bc_klo[lev]->array(), bc_khi[lev]->array(),
                                         geom[lev]);

          mol::compute_convective_rate(bx, ncomp, ccomp, dsdt, fx, fy, fz, geom[lev]);

        } // end enthalpy

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho*trac)) -- the update for tracer
        // ***************************************************************************
        if (l_advect_tracer) {

          const int scomp = 0; // starting component
          const int ncomp = 1; // number of components
          const int ccomp = 2; // convection (dsdt) component

          mol::compute_convective_fluxes(bx, ncomp, scomp, rhotrac,
                                         fx, fy, fz,
                                         ep_umac, ep_vmac, ep_wmac,
                                         bc_types, m_trac_g_bc_types,
                                         bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                         bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                         bc_klo[lev]->array(), bc_khi[lev]->array(),
                                         geom[lev]);

          mol::compute_convective_rate(bx, ncomp, ccomp, dsdt, fx, fy, fz, geom[lev]);


        } // end tracer

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho*X_g)) -- the update for species mass fraction
        // ***************************************************************************
        if (l_advect_species) {

          const int scomp = 0; // starting component
          const int ncomp = l_nspecies; // number of components
          const int ccomp = 0; // convection (dXdt) component

          mol::compute_convective_fluxes(bx, ncomp, scomp, rhoXgk,
                                         fx, fy, fz,
                                         ep_umac, ep_vmac, ep_wmac,
                                         bc_types, m_X_gk_bc_types,
                                         bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                         bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                         bc_klo[lev]->array(), bc_khi[lev]->array(),
                                         geom[lev]);

          mol::compute_convective_rate(bx, ncomp, ccomp, dXdt, fx, fy, fz, geom[lev]);


        } // end species

      }
    }// end MOL advection type

}
