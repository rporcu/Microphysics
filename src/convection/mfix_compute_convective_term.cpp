#include <mfix.H>
#include <MOL.H>
#include <Godunov.H>
#include <EBGodunov.H>
#include <Redistribution.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_EB_utils.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>

void mfix::init_advection ()
{
  if (advection_type() == AdvectionType::Godunov) {
    m_iconserv_velocity.resize(3, 1);
    m_iconserv_velocity_d.resize(3, 1);
  } else {
    m_iconserv_velocity.resize(3, 0);
    m_iconserv_velocity_d.resize(3, 0);
  }

  m_iconserv_density.resize(1, 1);
  m_iconserv_density_d.resize(1, 1);

  if (advect_enthalpy) {
    m_iconserv_enthalpy.resize(1, 1);
    m_iconserv_enthalpy_d.resize(1, 1);
  }

  if (advect_tracer) {
    m_iconserv_tracer.resize(ntrac, 1);
    m_iconserv_tracer_d.resize(ntrac, 1);
  }

  const int l_nspecies = FLUID::nspecies;

  if (advect_fluid_species) {
    m_iconserv_species.resize(l_nspecies, 1);
    m_iconserv_species_d.resize(l_nspecies, 1);
  }

}


//
// Compute the three components of the convection term
//
void
mfix::mfix_compute_convective_term (Vector< MultiFab*      >& conv_u_in,
                                    Vector< MultiFab*      >& conv_s_in,
                                    Vector< MultiFab*      >& conv_X_in,
                                    Vector< MultiFab*      > const& vel_forces,
                                    Vector< MultiFab*      > const& tra_forces,
                                    Vector< MultiFab const*> const& vel_in,
                                    Vector< MultiFab const*> const& ep_g_in,
                                    Vector< MultiFab const*> const& ro_g_in,
                                    Vector< MultiFab const*> const& h_g_in,
                                    Vector< MultiFab const*> const& trac_in,
                                    Vector< MultiFab const*> const& X_gk_in,
                                    Vector< MultiFab const*> const& txfr_in,
                                    Vector< MultiFab*      > const& ep_u_mac,
                                    Vector< MultiFab*      > const& ep_v_mac,
                                    Vector< MultiFab*      > const& ep_w_mac,
                                    Vector< MultiFab const*> const& rhs_mac,
                                    Vector< MultiFab      *> const& divtau_old,
                                    Real l_dt, Real time)
{
    BL_PROFILE("mfix::mfix_compute_convective_term");

    int ngmac = nghost_mac();

    const int l_nspecies = FLUID::nspecies;

    // We first compute the velocity forcing terms to be used in predicting
    //    to faces before the MAC projection
    if (advection_type() != AdvectionType::MOL) {

      bool include_pressure_gradient = !(m_use_mac_phi_in_godunov);
      bool include_drag_force = include_pressure_gradient && m_use_drag_in_godunov;
      compute_vel_forces(vel_forces, vel_in, ro_g_in, txfr_in,
         include_pressure_gradient, include_drag_force);

      if (m_godunov_include_diff_in_forcing)
        for (int lev = 0; lev <= finest_level; ++lev)
          MultiFab::Add(*vel_forces[lev], *m_leveldata[lev]->divtau_o, 0, 0, 3, 0);

      if (nghost_force() > 0)
        fillpatch_force(time, vel_forces, nghost_force());

    }

    // Do projection on all AMR levels in one shot -- note that the {u_mac, v_mac, w_mac}
    //    arrays returned from this call are in fact {ep * u_mac, ep * v_mac, ep * w_mac}
    //    on face CENTROIDS
    compute_MAC_projected_velocities(time, l_dt, vel_in, ep_u_mac, ep_v_mac, ep_w_mac,
                                     ep_g_in, ro_g_in, vel_forces, rhs_mac);


    // We now re-compute the velocity forcing terms including the pressure gradient,
    //    and compute the tracer forcing terms for the first time
    if (advection_type() != AdvectionType::MOL) {

      bool include_pressure_gradient = true;
      bool include_drag_force = true;
      compute_vel_forces(vel_forces, vel_in, ro_g_in, txfr_in,
         include_pressure_gradient, include_drag_force);

      if (m_godunov_include_diff_in_forcing)
        for (int lev = 0; lev <= finest_level; ++lev)
          MultiFab::Add(*vel_forces[lev], *m_leveldata[lev]->divtau_o, 0, 0, 3, 0);

      if (nghost_force() > 0)
        fillpatch_force(time, vel_forces, nghost_force());

#if 0
      // TODO
      if(advect_enthalpy){
        amrex::Abort("Enthalpy forces are not broken out yet.");
      }
      if(advect_fluid_species){
        amrex::Abort("Species forces are not broken out yet.");
      }

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

      MultiFab divu(vel_in[lev]->boxArray(),vel_in[lev]->DistributionMap(),1,4);
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
      for (MFIter mfi(*ep_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();
        mfix::compute_convective_term(bx, lev, l_dt, mfi,
                                      conv_u_in[lev]->array(mfi),
                                      conv_s_in[lev]->array(mfi),
                                      (advect_fluid_species) ? conv_X_in[lev]->array(mfi) : Array4<Real>{},
                                      advect_density, advect_enthalpy, advect_tracer, advect_fluid_species,
                                      vel_in[lev]->const_array(mfi),
                                      ep_g_in[lev]->const_array(mfi),
                                      ro_g_in[lev]->const_array(mfi),
                                      (advect_enthalpy) ?  h_g_in[lev]->const_array(mfi) : Array4<Real const>{},
                                      (advect_tracer) ? trac_in[lev]->const_array(mfi) : Array4<Real const>{},
                                      (advect_fluid_species) ? X_gk_in[lev]->const_array(mfi) : Array4<Real const>{},
                                      l_nspecies,
                                      divu.const_array(mfi),
                                      ep_u_mac[lev]->const_array(mfi),
                                      ep_v_mac[lev]->const_array(mfi),
                                      ep_w_mac[lev]->const_array(mfi),
                                      (!vel_forces.empty()) ? vel_forces[lev]->const_array(mfi)
                                      : Array4<Real const>{},
                                      (!tra_forces.empty()) ? tra_forces[lev]->const_array(mfi)
                                      : Array4<Real const>{});
      } // mfi

    } // lev

}





void
mfix::compute_convective_term (Box const& bx, int lev, const Real l_dt, MFIter const& mfi,
                               Array4<Real> const& dvdt,      // velocity
                               Array4<Real> const& dsdt,      // density, enthalpy, tracer
                               Array4<Real> const& dXdt,      // tracer
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
                               Array4<Real const> const& fvel,
                               Array4<Real const> const& ftra)
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

    FArrayBox rhohgfab, rhotracfab, rhoXgkfab;
    Elixir eli_rh, eli_rt, eli_rXgk;
    Array4<Real> rhohg, rhotrac, rhoXgk;

    // Make FABs holding (rho * hg), (rho * tracer), and (rho * Xgk)
    // that are the same size as the original FABs
    if (l_advect_enthalpy) {
      Box rhohg_box(hg);
      rhohgfab.resize(rhohg_box, 1);
      eli_rh = rhohgfab.elixir();
      rhohg = rhohgfab.array();
      amrex::ParallelFor(rhohg_box, [rho, hg, rhohg]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        rhohg(i,j,k) = rho(i,j,k) * hg(i,j,k);
      });
    }

    if (l_advect_tracer) {
      Box rhotrac_box(tra);
      rhotracfab.resize(rhotrac_box, ntrac);
      eli_rt = rhotracfab.elixir();
      rhotrac = rhotracfab.array();
      amrex::ParallelFor(rhotrac_box, ntrac, [rho, tra, rhotrac]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        rhotrac(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
      });
    }

    if (l_advect_species) {
      Box rhoXgk_box(Xgk);
      rhoXgkfab.resize(rhoXgk_box, l_nspecies);
      eli_rXgk = rhoXgkfab.elixir();
      rhoXgk = rhoXgkfab.array();
      amrex::ParallelFor(rhoXgk_box, l_nspecies, [rho, Xgk, rhoXgk]
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

    if (advection_type() == AdvectionType::Godunov) {

      int n_tmp_fac  = 14;

      // n_tmp_grow is used to create tmpfab which is passed in to compute_godunov_advection
      // (both regular and EB) as pointer "p" and is used to hold Imx/Ipx etc ...
      int n_tmp_grow;
      if (m_redistribution_type == "StateRedist")
        n_tmp_grow = 5;
      else
        n_tmp_grow = 4;

      FArrayBox tmpfab(amrex::grow(bx,n_tmp_grow), nmaxcomp*n_tmp_fac+1);
      Elixir eli = tmpfab.elixir();

      Box gbx = bx;
      if (!regular) {
        // We need 3 if we are doing state redistribution
        if (m_redistribution_type == "StateRedist")
          gbx.grow(3);
        else if (m_redistribution_type == "FluxRedist")
          gbx.grow(2);
        else if (m_redistribution_type == "NoRedist")
          gbx.grow(1);
        else
          amrex::Abort("Dont know this redistribution type");
      }
      // This one holds the convective term on a grown region so we can redistribute
      FArrayBox dUdt_tmpfab(gbx,nmaxcomp);
      Array4<Real> dUdt_tmp = dUdt_tmpfab.array();
      Elixir eli_du = dUdt_tmpfab.elixir();

      FArrayBox scratch_fab(gbx,3*nmaxcomp);
      Array4<Real> scratch = scratch_fab.array();
      Elixir eli_scratch = scratch_fab.elixir();

      if (!regular) {

        // ***************************************************************************
        // Compute div (ep_u_mac * (u)) -- the update for velocity
        // ***************************************************************************
        {

          const int ncomp = 3; // number of components
          const int ccomp = 0; // convection (dsdt) component

          ebgodunov::compute_godunov_advection(gbx, ncomp,
                                               dUdt_tmp, vel,
                                               ep_umac, ep_vmac, ep_wmac,
                                               fvel, divu, l_dt,
                                               get_velocity_bcrec(),
                                               get_velocity_bcrec_device_ptr(),
                                               get_velocity_iconserv_device_ptr(),
                                               tmpfab.dataPtr(), flag,
                                               apx, apy, apz, vfrac,
                                               fcx, fcy, fcz, ccc,
                                               geom[lev], true); // is_velocity

          redistribution::redistribute_eb(bx, ncomp, ccomp, dvdt, dUdt_tmp, vel, scratch,
                                          ep_g, flag, apx, apy, apz, vfrac,
                                          fcx, fcy, fcz, ccc, geom[lev],
                                          l_dt, m_redistribution_type);
        } //end velocity


        // ***************************************************************************
        // Compute div (ep_u_mac * (rho)) -- the update for density
        // ***************************************************************************
        if (l_advect_density) {

          const int ncomp = 1; // number of components
          const int ccomp = 0; // convection (dsdt) component

          ebgodunov::compute_godunov_advection(gbx, ncomp,
                                               dUdt_tmp, rho,
                                               ep_umac, ep_vmac, ep_wmac,
                                               {}, divu, l_dt,
                                               get_density_bcrec(),
                                               get_density_bcrec_device_ptr(),
                                               get_density_iconserv_device_ptr(),
                                               tmpfab.dataPtr(), flag,
                                               apx, apy, apz, vfrac,
                                               fcx, fcy, fcz, ccc,
                                               geom[lev]);

          redistribution::redistribute_eb(bx, ncomp, ccomp, dsdt, dUdt_tmp, rho, scratch,
                                          ep_g, flag, apx, apy, apz, vfrac,
                                          fcx, fcy, fcz, ccc, geom[lev],
                                          l_dt, m_redistribution_type);
        } // end density

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho*h_g)) -- the update for enthalpy
        // ***************************************************************************
        if (l_advect_enthalpy) {

          const int ncomp = 1; // number of components
          const int ccomp = 1; // convection (dsdt) component

          ebgodunov::compute_godunov_advection(gbx, ncomp,
                                               dUdt_tmp, rhohg,
                                               ep_umac, ep_vmac, ep_wmac,
                                               {}, divu, l_dt,
                                               get_enthalpy_bcrec(),
                                               get_enthalpy_bcrec_device_ptr(),
                                               get_enthalpy_iconserv_device_ptr(),
                                               tmpfab.dataPtr(), flag,
                                               apx, apy, apz, vfrac,
                                               fcx, fcy, fcz, ccc,
                                               geom[lev]);

          redistribution::redistribute_eb(bx, ncomp, ccomp, dsdt, dUdt_tmp, rhohg, scratch,
                                          ep_g, flag, apx, apy, apz, vfrac,
                                          fcx, fcy, fcz, ccc, geom[lev],
                                          l_dt, m_redistribution_type);
        } // end enthalpy

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho*trac)) -- the update for tracer
        // ***************************************************************************
        if (l_advect_tracer) {

          const int ncomp = 1; // number of components
          const int ccomp = 2; // convection (dsdt) component

          ebgodunov::compute_godunov_advection(gbx, ncomp,
                                               dUdt_tmp, rhotrac,
                                               ep_umac, ep_vmac, ep_wmac,
                                               {}, divu, l_dt,
                                               get_tracer_bcrec(),
                                               get_tracer_bcrec_device_ptr(),
                                               get_tracer_iconserv_device_ptr(),
                                               tmpfab.dataPtr(), flag,
                                               apx, apy, apz, vfrac,
                                               fcx, fcy, fcz, ccc,
                                               geom[lev]);

          redistribution::redistribute_eb(bx, ncomp, ccomp, dsdt, dUdt_tmp, rhotrac, scratch,
                                          ep_g, flag, apx, apy, apz, vfrac,
                                          fcx, fcy, fcz, ccc, geom[lev],
                                          l_dt, m_redistribution_type);
        } // end tracer

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho*X_g)) -- the update for species mass fraction
        // ***************************************************************************
        if (l_advect_species) {

          const int ncomp = l_nspecies; // number of components
          const int ccomp = 0; // convection (dXdt) component

          ebgodunov::compute_godunov_advection(gbx, ncomp,
                                               dUdt_tmp, rhoXgk,
                                               ep_umac, ep_vmac, ep_wmac,
                                               {}, divu, l_dt,
                                               get_species_bcrec(),
                                               get_species_bcrec_device_ptr(),
                                               get_species_iconserv_device_ptr(),
                                               tmpfab.dataPtr(), flag,
                                               apx, apy, apz, vfrac,
                                               fcx, fcy, fcz, ccc,
                                               geom[lev]);

          redistribution::redistribute_eb(bx, ncomp, ccomp, dXdt, dUdt_tmp, rhoXgk, scratch,
                                          ep_g, flag, apx, apy, apz, vfrac,
                                          fcx, fcy, fcz, ccc, geom[lev],
                                          l_dt, m_redistribution_type);
        } // end species

        Gpu::streamSynchronize();

      } else {

        // ***************************************************************************
        // Compute div (ep_u_mac * (u)) -- the update for velocity
        // ***************************************************************************
        {

          const int ncomp = 3; // number of components
          const int ccomp = 0; // convection (dsdt) component

          godunov::compute_godunov_advection(bx, ncomp, ccomp, dvdt, vel,
                                             ep_umac, ep_vmac, ep_wmac,
                                             fvel, divu, l_dt,
                                             get_velocity_bcrec_device_ptr(),
                                             get_velocity_iconserv_device_ptr(),
                                             tmpfab.dataPtr(),m_godunov_ppm,
                                             m_godunov_use_forces_in_trans,
                                             geom[lev], true);

        } //end velocity


        // ***************************************************************************
        // Compute div (ep_u_mac * (rho)) -- the update for density
        // ***************************************************************************
        if (l_advect_density) {

          const int ncomp = 1; // number of components
          const int ccomp = 0; // convection (dsdt) component

          godunov::compute_godunov_advection(bx, ncomp, ccomp, dsdt, rho,
                                             ep_umac, ep_vmac, ep_wmac,
                                             fvel, divu, l_dt,
                                             get_density_bcrec_device_ptr(),
                                             get_density_iconserv_device_ptr(),
                                             tmpfab.dataPtr(),m_godunov_ppm,
                                             m_godunov_use_forces_in_trans,
                                             geom[lev], true);


        } // end density

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho*h_g)) -- the update for enthalpy
        // ***************************************************************************
        if (l_advect_enthalpy) {

          const int ncomp = 1; // number of components
          const int ccomp = 1; // convection (dsdt) component

          godunov::compute_godunov_advection(bx, ncomp, ccomp, dsdt, rhohg,
                                             ep_umac, ep_vmac, ep_wmac,
                                             fvel, divu, l_dt,
                                             get_enthalpy_bcrec_device_ptr(),
                                             get_enthalpy_iconserv_device_ptr(),
                                             tmpfab.dataPtr(),m_godunov_ppm,
                                             m_godunov_use_forces_in_trans,
                                             geom[lev], true);

        } // end enthalpy

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho*trac)) -- the update for tracer
        // ***************************************************************************
        if (l_advect_tracer) {

          const int ncomp = 1; // number of components
          const int ccomp = 2; // convection (dsdt) component

          godunov::compute_godunov_advection(bx, ncomp, ccomp, dsdt, rhotrac,
                                             ep_umac, ep_vmac, ep_wmac,
                                             fvel, divu, l_dt,
                                             get_tracer_bcrec_device_ptr(),
                                             get_tracer_iconserv_device_ptr(),
                                             tmpfab.dataPtr(),m_godunov_ppm,
                                             m_godunov_use_forces_in_trans,
                                             geom[lev], true);
        } // end tracer

        // ***************************************************************************
        // Compute div (ep_u_mac * (rho*X_g)) -- the update for species mass fraction
        // ***************************************************************************
        if (l_advect_species) {

          const int ncomp = l_nspecies; // number of components
          const int ccomp = 0; // convection (dXdt) component

          godunov::compute_godunov_advection(bx, ncomp, ccomp, dXdt, rhoXgk,
                                             ep_umac, ep_vmac, ep_wmac,
                                             fvel, divu, l_dt,
                                             get_species_bcrec_device_ptr(),
                                             get_species_iconserv_device_ptr(),
                                             tmpfab.dataPtr(),m_godunov_ppm,
                                             m_godunov_use_forces_in_trans,
                                             geom[lev], true);
        } // end species

        Gpu::streamSynchronize();
      }


      /**************************************************************************
       *                                                                        *
       *                       Method of Lines (MOL)                            *
       *                                                                        *
       **************************************************************************/
    } else if (advection_type() == AdvectionType::MOL) {

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

        Array4<Real> scratch = tmpfab.array(0);
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

          redistribution::redistribute_eb(bx, ncomp, ccomp, dvdt, dUdt_tmp, vel, scratch,
                                          ep_g, flag, apx, apy, apz, vfrac,
                                          fcx, fcy, fcz, ccc, geom[lev],
                                          l_dt, m_redistribution_type);


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

          redistribution::redistribute_eb(bx, ncomp, ccomp, dsdt, dUdt_tmp, rho, scratch,
                                          ep_g, flag, apx, apy, apz, vfrac,
                                          fcx, fcy, fcz, ccc, geom[lev],
                                          l_dt, m_redistribution_type);
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

          redistribution::redistribute_eb(bx, ncomp, ccomp, dsdt, dUdt_tmp, rhohg, scratch,
                                          ep_g, flag, apx, apy, apz, vfrac,
                                          fcx, fcy, fcz, ccc, geom[lev],
                                          l_dt, m_redistribution_type);
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

          redistribution::redistribute_eb(bx, ncomp, ccomp, dsdt, dUdt_tmp, rhotrac, scratch,
                                          ep_g, flag, apx, apy, apz, vfrac,
                                          fcx, fcy, fcz, ccc, geom[lev],
                                          l_dt, m_redistribution_type);
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

          redistribution::redistribute_eb(bx, ncomp, ccomp, dXdt, dUdt_tmp, rhoXgk, scratch,
                                          ep_g, flag, apx, apy, apz, vfrac,
                                          fcx, fcy, fcz, ccc, geom[lev],
                                          l_dt, m_redistribution_type);

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

      Gpu::streamSynchronize();

    }// end MOL advection type

}
