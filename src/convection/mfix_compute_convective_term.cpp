#include <mfix.H>
#include <mfix_algorithm.H>
#include <hydro_utils.H>
#include <hydro_redistribution.H>

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
    m_iconserv_velocity.resize(3, 1);
    m_iconserv_velocity_d.resize(3, 1);
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

  const int l_nspecies = fluid.nspecies;

  if (advect_fluid_species) {
    m_iconserv_species.resize(l_nspecies, 1);
    m_iconserv_species_d.resize(l_nspecies, 1);
  }

}

//
// Compute the convection terms at all levels for all variables
//
void
mfix::mfix_compute_convective_term (Vector< MultiFab*      >& conv_u,  // velocity (3)
                                    Vector< MultiFab*      >& conv_s,  // density, tracers, enthalpy (2 + ntrac)
                                    Vector< MultiFab*      >& conv_X,  // species (nspecies)
                                    Vector< MultiFab*      > const& vel_forces,
                                    Vector< MultiFab*      > const& tra_forces,
                                    Vector< MultiFab const*> const& vel_in,
                                    Vector< MultiFab*      > const& ep_g_in,
                                    Vector< MultiFab const*> const& ro_g_in,
                                    Vector< MultiFab const*> const& h_g_in,
                                    Vector< MultiFab const*> const& trac_in,
                                    Vector< MultiFab const*> const& X_gk_in,
                                    Vector< MultiFab const*> const& txfr_in,
                                    Vector< MultiFab*      > const& ep_u_mac,
                                    Vector< MultiFab*      > const& ep_v_mac,
                                    Vector< MultiFab*      > const& ep_w_mac,
                                    Real l_dt, Real time)
{
    BL_PROFILE("mfix::mfix_compute_convective_term");

    const int ngmac = nghost_mac();

    const int l_nspecies = fluid.nspecies;

    int flux_comp, num_comp;

    bool fluxes_are_area_weighted = false;

    std::string advection_string;
    if (advection_type() == AdvectionType::Godunov)
        advection_string = "Godunov";
    else
        advection_string = "MOL";

    // Make one flux MF at each level to hold all the fluxes (velocity, density, tracers)
    // Note that we allocate data for all the possible variables but don't necessarily
    //      use that space if advect_density, etc not true
    int n_flux_comp = AMREX_SPACEDIM + 2 + ntrac + l_nspecies;

    // This will hold state on faces
    Vector<MultiFab> face_x(finest_level+1);
    Vector<MultiFab> face_y(finest_level+1);
    Vector<MultiFab> face_z(finest_level+1);

    // This will hold fluxes on faces
    Vector<MultiFab> flux_x(finest_level+1);
    Vector<MultiFab> flux_y(finest_level+1);
    Vector<MultiFab> flux_z(finest_level+1);

    Vector<Array<MultiFab*,AMREX_SPACEDIM> > fluxes(finest_level+1);
    Vector<Array<MultiFab*,AMREX_SPACEDIM> >  faces(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
       face_x[lev].define(ep_u_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),EBFactory(lev));
       face_y[lev].define(ep_v_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),EBFactory(lev));
       face_z[lev].define(ep_w_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),EBFactory(lev));

       flux_x[lev].define(ep_u_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),EBFactory(lev));
       flux_y[lev].define(ep_v_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),EBFactory(lev));
       flux_z[lev].define(ep_w_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),EBFactory(lev));

       faces[lev][0] = &face_x[lev];
       faces[lev][1] = &face_y[lev];
       faces[lev][2] = &face_z[lev];

       fluxes[lev][0] = &flux_x[lev];
       fluxes[lev][1] = &flux_y[lev];
       fluxes[lev][2] = &flux_z[lev];

       face_x[lev].setVal(0.);  face_y[lev].setVal(0.);  face_z[lev].setVal(0.);
       flux_x[lev].setVal(0.);  flux_y[lev].setVal(0.);  flux_z[lev].setVal(0.);
    }

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

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        if (ngmac > 0) {
            AMREX_D_TERM(ep_u_mac[lev]->FillBoundary(geom[lev].periodicity());,
                         ep_v_mac[lev]->FillBoundary(geom[lev].periodicity());,
                         ep_w_mac[lev]->FillBoundary(geom[lev].periodicity()););
        }

        MultiFab divu(vel_in[lev]->boxArray(),vel_in[lev]->DistributionMap(),1,4);
        divu.setVal(0.);
        Array<MultiFab const*, AMREX_SPACEDIM> u;
        AMREX_D_TERM(u[0] = ep_u_mac[lev];,
                     u[1] = ep_v_mac[lev];,
                     u[2] = ep_w_mac[lev];);

        const auto& ebfact = EBFactory(lev);

        if (!ebfact.isAllRegular())
            amrex::EB_computeDivergence(divu,u,geom[lev],true);
        else
            amrex::computeDivergence(divu,u,geom[lev]);

        divu.FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*ro_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();

            Array4<Real const> const& divu_arr = divu.const_array(mfi);

            // ************************************************************************
            // Velocity
            // ************************************************************************
            {
                int face_comp = 0;
                int ncomp = AMREX_SPACEDIM;
                bool is_velocity = true;
                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi, vel_in[lev]->const_array(mfi),
                                         AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                      flux_y[lev].array(mfi,face_comp),
                                                      flux_z[lev].array(mfi,face_comp)),
                                         AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                      face_y[lev].array(mfi,face_comp),
                                                      face_z[lev].array(mfi,face_comp)),
                                         AMREX_D_DECL(ep_u_mac[lev]->const_array(mfi),
                                                      ep_v_mac[lev]->const_array(mfi),
                                                      ep_w_mac[lev]->const_array(mfi)),
                                         divu_arr,
                                         (!vel_forces.empty()) ? vel_forces[lev]->const_array(mfi) : Array4<Real const>{},
                                         geom[lev], l_dt,
                                         get_hydro_velocity_bcrec(),
                                         get_hydro_velocity_bcrec_device_ptr(),
                                         get_velocity_iconserv_device_ptr(),
                                         ebfact,
                                         m_godunov_ppm, m_godunov_use_forces_in_trans,
                                         is_velocity, fluxes_are_area_weighted,
                                         advection_string);
            }

            // ************************************************************************
            // Density
            // ************************************************************************
            if (advect_density)
            {
                int face_comp = AMREX_SPACEDIM;
                int ncomp = 1;
                bool is_velocity = false;
                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi,
                                                        ro_g_in[lev]->const_array(mfi),
                                          AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                       flux_y[lev].array(mfi,face_comp),
                                                       flux_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                       face_y[lev].array(mfi,face_comp),
                                                       face_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(ep_u_mac[lev]->const_array(mfi),
                                                       ep_v_mac[lev]->const_array(mfi),
                                                       ep_w_mac[lev]->const_array(mfi)),
                                          divu_arr, Array4<Real const>{},
                                          geom[lev], l_dt,
                                          get_density_bcrec(),
                                          get_density_bcrec_device_ptr(),
                                          get_density_iconserv_device_ptr(),
                                          ebfact,
                                          m_godunov_ppm, m_godunov_use_forces_in_trans,
                                          is_velocity, fluxes_are_area_weighted,
                                          advection_string);
            }

            // ************************************************************************
            // (Rho*Enthalpy)
            // ************************************************************************
            // Make a FAB holding (rho * enthalpy) that is the same size as the original enthalpy FAB
            FArrayBox rhohfab;
            if (advect_enthalpy)
            {
                Box rhoh_box = Box((*h_g_in[lev])[mfi].box());
                Array4<Real> rhoh;
                Array4<Real const>   h =  h_g_in[lev]->const_array(mfi);
                Array4<Real const> rho = ro_g_in[lev]->const_array(mfi);
                rhohfab.resize(rhoh_box, 1);
                Elixir eli_rh = rhohfab.elixir();
                rhoh   = rhohfab.array();
                amrex::ParallelFor(rhoh_box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rhoh(i,j,k) = rho(i,j,k) * h(i,j,k);
                });

                int face_comp = 4;
                int ncomp = 1;
                bool is_velocity = false;

                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi, rhoh,
                                          AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                       flux_y[lev].array(mfi,face_comp),
                                                       flux_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                       face_y[lev].array(mfi,face_comp),
                                                       face_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(ep_u_mac[lev]->const_array(mfi),
                                                       ep_v_mac[lev]->const_array(mfi),
                                                       ep_w_mac[lev]->const_array(mfi)),
                                          divu_arr,
                                          Array4<Real const>{}, // enthalpy forces not defined
                                          geom[lev], l_dt,
                                          get_enthalpy_bcrec(),
                                          get_enthalpy_bcrec_device_ptr(),
                                          get_enthalpy_iconserv_device_ptr(),
                                          ebfact,
                                          m_godunov_ppm, m_godunov_use_forces_in_trans,
                                          is_velocity, fluxes_are_area_weighted,
                                          advection_string);
            }

            // ************************************************************************
            // (Rho*Tracer)
            // ************************************************************************
            // Make a FAB holding (rho * tracer) that is the same size as the original tracer FAB
            FArrayBox rhotracfab;
            if (advect_tracer && (ntrac>0)) {

                Box rhotrac_box = Box((*trac_in[lev])[mfi].box());
                Array4<Real const> tra = trac_in[lev]->const_array(mfi);
                Array4<Real const> rho = ro_g_in[lev]->const_array(mfi);
                rhotracfab.resize(rhotrac_box, ntrac);
                Elixir eli_rt  = rhotracfab.elixir();
                Array4<Real> rhotrac = rhotracfab.array();
                amrex::ParallelFor(rhotrac_box, ntrac,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhotrac(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
                });

                int face_comp = 5;
                int ncomp = ntrac;
                bool is_velocity = false;

                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi, rhotrac,
                                          AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                       flux_y[lev].array(mfi,face_comp),
                                                       flux_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                       face_y[lev].array(mfi,face_comp),
                                                       face_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(ep_u_mac[lev]->const_array(mfi),
                                                       ep_v_mac[lev]->const_array(mfi),
                                                       ep_w_mac[lev]->const_array(mfi)),
                                          divu_arr,
                                          (!tra_forces.empty()) ? tra_forces[lev]->const_array(mfi) : Array4<Real const>{},
                                          geom[lev], l_dt,
                                          get_tracer_bcrec(),
                                          get_tracer_bcrec_device_ptr(),
                                          get_tracer_iconserv_device_ptr(),
                                          ebfact,
                                          m_godunov_ppm, m_godunov_use_forces_in_trans,
                                          is_velocity, fluxes_are_area_weighted,
                                          advection_string);
            }

            // ************************************************************************
            // (Rho*Species)
            // ************************************************************************
            // Make a FAB holding (rho * species) that is the same size as the original tracer FAB
            FArrayBox rhoXfab;
            if (advect_fluid_species && (l_nspecies>0))
            {
                Box rhoX_box = Box((*X_gk_in[lev])[mfi].box());
                Array4<Real const>   X = X_gk_in[lev]->const_array(mfi);
                Array4<Real const> rho = ro_g_in[lev]->const_array(mfi);
                rhoXfab.resize(rhoX_box, l_nspecies);
                Elixir eli_rX  = rhoXfab.elixir();
                Array4<Real> rhoX = rhoXfab.array();
                amrex::ParallelFor(rhoX_box, l_nspecies,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhoX(i,j,k,n) = rho(i,j,k) * X(i,j,k,n);
                });

                int face_comp = 5+ntrac;
                int ncomp = l_nspecies;
                bool is_velocity = false;

                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi, rhoX,
                                          AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                       flux_y[lev].array(mfi,face_comp),
                                                       flux_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                       face_y[lev].array(mfi,face_comp),
                                                       face_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(ep_u_mac[lev]->const_array(mfi),
                                                       ep_v_mac[lev]->const_array(mfi),
                                                       ep_w_mac[lev]->const_array(mfi)),
                                          divu_arr,
                                          Array4<Real const>{}, // species forces not defined
                                          geom[lev], l_dt,
                                          get_species_bcrec(),
                                          get_species_bcrec_device_ptr(),
                                          get_species_iconserv_device_ptr(),
                                          ebfact,
                                          m_godunov_ppm, m_godunov_use_forces_in_trans,
                                          is_velocity, fluxes_are_area_weighted,
                                          advection_string);
            }
        } // mfi
    } // lev

    // In order to enforce conservation across coarse-fine boundaries we must be sure to average down the fluxes
    //    before we use them.  Note we also need to average down the face states if we are going to do
    //    convective differencing
    for (int lev = finest_level; lev > 0; --lev)
    {
        IntVect rr  = geom[lev].Domain().size() / geom[lev-1].Domain().size();
        EB_average_down_faces(GetArrOfConstPtrs( faces[lev]),  faces[lev-1], rr, geom[lev-1]);
        EB_average_down_faces(GetArrOfConstPtrs(fluxes[lev]), fluxes[lev-1], rr, geom[lev-1]);
    }

    MultiFab dXdt_tmp;

    int ngrow = 4;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        MultiFab dvdt_tmp(vel_in[lev]->boxArray(),dmap[lev],AMREX_SPACEDIM,ngrow,MFInfo(),EBFactory(lev));
        MultiFab dsdt_tmp(vel_in[lev]->boxArray(),dmap[lev],2+ntrac       ,ngrow,MFInfo(),EBFactory(lev));
        if (advect_fluid_species && l_nspecies > 0)
            dXdt_tmp.define(vel_in[lev]->boxArray(),dmap[lev],l_nspecies  ,ngrow,MFInfo(),EBFactory(lev));

        // Must initialize to zero because not all values may be set, e.g. outside the domain.
        dvdt_tmp.setVal(0.);
        dsdt_tmp.setVal(0.);
        if (advect_fluid_species && l_nspecies > 0)
            dXdt_tmp.setVal(0.);

        // We need to make sure the ghost cells are filled if using FluxRedistribution because
        // ep_g is used as the weights
        ep_g_in[lev]->FillBoundary();

        const auto& ebfact = EBFactory(lev);

        Real mult = -1.0;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*conv_u[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();

            EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
            bool regular = (flagfab.getType(amrex::grow(bx,2)) == FabType::regular);

            auto const& vfrac = ebfact.getVolFrac().const_array(mfi);

            if (flagfab.getType(bx) != FabType::covered)
            {
                if (!regular)
                {
                    {
                        int flux_comp = 0;
                        int num_comp = 3;
                        HydroUtils::EB_ComputeDivergence(bx, dvdt_tmp.array(mfi),
                                                         AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                                      flux_y[lev].const_array(mfi,flux_comp),
                                                                      flux_z[lev].const_array(mfi,flux_comp)),
                                                         vfrac, num_comp, geom[lev], mult, fluxes_are_area_weighted);
                    }

                    {
                       int flux_comp = 3;
                       int  num_comp = 2+ntrac;
                       HydroUtils::EB_ComputeDivergence(bx, dsdt_tmp.array(mfi),
                                                        AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                                     flux_y[lev].const_array(mfi,flux_comp),
                                                                     flux_z[lev].const_array(mfi,flux_comp)),
                                                        vfrac, num_comp, geom[lev], mult, fluxes_are_area_weighted);
                    }

                    if (advect_fluid_species && l_nspecies > 0)
                    {
                        int flux_comp = 5+ntrac;
                        int  num_comp = l_nspecies;
                        HydroUtils::EB_ComputeDivergence(bx, dXdt_tmp.array(mfi),
                                                         AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                                      flux_y[lev].const_array(mfi,flux_comp),
                                                                      flux_z[lev].const_array(mfi,flux_comp)),
                                                         vfrac, num_comp, geom[lev], mult, fluxes_are_area_weighted);
                    }
                }
                else
                {
                    // Velocity
                    {
                        int flux_comp = 0;
                        int num_comp = 3;
                        HydroUtils::ComputeDivergence(bx, dvdt_tmp.array(mfi),
                                                  AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                               flux_y[lev].const_array(mfi,flux_comp),
                                                               flux_z[lev].const_array(mfi,flux_comp)),
                                                  AMREX_D_DECL(face_x[lev].const_array(mfi,flux_comp),
                                                               face_y[lev].const_array(mfi,flux_comp),
                                                               face_z[lev].const_array(mfi,flux_comp)),
                                                  AMREX_D_DECL(ep_u_mac[lev]->const_array(mfi),
                                                               ep_v_mac[lev]->const_array(mfi),
                                                               ep_w_mac[lev]->const_array(mfi)),
                                                  num_comp, geom[lev],
                                                  get_velocity_iconserv_device_ptr(), mult,
                                                  fluxes_are_area_weighted);
                    }

                    // Density, enthalpy, tracers
                    // NOTE: this isn't fully general because we pass in the density "iconserv"
                    //       flag but in this case density, (rho h) and (rho T) are all viewed as conservative
                    {
                        int flux_comp = 3;
                        int  num_comp = 2 + ntrac;
                        HydroUtils::ComputeDivergence(bx, dsdt_tmp.array(mfi),
                                                  AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                               flux_y[lev].const_array(mfi,flux_comp),
                                                               flux_z[lev].const_array(mfi,flux_comp)),
                                                  AMREX_D_DECL(face_x[lev].const_array(mfi,flux_comp),
                                                               face_y[lev].const_array(mfi,flux_comp),
                                                               face_z[lev].const_array(mfi,flux_comp)),
                                                  AMREX_D_DECL(ep_u_mac[lev]->const_array(mfi),
                                                               ep_v_mac[lev]->const_array(mfi),
                                                               ep_w_mac[lev]->const_array(mfi)),
                                                  num_comp, geom[lev],
                                                  get_density_iconserv_device_ptr(), mult,
                                                  fluxes_are_area_weighted);
                    }

                    // Species
                    if (advect_fluid_species && l_nspecies > 0)
                    {
                        int flux_comp = 5+ntrac;
                        int  num_comp = l_nspecies;
                        HydroUtils::ComputeDivergence(bx, dXdt_tmp.array(mfi),
                                                  AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                               flux_y[lev].const_array(mfi,flux_comp),
                                                               flux_z[lev].const_array(mfi,flux_comp)),
                                                  AMREX_D_DECL(face_x[lev].const_array(mfi,flux_comp),
                                                               face_y[lev].const_array(mfi,flux_comp),
                                                               face_z[lev].const_array(mfi,flux_comp)),
                                                  AMREX_D_DECL(ep_u_mac[lev]->const_array(mfi),
                                                               ep_v_mac[lev]->const_array(mfi),
                                                               ep_w_mac[lev]->const_array(mfi)),
                                                  num_comp, geom[lev],
                                                  get_species_iconserv_device_ptr(), mult,
                                                  fluxes_are_area_weighted);
                    } // species
                }
            }
        } // mfi

        // We only filled these on the valid cells so we fill same-level interior ghost cells here.
        // (We don't need values outside the domain or at a coarser level so we can call just FillBoundary)
        dvdt_tmp.FillBoundary(geom[lev].periodicity());
        dsdt_tmp.FillBoundary(geom[lev].periodicity());
        if (advect_fluid_species && l_nspecies > 0)
            dXdt_tmp.FillBoundary(geom[lev].periodicity());

        int max_ncomp = std::max(std::max(AMREX_SPACEDIM,l_nspecies),2+ntrac);

        for (MFIter mfi(*ro_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.tilebox();

          EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
          bool regular = (flagfab.getType(amrex::grow(bx,2)) == FabType::regular);

          if (flagfab.getType(bx) != FabType::covered)
          {
           if (!regular)
           {
            // Make a FAB holding (rho * tracer) that is the same size as the original tracer FAB
            Box tmp_box = Box((*ro_g_in[lev])[mfi].box());
            Array4<Real const> rho = ro_g_in[lev]->const_array(mfi);

            FArrayBox scratchfab(amrex::grow(bx,ngrow), max_ncomp);
            Elixir eli_scratch = scratchfab.elixir();

            auto const& vfrac = ebfact.getVolFrac().const_array(mfi);
            auto const& ccc   = ebfact.getCentroid().const_array(mfi);

            auto const& apx = ebfact.getAreaFrac()[0]->const_array(mfi);
            auto const& apy = ebfact.getAreaFrac()[1]->const_array(mfi);
            auto const& apz = ebfact.getAreaFrac()[2]->const_array(mfi);

            auto const& fcx = ebfact.getFaceCent()[0]->const_array(mfi);
            auto const& fcy = ebfact.getFaceCent()[1]->const_array(mfi);
            auto const& fcz = ebfact.getFaceCent()[2]->const_array(mfi);

            // Velocity
            {
                int ncomp = 3;

                auto const& bc_vel = get_hydro_velocity_bcrec_device_ptr();

                Redistribution::Apply( bx, ncomp, conv_u[lev]->array(mfi), dvdt_tmp.array(mfi),
                                       vel_in[lev]->const_array(mfi),
                                       (m_redistribution_type == "StateRedist") ? scratchfab.array()
                                                                                : ep_g_in[lev]->array(mfi),
                                       flagfab.const_array(),
                                       apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                       bc_vel, geom[lev], l_dt, m_redistribution_type);
            }

            // Density
            if (advect_density)
            {
                int ncomp = 1;

                auto const& bc_den = get_density_bcrec_device_ptr();

                Redistribution::Apply( bx, ncomp, conv_s[lev]->array(mfi,0), dsdt_tmp.array(mfi,0),
                                       ro_g_in[lev]->const_array(mfi),
                                       (m_redistribution_type == "StateRedist") ? scratchfab.array()
                                                                                : ep_g_in[lev]->array(mfi),
                                       flagfab.const_array(),
                                       apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                       bc_den, geom[lev], l_dt, m_redistribution_type);
            }


            // Enthalpy
            if (advect_enthalpy)
            {
                int ncomp = 1;
                Array4<Real const>   h =  h_g_in[lev]->const_array(mfi);
                FArrayBox rhohfab;
                rhohfab.resize(tmp_box, 1);
                Elixir eli_rh  = rhohfab.elixir();
                Array4<Real> rhoh = rhohfab.array();
                amrex::ParallelFor(tmp_box, ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhoh(i,j,k,n) = rho(i,j,k) * h(i,j,k,n);
                });

                auto const& bc_rh = get_enthalpy_bcrec_device_ptr();

                Redistribution::Apply( bx, ncomp, conv_s[lev]->array(mfi,1), dsdt_tmp.array(mfi,1), rhoh,
                                       (m_redistribution_type == "StateRedist") ? scratchfab.array()
                                                                                : ep_g_in[lev]->array(mfi),
                                       flagfab.const_array(),
                                       apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                       bc_rh, geom[lev], l_dt, m_redistribution_type);
            }

            // Tracers
            if (advect_tracer && (ntrac > 0)) {
                int ncomp = ntrac;

                Array4<Real const> tra = trac_in[lev]->const_array(mfi);
                FArrayBox rhotracfab;
                rhotracfab.resize(tmp_box, ntrac);
                Elixir eli_rt  = rhotracfab.elixir();
                Array4<Real> rhotrac = rhotracfab.array();
                amrex::ParallelFor(tmp_box, ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhotrac(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
                });

                auto const& bc_rt = get_tracer_bcrec_device_ptr();

                Redistribution::Apply( bx, ncomp, conv_s[lev]->array(mfi,2), dsdt_tmp.array(mfi,2), rhotrac,
                                       (m_redistribution_type == "StateRedist") ? scratchfab.array()
                                                                                : ep_g_in[lev]->array(mfi),
                                       flagfab.const_array(),
                                       apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                       bc_rt, geom[lev], l_dt, m_redistribution_type);
            }

            if (advect_fluid_species && (l_nspecies > 0))
            {
                int ncomp = l_nspecies;

                Array4<Real const>   X = X_gk_in[lev]->const_array(mfi);
                FArrayBox rhoXfab;
                rhoXfab.resize(tmp_box, l_nspecies);
                Elixir eli_rX  = rhoXfab.elixir();
                Array4<Real> rhoX = rhoXfab.array();
                amrex::ParallelFor(tmp_box, ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhoX(i,j,k,n) = rho(i,j,k) * X(i,j,k,n);
                });

                auto const& bc_rX = get_species_bcrec_device_ptr();

                Redistribution::Apply( bx, ncomp, conv_X[lev]->array(mfi), dXdt_tmp.array(mfi), rhoX,
                                       (m_redistribution_type == "StateRedist") ? scratchfab.array()
                                                                                : ep_g_in[lev]->array(mfi),
                                       flagfab.const_array(),
                                       apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                       bc_rX, geom[lev], l_dt, m_redistribution_type);
            }

         } // test on if regular
         else
         {
             Array4<Real      > conv_u_arr = conv_u[lev]->array(mfi);
             Array4<Real const>   dvdt_arr = dvdt_tmp.array(mfi);
             int ncomp = AMREX_SPACEDIM;
             amrex::ParallelFor(bx,ncomp,
             [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
             {
                 conv_u_arr(i,j,k,n) = dvdt_arr(i,j,k,n);
             });

             Array4<Real      > conv_s_arr = conv_s[lev]->array(mfi);
             Array4<Real const>   dsdt_arr = dsdt_tmp.array(mfi);
             int ncomp = conv_s[lev]->nComp();
             amrex::ParallelFor(bx,ncomp,
             [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
             {
                 conv_s_arr(i,j,k,n) = dsdt_arr(i,j,k,n);
             });

             if (advect_fluid_species && (l_nspecies > 0))
             {
                 Array4<Real      > conv_X_arr = conv_X[lev]->array(mfi);
                 Array4<Real const>   dXdt_arr = dXdt_tmp.array(mfi);
                 int ncomp = conv_X[lev]->nComp();
                 amrex::ParallelFor(bx,ncomp,
                 [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                 {
                     conv_X_arr(i,j,k,n) = dXdt_arr(i,j,k,n);
                 });
             }
         }
        } // test on if covered
       } // mfi
    } // lev
}

