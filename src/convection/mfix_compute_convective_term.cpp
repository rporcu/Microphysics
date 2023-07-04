#include <mfix.H>
#include <mfix_algorithm.H>
#include <hydro_utils.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_EB_utils.H>
#include <AMReX_EB_Redistribution.H>

#include <mfix_eb.H>
#include <mfix_fluid.H>
#include <mfix_species.H>

void mfix::init_advection ()
{
  // Update velocity with convective differencing
  m_iconserv_velocity.resize(3, 0);
  m_iconserv_velocity_d.resize(3, 0);

  // Update density with conservative differencing
  m_iconserv_density.resize(1, 1);
  m_iconserv_density_d.resize(1, 1);

  // Update (rho h) with conservative differencing
  if (fluid.solve_enthalpy()) {
    m_iconserv_enthalpy.resize(1, 1);
    m_iconserv_enthalpy_d.resize(1, 1);
  }

  // Update (rho T) with conservative differencing
  if (fluid.solve_tracer()) {
    m_iconserv_tracer.resize(ntrac, 1);
    m_iconserv_tracer_d.resize(ntrac, 1);
  }

  // Update (rho X) with conservative differencing
  if (fluid.solve_species()) {
    m_iconserv_species.resize(fluid.nspecies(), 1);
    m_iconserv_species_d.resize(fluid.nspecies(), 1);
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
                                    Vector< MultiFab const*> const& eb_vel,
                                    Vector< MultiFab const*> const& eb_scalars,
                                    Vector< MultiFab const*> const& eb_species,
                                    Vector< MultiFab*      > const& ep_u_mac,
                                    Vector< MultiFab*      > const& ep_v_mac,
                                    Vector< MultiFab*      > const& ep_w_mac,
                                    Real l_dt, Real time)
{
    BL_PROFILE("mfix::mfix_compute_convective_term");

    int flux_comp, num_comp;

    bool fluxes_are_area_weighted = false;
    bool knownFaceStates          = false; // HydroUtils always recompute face states

    std::string advection_string;
    if (advection_type() == AdvectionType::Godunov)
        advection_string = "Godunov";
    else
        advection_string = "MOL";

    // Make one flux MF at each level to hold all the fluxes (velocity, density, tracers)
    // Note that we allocate data for all the possible variables but don't necessarily
    //      use that space if fluid.solve_density, etc not true
    int n_flux_comp = AMREX_SPACEDIM + 2 + ntrac + fluid.nspecies();

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

    Vector<MultiFab> divu(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
       face_x[lev].define(ep_u_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),EBFactory(lev));
       face_y[lev].define(ep_v_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),EBFactory(lev));
       face_z[lev].define(ep_w_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),EBFactory(lev));

       flux_x[lev].define(ep_u_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),EBFactory(lev));
       flux_y[lev].define(ep_v_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),EBFactory(lev));
       flux_z[lev].define(ep_w_mac[lev]->boxArray(),dmap[lev],n_flux_comp,0,MFInfo(),EBFactory(lev));

       divu[lev].define(vel_in[lev]->boxArray(),dmap[lev],1,4,MFInfo(),EBFactory(lev));

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
      if(fluid.solve_enthalpy()){
        amrex::Abort("Enthalpy forces are not broken out yet.");
      }
      if(solve_species){
        amrex::Abort("Species forces are not broken out yet.");
      }

      // Note this is forcing for (rho s), not for s
      if (fluid.solve_tracer) {
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
        divu[lev].setVal(0.);
        Array<MultiFab const*, AMREX_SPACEDIM> u;
        AMREX_D_TERM(u[0] = ep_u_mac[lev];,
                     u[1] = ep_v_mac[lev];,
                     u[2] = ep_w_mac[lev];);

        const auto& ebfact = EBFactory(lev);

        if (!ebfact.isAllRegular()) {
            if (m_embedded_boundaries.has_flow())
               amrex::EB_computeDivergence(divu[lev],u,geom[lev],true, *eb_vel[lev] );
            else
               amrex::EB_computeDivergence(divu[lev],u,geom[lev],true);
        } else {
            amrex::computeDivergence(divu[lev],u,geom[lev]);
        }

        divu[lev].FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel private (flux_comp, num_comp) if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*ro_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();

            Array4<Real const> const& divu_arr = divu[lev].const_array(mfi);

            // ************************************************************************
            // Velocity
            // ************************************************************************
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
                                     knownFaceStates,
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
                                     (m_embedded_boundaries.has_flow()) ? eb_vel[lev]->const_array(mfi) : Array4<Real const>{},
                                     m_godunov_ppm, m_godunov_use_forces_in_trans,
                                     is_velocity, fluxes_are_area_weighted,
                                     advection_string);


            // ************************************************************************
            // Density
            // ************************************************************************
            const int use_species_advection = fluid.isMixture() && fluid.solve_species();

            if (fluid.solve_density() && !use_species_advection)
            {

                face_comp = AMREX_SPACEDIM;
                ncomp = 1;
                is_velocity = false;
                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi,
                                                        ro_g_in[lev]->const_array(mfi),
                                          AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                       flux_y[lev].array(mfi,face_comp),
                                                       flux_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                       face_y[lev].array(mfi,face_comp),
                                                       face_z[lev].array(mfi,face_comp)),
                                          knownFaceStates,
                                          AMREX_D_DECL(ep_u_mac[lev]->const_array(mfi),
                                                       ep_v_mac[lev]->const_array(mfi),
                                                       ep_w_mac[lev]->const_array(mfi)),
                                          divu_arr, Array4<Real const>{},
                                          geom[lev], l_dt,
                                          get_density_bcrec(),
                                          get_density_bcrec_device_ptr(),
                                          get_density_iconserv_device_ptr(),
                                          ebfact,
                                          (m_embedded_boundaries.has_flow()) ? eb_scalars[lev]->const_array(mfi,0) : Array4<Real const>{},
                                          m_godunov_ppm, m_godunov_use_forces_in_trans,
                                          is_velocity, fluxes_are_area_weighted,
                                          advection_string);
            }

            // ************************************************************************
            // (Rho*Enthalpy)
            // ************************************************************************
            // Make a FAB holding (rho * enthalpy) that is the same size as the original enthalpy FAB
            FArrayBox rhohfab;
            if (fluid.solve_enthalpy())
            {
                Box rhoh_box = Box((*h_g_in[lev])[mfi].box());
                Array4<Real> rhoh;
                Array4<Real const>   h =  h_g_in[lev]->const_array(mfi);
                Array4<Real const> rho = ro_g_in[lev]->const_array(mfi);
                rhohfab.resize(rhoh_box, 1, The_Async_Arena());
                rhoh   = rhohfab.array();
                amrex::ParallelFor(rhoh_box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rhoh(i,j,k) = rho(i,j,k) * h(i,j,k);
                });

                face_comp = 4;
                ncomp = 1;
                is_velocity = false;

                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi, rhoh,
                                          AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                       flux_y[lev].array(mfi,face_comp),
                                                       flux_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                       face_y[lev].array(mfi,face_comp),
                                                       face_z[lev].array(mfi,face_comp)),
                                          knownFaceStates,
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
                                          (m_embedded_boundaries.has_flow()) ? eb_scalars[lev]->const_array(mfi,1) : Array4<Real const>{},
                                          m_godunov_ppm, m_godunov_use_forces_in_trans,
                                          is_velocity, fluxes_are_area_weighted,
                                          advection_string);
            }

            // ************************************************************************
            // (Rho*Tracer)
            // ************************************************************************
            // Make a FAB holding (rho * tracer) that is the same size as the original tracer FAB
            FArrayBox rhotracfab;
            if (fluid.solve_tracer() && (ntrac>0)) {

                Box rhotrac_box = Box((*trac_in[lev])[mfi].box());
                Array4<Real const> tra = trac_in[lev]->const_array(mfi);
                Array4<Real const> rho = ro_g_in[lev]->const_array(mfi);
                rhotracfab.resize(rhotrac_box, ntrac, The_Async_Arena());
                Array4<Real> rhotrac = rhotracfab.array();
                amrex::ParallelFor(rhotrac_box, ntrac,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhotrac(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
                });

                face_comp = 5;
                ncomp = ntrac;
                is_velocity = false;

                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi, rhotrac,
                                          AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                       flux_y[lev].array(mfi,face_comp),
                                                       flux_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                       face_y[lev].array(mfi,face_comp),
                                                       face_z[lev].array(mfi,face_comp)),
                                          knownFaceStates,
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
                                          (m_embedded_boundaries.has_flow()) ? eb_scalars[lev]->const_array(mfi,2) : Array4<Real const>{},
                                          m_godunov_ppm, m_godunov_use_forces_in_trans,
                                          is_velocity, fluxes_are_area_weighted,
                                          advection_string);
            }

            // ************************************************************************
            // (Rho*Species)
            // ************************************************************************
            // Make a FAB holding (rho * species) that is the same size as the original tracer FAB
            FArrayBox rhoXfab;
            if (fluid.solve_species() && (fluid.nspecies()>0))
            {
                Box rhoX_box = Box((*X_gk_in[lev])[mfi].box());
                Array4<Real const>   X = X_gk_in[lev]->const_array(mfi);
                Array4<Real const> rho = ro_g_in[lev]->const_array(mfi);
                rhoXfab.resize(rhoX_box, fluid.nspecies(), The_Async_Arena());
                Array4<Real> rhoX = rhoXfab.array();
                amrex::ParallelFor(rhoX_box, fluid.nspecies(),
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhoX(i,j,k,n) = rho(i,j,k) * X(i,j,k,n);
                });

                face_comp = 5+ntrac;
                ncomp = fluid.nspecies();
                is_velocity = false;

                HydroUtils::ComputeFluxesOnBoxFromState( bx, ncomp, mfi, rhoX,
                                          AMREX_D_DECL(flux_x[lev].array(mfi,face_comp),
                                                       flux_y[lev].array(mfi,face_comp),
                                                       flux_z[lev].array(mfi,face_comp)),
                                          AMREX_D_DECL(face_x[lev].array(mfi,face_comp),
                                                       face_y[lev].array(mfi,face_comp),
                                                       face_z[lev].array(mfi,face_comp)),
                                          knownFaceStates,
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
                                          (m_embedded_boundaries.has_flow()) ? eb_species[lev]->const_array(mfi) : Array4<Real const>{},
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
        if (fluid.solve_species() && fluid.nspecies() > 0)
            dXdt_tmp.define(vel_in[lev]->boxArray(),dmap[lev], fluid.nspecies() ,ngrow,MFInfo(),EBFactory(lev));

        // Must initialize to zero because not all values may be set, e.g. outside the domain.
        dvdt_tmp.setVal(0.);
        dsdt_tmp.setVal(0.);
        if (fluid.solve_species() && fluid.nspecies() > 0)
            dXdt_tmp.setVal(0.);

        const auto& ebfact = EBFactory(lev);

#ifdef _OPENMP
#pragma omp parallel private (flux_comp, num_comp) if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*conv_u[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();

            Real mult = -1.0;

            EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
            bool regular_bxg2 = (flagfab.getType(amrex::grow(bx,2)) == FabType::regular);

            auto const& vfrac = ebfact.getVolFrac().const_array(mfi);

            if (flagfab.getType(bx) != FabType::covered)
            {
                if (!regular_bxg2)
                {
                    flux_comp = 0;
                    num_comp = 3;
                    HydroUtils::EB_ComputeDivergence(bx, dvdt_tmp.array(mfi),
                                                     AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                                  flux_y[lev].const_array(mfi,flux_comp),
                                                                  flux_z[lev].const_array(mfi,flux_comp)),
                                                     vfrac, num_comp, geom[lev], mult, fluxes_are_area_weighted,
                                                     (m_embedded_boundaries.has_flow()) ? eb_vel[lev]->const_array(mfi) : Array4<Real const>{},
                                                     (m_embedded_boundaries.has_flow()) ? eb_vel[lev]->const_array(mfi) : Array4<Real const>{},
                                                     flagfab.const_array(),
                                                     ebfact.getBndryArea().const_array(mfi),
                                                     ebfact.getBndryNormal().const_array(mfi));


                    flux_comp = 3;
                    num_comp = 2+ntrac;
                    HydroUtils::EB_ComputeDivergence(bx, dsdt_tmp.array(mfi),
                                                     AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                                  flux_y[lev].const_array(mfi,flux_comp),
                                                                  flux_z[lev].const_array(mfi,flux_comp)),
                                                     vfrac, num_comp, geom[lev], mult, fluxes_are_area_weighted,
                                                     (m_embedded_boundaries.has_flow()) ? eb_vel[lev]->const_array(mfi) : Array4<Real const>{},
                                                     (m_embedded_boundaries.has_flow()) ? eb_scalars[lev]->const_array(mfi) : Array4<Real const>{},
                                                     flagfab.const_array(),
                                                     ebfact.getBndryArea().const_array(mfi),
                                                     ebfact.getBndryNormal().const_array(mfi));

                    flux_comp = 5+ntrac;
                    num_comp = fluid.nspecies();
                    if (fluid.solve_species() && fluid.nspecies() > 0)
                    {
                        HydroUtils::EB_ComputeDivergence(bx, dXdt_tmp.array(mfi),
                                                         AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                                      flux_y[lev].const_array(mfi,flux_comp),
                                                                      flux_z[lev].const_array(mfi,flux_comp)),
                                                         vfrac, num_comp, geom[lev], mult, fluxes_are_area_weighted,
                                                        (m_embedded_boundaries.has_flow()) ? eb_vel[lev]->const_array(mfi) : Array4<Real const>{},
                                                        (m_embedded_boundaries.has_flow()) ? eb_species[lev]->const_array(mfi) : Array4<Real const>{},
                                                        flagfab.const_array(),
                                                        ebfact.getBndryArea().const_array(mfi),
                                                        ebfact.getBndryNormal().const_array(mfi));
                    }
                }
                else
                {
                    // Velocity
                    flux_comp = 0;
                    num_comp = 3;
                    HydroUtils::ComputeDivergence(bx, dvdt_tmp.array(mfi),
                                              AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                           flux_y[lev].const_array(mfi,flux_comp),
                                                           flux_z[lev].const_array(mfi,flux_comp)),
                                              num_comp, geom[lev], mult,
                                              fluxes_are_area_weighted);

                    // Density, enthalpy, tracers
                    // NOTE: this isn't fully general because we pass in the density "iconserv"
                    //       flag but in this case density, (rho h) and (rho T) are all viewed as conservative
                    flux_comp = 3;
                    num_comp = 2 + ntrac;
                    HydroUtils::ComputeDivergence(bx, dsdt_tmp.array(mfi),
                                              AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                           flux_y[lev].const_array(mfi,flux_comp),
                                                           flux_z[lev].const_array(mfi,flux_comp)),
                                              num_comp, geom[lev], mult,
                                              fluxes_are_area_weighted);

                    // Species
                    flux_comp = 5+ntrac;
                     num_comp = fluid.nspecies();
                    if (fluid.solve_species() && fluid.nspecies() > 0)
                    {
                        HydroUtils::ComputeDivergence(bx, dXdt_tmp.array(mfi),
                                                  AMREX_D_DECL(flux_x[lev].const_array(mfi,flux_comp),
                                                               flux_y[lev].const_array(mfi,flux_comp),
                                                               flux_z[lev].const_array(mfi,flux_comp)),
                                                  num_comp, geom[lev], mult,
                                                  fluxes_are_area_weighted);
                    } // species
                } // regular

                // ***************************************************************
                // If convective, we define u dot grad u = div (u u) - u div(u)
                // ***************************************************************
                auto const& q           = vel_in[lev]->const_array(mfi,0);
                auto const& divu_arr    =   divu[lev].const_array(mfi);
                auto        upd_arr     = dvdt_tmp.array(mfi);
                AMREX_D_TERM(auto const& q_on_face_x  = face_x[lev].const_array(mfi);,
                             auto const& q_on_face_y  = face_y[lev].const_array(mfi);,
                             auto const& q_on_face_z  = face_z[lev].const_array(mfi););

                int const* iconserv_ptr = get_velocity_iconserv_device_ptr();
                if (advection_type() == AdvectionType::MOL)
                {
                    // Here we use q at the same time as the velocity
                    amrex::ParallelFor(bx, AMREX_SPACEDIM, [=]
                    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (!iconserv_ptr[n])
                        {
                            // Note that because we define upd_arr as MINUS div(u u), here we add u div (u)
                            upd_arr(i,j,k,n) += q(i,j,k,n)*divu_arr(i,j,k);
                        }
                    });

                } else if (advection_type() == AdvectionType::Godunov) {

                    // Here we want to use q predicted to t^{n+1/2}
                    //bool regular = (flagfab.getType(bx) == FabType::regular);
                    if (flagfab.getType(bx) == FabType::regular)
                    {
                        amrex::ParallelFor(bx, AMREX_SPACEDIM, [=]
                        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                        {
                            if (!iconserv_ptr[n])
                            {
                                Real qavg  = (q_on_face_x(i,j,k,n) + q_on_face_x(i+1,j,k,n) +
                                              q_on_face_y(i,j,k,n) + q_on_face_y(i,j+1,k,n) +
                                              q_on_face_z(i,j,k,n) + q_on_face_z(i,j,k+1,n) ) / 6.;

                                // Note that because we define upd_arr as MINUS div(u u), here we add u div (u)
                                upd_arr(i,j,k,n) += qavg*divu_arr(i,j,k);
                            }
                        });
                    }
                    else
                    {
                        auto const& apx_arr = ebfact.getAreaFrac()[0]->const_array(mfi);
                        auto const& apy_arr = ebfact.getAreaFrac()[1]->const_array(mfi);
                        auto const& apz_arr = ebfact.getAreaFrac()[2]->const_array(mfi);

                        amrex::ParallelFor(bx, AMREX_SPACEDIM, [=]
                        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                        {
                            if (!iconserv_ptr[n] && vfrac(i,j,k) > 0.)
                            {
                                Real qavg  = ( apx_arr(i,j,k)*q_on_face_x(i,j,k,n) + apx_arr(i+1,j,k)*q_on_face_x(i+1,j,k,n) +
                                               apy_arr(i,j,k)*q_on_face_y(i,j,k,n) + apy_arr(i,j+1,k)*q_on_face_y(i,j+1,k,n) +
                                               apz_arr(i,j,k)*q_on_face_z(i,j,k,n) + apz_arr(i,j,k+1)*q_on_face_z(i,j,k+1,n) ) /
                                             ( apx_arr(i,j,k) + apx_arr(i+1,j,k) + apy_arr(i,j,k) + apy_arr(i,j+1,k)
                                              +apz_arr(i,j,k) + apz_arr(i,j,k+1) );

                                // Note that because we define upd_arr as MINUS div(u u), here we add u div (u)
                                upd_arr(i,j,k,n) += qavg*divu_arr(i,j,k);
                            }
                        });
                    } // regular?
                } // Godunov
            } // not covered
        } // mfi

        // We only filled these on the valid cells so we fill same-level interior ghost cells here.
        // (We don't need values outside the domain or at a coarser level so we can call just FillBoundary)
        dvdt_tmp.FillBoundary(geom[lev].periodicity());
        dsdt_tmp.FillBoundary(geom[lev].periodicity());
        if (fluid.solve_species() && fluid.nspecies() > 0)
            dXdt_tmp.FillBoundary(geom[lev].periodicity());

        int max_ncomp = std::max(std::max(AMREX_SPACEDIM, fluid.nspecies()), 2+ntrac);

        for (MFIter mfi(*ro_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.tilebox();

          EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
          bool regular = (flagfab.getType(amrex::grow(bx,4)) == FabType::regular);

          if (flagfab.getType(bx) != FabType::covered)
          {
           if (!regular)
           {
            // Make a FAB holding (rho * tracer) that is the same size as the original tracer FAB
            Box tmp_box = Box((*ro_g_in[lev])[mfi].box());
            Array4<Real const> rho = ro_g_in[lev]->const_array(mfi);

            FArrayBox scratchfab(amrex::grow(bx,ngrow), max_ncomp, The_Async_Arena());

            auto const& vfrac = ebfact.getVolFrac().const_array(mfi);
            auto const& ccc   = ebfact.getCentroid().const_array(mfi);

            auto const& apx = ebfact.getAreaFrac()[0]->const_array(mfi);
            auto const& apy = ebfact.getAreaFrac()[1]->const_array(mfi);
            auto const& apz = ebfact.getAreaFrac()[2]->const_array(mfi);

            auto const& fcx = ebfact.getFaceCent()[0]->const_array(mfi);
            auto const& fcy = ebfact.getFaceCent()[1]->const_array(mfi);
            auto const& fcz = ebfact.getFaceCent()[2]->const_array(mfi);

            Array4<const Real> epg = ep_g_in[lev]->array(mfi);

            // Velocity
            int ncomp = 3;

            auto const& bc_vel = get_hydro_velocity_bcrec_device_ptr();

            ApplyRedistribution(bx, ncomp, conv_u[lev]->array(mfi), dvdt_tmp.array(mfi),
                                vel_in[lev]->const_array(mfi),
                                (m_redistribution_type == "StateRedist") ? scratchfab.array()
                                                                         : ep_g_in[lev]->array(mfi),
                                flagfab.const_array(),
                                apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                bc_vel, geom[lev], l_dt,
                                m_redistribution_type, true, 2, Real(0.5), epg);

            // Density
            const int use_species_advection = fluid.isMixture() && fluid.solve_species();

            if (fluid.solve_density() && !use_species_advection)
            {
                ncomp = 1;

                auto const& bc_den = get_density_bcrec_device_ptr();

                ApplyRedistribution(bx, ncomp, conv_s[lev]->array(mfi,0), dsdt_tmp.array(mfi,0),
                                    ro_g_in[lev]->const_array(mfi),
                                    (m_redistribution_type == "StateRedist") ? scratchfab.array()
                                                                             : ep_g_in[lev]->array(mfi),
                                    flagfab.const_array(),
                                    apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                    bc_den, geom[lev], l_dt,
                                    m_redistribution_type, true, 2, Real(0.5), epg);
            }

            // Enthalpy
            if (fluid.solve_enthalpy())
            {
                ncomp = 1;
                Array4<Real const>   h =  h_g_in[lev]->const_array(mfi);
                FArrayBox rhohfab;
                rhohfab.resize(tmp_box, 1, The_Async_Arena());
                Array4<Real> rhoh = rhohfab.array();
                amrex::ParallelFor(tmp_box, ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhoh(i,j,k,n) = rho(i,j,k) * h(i,j,k,n);
                });

                auto const& bc_rh = get_enthalpy_bcrec_device_ptr();

                ApplyRedistribution(bx, ncomp, conv_s[lev]->array(mfi,1), dsdt_tmp.array(mfi,1), rhoh,
                                    (m_redistribution_type == "StateRedist") ? scratchfab.array()
                                                                             : ep_g_in[lev]->array(mfi),
                                    flagfab.const_array(),
                                    apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                    bc_rh, geom[lev], l_dt,
                                    m_redistribution_type, true, 2, Real(0.5), epg);
            }

            // Tracers
            if (fluid.solve_tracer() && (ntrac > 0)) {
                ncomp = ntrac;

                Array4<Real const> tra = trac_in[lev]->const_array(mfi);
                FArrayBox rhotracfab;
                rhotracfab.resize(tmp_box, ntrac, The_Async_Arena());
                Array4<Real> rhotrac = rhotracfab.array();
                amrex::ParallelFor(tmp_box, ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhotrac(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
                });

                auto const& bc_rt = get_tracer_bcrec_device_ptr();

                ApplyRedistribution(bx, ncomp, conv_s[lev]->array(mfi,2), dsdt_tmp.array(mfi,2), rhotrac,
                                    (m_redistribution_type == "StateRedist") ? scratchfab.array()
                                                                             : ep_g_in[lev]->array(mfi),
                                    flagfab.const_array(),
                                    apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                    bc_rt, geom[lev], l_dt,
                                    m_redistribution_type, true, 2, Real(0.5), epg);
            }

            if (fluid.solve_species() && (fluid.nspecies() > 0))
            {

                ncomp = fluid.nspecies();

                Array4<Real const>   X = X_gk_in[lev]->const_array(mfi);
                FArrayBox rhoXfab;
                rhoXfab.resize(tmp_box, fluid.nspecies(), The_Async_Arena());
                Array4<Real> rhoX = rhoXfab.array();
                amrex::ParallelFor(tmp_box, ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhoX(i,j,k,n) = rho(i,j,k) * X(i,j,k,n);
                });

                auto const& bc_rX = get_species_bcrec_device_ptr();

                ApplyRedistribution(bx, ncomp, conv_X[lev]->array(mfi), dXdt_tmp.array(mfi), rhoX,
                                    (m_redistribution_type == "StateRedist") ? scratchfab.array()
                                                                             : ep_g_in[lev]->array(mfi),
                                    flagfab.const_array(),
                                    apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                                    bc_rX, geom[lev], l_dt,
                                    m_redistribution_type, true, 2, Real(0.5), epg);

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
             ncomp = conv_s[lev]->nComp();
             amrex::ParallelFor(bx,ncomp,
             [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
             {
                 conv_s_arr(i,j,k,n) = dsdt_arr(i,j,k,n);
             });

             if (fluid.solve_species() && (fluid.nspecies() > 0))
             {
                 Array4<Real      > conv_X_arr = conv_X[lev]->array(mfi);
                 Array4<Real const>   dXdt_arr = dXdt_tmp.array(mfi);
                 ncomp = conv_X[lev]->nComp();
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


  if (m_rw->report_mass_balance) {

    flux_comp = 5+ntrac;
    num_comp = fluid.nspecies();

    m_rw->ComputeMassFlux(GetVecOfConstPtrs(flux_x),
                            GetVecOfConstPtrs(flux_y),
                            GetVecOfConstPtrs(flux_z),
                            flux_comp, num_comp, fluxes_are_area_weighted,
                            m_embedded_boundaries.has_flow(), eb_vel, eb_species,
                            l_dt);
  }

}
