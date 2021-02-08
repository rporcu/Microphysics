#include <MOL.H>
#include <Godunov.H>

#include <mfix_bc_parms.H>
#include <mfix_mf_helpers.H>
#include <mfix.H>

#include <AMReX_MacProjector.H>

//
// Computes the following decomposition:
//
//    u* = u + c*grad(phi)/ro with  div(ep u) = 0
//
// Notes:
//
//  phi is computed by solving
//
//       div(ep*grad(phi)/ro) = div(ep u*) - S
//
//  This method returns the MAC velocity with up-to-date BCs in place
//
void
mfix::compute_MAC_projected_velocities (amrex::Real time, const amrex::Real l_dt,
                                        amrex::Vector< amrex::MultiFab const*> const& vel_in,
                                        amrex::Vector< amrex::MultiFab*      > const& ep_u_mac,
                                        amrex::Vector< amrex::MultiFab*      > const& ep_v_mac,
                                        amrex::Vector< amrex::MultiFab*      > const& ep_w_mac,
                                        amrex::Vector< amrex::MultiFab const*> const& ep_g_in,
                                        amrex::Vector< amrex::MultiFab const*> const& ro_g_in,
                                        amrex::Vector< amrex::MultiFab      *> const& vel_forces,
                                        amrex::Vector< amrex::MultiFab const*> const& rhs_mac_in)
{
  BL_PROFILE("mfix::compute_MAC_projected_velocities()");

  if (m_verbose) {
    Print() << "MAC Projection:\n";
  }

  auto mac_phi = get_mac_phi();

  // ro_face and ep_face are temporary, no need to keep it outside this routine
  Vector< Array<MultiFab*,3> > ro_face(finest_level+1);
  Vector< Array<MultiFab*,3> > ep_face(finest_level+1);

  for ( int lev=0; lev <= finest_level; ++lev )
  {
    // ep_g_in[lev]->FillBoundary(geom[lev].periodicity());
    // ro_g_in[lev]->FillBoundary(geom[lev].periodicity());

    ep_face[lev][0] = new MultiFab(ep_u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]);
    ep_face[lev][1] = new MultiFab(ep_v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]);
    ep_face[lev][2] = new MultiFab(ep_w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]);

    ep_face[lev][0]->setVal(covered_val);
    ep_face[lev][1]->setVal(covered_val);
    ep_face[lev][2]->setVal(covered_val);

    ro_face[lev][0] = new MultiFab(ep_u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]);
    ro_face[lev][1] = new MultiFab(ep_v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]);
    ro_face[lev][2] = new MultiFab(ep_w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]);

    // Define ep and rho on face centroids (using interpolation from cell centroids)
    // The only use of bcs in this call is to test on whether a domain boundary is ext_dir
    // average_cellcenter_to_face(ro_face[lev], *ro_g_in[lev], geom[lev]);
    // average_cellcenter_to_face(ep_face[lev], *ep_g_in[lev], geom[lev]);
    EB_interp_CellCentroid_to_FaceCentroid (*ro_g_in[lev], ro_face[lev], 0, 0, 1, geom[lev], bcs_s);
    EB_interp_CellCentroid_to_FaceCentroid (*ep_g_in[lev], ep_face[lev], 0, 0, 1, geom[lev], bcs_s);

    // These will be reused to predict velocites (ep*u) on faces
    ep_face[lev][0]->FillBoundary();
    ep_face[lev][1]->FillBoundary();
    ep_face[lev][2]->FillBoundary();

    // Compute ep_face into bcoeff
    MultiFab::Copy(*bcoeff[lev][0], *(ep_face[lev][0]), 0, 0, 1, 0);
    MultiFab::Copy(*bcoeff[lev][1], *(ep_face[lev][1]), 0, 0, 1, 0);
    MultiFab::Copy(*bcoeff[lev][2], *(ep_face[lev][2]), 0, 0, 1, 0);

    // Compute beta coefficients for div(beta*grad(phi)) = RHS:  beta = ep / ro
    MultiFab::Divide(*bcoeff[lev][0], *(ro_face[lev][0]), 0, 0, 1, 0);
    MultiFab::Divide(*bcoeff[lev][1], *(ro_face[lev][1]), 0, 0, 1, 0);
    MultiFab::Divide(*bcoeff[lev][2], *(ro_face[lev][2]), 0, 0, 1, 0);
  }


  Vector< Array <MultiFab const*, 3>> const_bcoeff;
  const_bcoeff.reserve(bcoeff.size());
  for (const auto& x : bcoeff) const_bcoeff.push_back(GetArrOfConstPtrs(x));

  //
  // Initialize (or redefine the beta in) the MacProjector
  //

  if (macproj->needInitialization())
    {
      LPInfo lp_info;
      // If we want to set max_coarsening_level we have to send it in to the constructor
      lp_info.setMaxCoarseningLevel(mac_mg_max_coarsening_level);
      macproj->initProjector(lp_info, const_bcoeff);
      macproj->setDomainBC(BC::ppe_lobc, BC::ppe_hibc);
    } else {
    macproj->updateBeta(const_bcoeff);
  }


  Vector<Array<MultiFab,AMREX_SPACEDIM> > m_fluxes;
  m_fluxes.resize(finest_level+1);
  for (int lev=0; lev <= finest_level; ++lev) {
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      m_fluxes[lev][idim].define(
         amrex::convert(grids[lev], IntVect::TheDimensionVector(idim)),
         dmap[lev], 1, 0, MFInfo(), EBFactory(lev));
    }
  }

  if (m_use_mac_phi_in_godunov) {
    macproj->getFluxes(amrex::GetVecOfArrOfPtrs(m_fluxes), mac_phi, MLMG::Location::FaceCentroid);
  } else {
    for (int lev=0; lev <= finest_level; ++lev)
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        m_fluxes[lev][idim].setVal(0.);
  }



  // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
  //    arrays returned from this call are on face CENTROIDS and have been
  //    multiplied by the phasic voluem fraction {ep * u_mac, ep * v_mac, ep * w_mac}
  for (int lev = 0; lev < nlev; ++lev) {

    mac_phi[lev]->FillBoundary(geom[lev].periodicity());

    const EBFArrayBoxFactory* ebfact = &EBFactory(lev);

    // We need this to avoid FPE
    ep_u_mac[lev]->setVal(covered_val);
    ep_v_mac[lev]->setVal(covered_val);
    ep_w_mac[lev]->setVal(covered_val);

    if (advection_type() == AdvectionType::Godunov) {


      if (ebfact->isAllRegular()) {
#if 0
        godunov::predict_godunov(time,
                                 *ep_u_mac[lev], *ep_v_mac[lev], *ep_w_mac[lev],
                                 *vel[lev], *vel_forces[lev],
                                 get_velocity_bcrec(), get_velocity_bcrec_device_ptr(),
                                 geom[lev], l_dt, m_godunov_ppm, m_godunov_use_forces_in_trans,
                                 m_fluxes[lev][0], m_fluxes[lev][1], m_fluxes[lev][2],
                                 m_use_mac_phi_in_godunov);
#endif
      } else {
#if 0
        ebgodunov::predict_godunov(time,
                                   *ep_u_mac[lev], *ep_v_mac[lev], *ep_w_mac[lev],
                                   *vel_in[lev], *vel_forces[lev],
                                   get_velocity_bcrec(), get_velocity_bcrec_device_ptr(),
                                   ebfact, geom[lev], l_dt,
                                   m_fluxes[lev][0], m_fluxes[lev][1], m_fluxes[lev][2],
                                   m_use_mac_phi_in_godunov);
#endif

      }


    } else if (advection_type() == AdvectionType::MOL) {

      mol::predict_vels_on_faces(lev,
                                 *ep_u_mac[lev],   *ep_v_mac[lev],   *ep_w_mac[lev],
                                 *ep_face[lev][0], *ep_face[lev][1], *ep_face[lev][2],
                                 *vel_in[lev], m_vel_g_bc_types,
                                 bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                 bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                 bc_klo[lev]->array(), bc_khi[lev]->array(),
                                 ebfact, geom);
    }
  }


  // Check that everything is consistent with amrcore
  // update_internals();

  // Setup for solve
  Vector<Array<MultiFab*,AMREX_SPACEDIM> > mac_vec(finest_level+1);

  // We only need this MF if were are verbose.
  Vector<MultiFab> tmp_div;

  if (m_verbose) {
    Print() << " >> Before projection\n" ;
    for (int lev = 0; lev <= finest_level; ++lev) {
      tmp_div.emplace_back(grids[lev], dmap[lev], 1, mfix::nghost,
                           MFInfo(), EBFactory(lev));
    }
  }

  for ( int lev=0; lev <= finest_level; ++lev )
  {

    // Store (ep * u) in temporaries
    (mac_vec[lev])[0] = ep_u_mac[lev];
    (mac_vec[lev])[1] = ep_v_mac[lev];
    (mac_vec[lev])[2] = ep_w_mac[lev];

    for (int i=0; i < 3; ++i)
      (mac_vec[lev])[i]->FillBoundary(geom[lev].periodicity());

    if (m_verbose)
    {
      bool already_on_centroid = true;
      EB_computeDivergence(tmp_div[lev], GetArrOfConstPtrs(mac_vec[lev]),
          geom[lev], already_on_centroid);

      Print() << "  * On level "<< lev << " max(abs(diveu)) = "
              << tmp_div[lev].norm0(0,0,false,true) << "\n";
    }
  }

  for ( int lev=0; lev <= finest_level; ++lev )
  {
    delete ep_face[lev][0];
    delete ep_face[lev][1];
    delete ep_face[lev][2];

    delete ro_face[lev][0];
    delete ro_face[lev][1];
    delete ro_face[lev][2];
  }

  macproj->setUMAC(mac_vec);
  macproj->setDivU(rhs_mac_in);

  if (m_steady_state)
  {
    // Solve using mac_phi as an initial guess -- note that mac_phi is
    //       stored from iteration to iteration
    macproj->project(mac_phi,mac_mg_rtol, mac_mg_atol);

  }
  else
  {
    // Solve with initial guess of zero
    macproj->project(mac_mg_rtol, mac_mg_atol);
  }

  // Get MAC velocities at face CENTER by dividing solution by ep at faces
  if (m_verbose)
    Print() << " >> After projection\n" ;

  for ( int lev=0; lev <= finest_level ; ++lev )
  {
    if (m_verbose)
    {
      mac_vec[lev][0]->FillBoundary(geom[lev].periodicity());
      mac_vec[lev][1]->FillBoundary(geom[lev].periodicity());
      mac_vec[lev][2]->FillBoundary(geom[lev].periodicity());

      bool already_on_centroid = true;
      EB_computeDivergence(tmp_div[lev], GetArrOfConstPtrs(mac_vec[lev]),
                           geom[lev], already_on_centroid);

      Print() << "  * On level "<< lev << " max(abs(diveu)) = "
              << tmp_div[lev].norm0(0,0,false,true) << "\n";
    }

    // Set bcs on (ep * u_mac)
    set_MAC_velocity_bcs(lev, rhs_mac_in, ep_u_mac, ep_v_mac, ep_w_mac, time);

    ep_u_mac[lev]->FillBoundary(geom[lev].periodicity());
    ep_v_mac[lev]->FillBoundary(geom[lev].periodicity());
    ep_w_mac[lev]->FillBoundary(geom[lev].periodicity());
  }
}
