#include <hydro_utils.H>
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
mfix::compute_MAC_projected_velocities (Real time, const amrex::Real l_dt,
                                        Vector< MultiFab const*> const& vel_in,
                                        Vector< MultiFab*      > const& ep_u_mac,
                                        Vector< MultiFab*      > const& ep_v_mac,
                                        Vector< MultiFab*      > const& ep_w_mac,
                                        Vector< MultiFab const*> const& ep_g_in,
                                        Vector< MultiFab const*> const& ro_g_in,
                                        Vector< MultiFab const*> const& txfr_in,
                                        Vector< MultiFab      *> const& vel_forces,
                                        Vector< MultiFab const*> const& rhs_mac_in)
{
  BL_PROFILE("mfix::compute_MAC_projected_velocities()");

  if (m_verbose) {
    Print() << "MAC Projection:\n";
  }

  // We first compute the velocity forcing terms to be used in predicting
  //    to faces before the MAC projection
  if (advection_type() != AdvectionType::MOL)
  {
    bool include_pressure_gradient = !(m_use_mac_phi_in_godunov);
    bool include_drag_force = include_pressure_gradient && m_use_drag_in_godunov;
    compute_vel_forces(vel_forces, vel_in, ro_g_in, txfr_in, include_pressure_gradient, include_drag_force);

    if (m_godunov_include_diff_in_forcing)
      for (int lev = 0; lev <= finest_level; ++lev)
        MultiFab::Add(*vel_forces[lev], *m_leveldata[lev]->divtau_o, 0, 0, 3, 0);

    if (nghost_force() > 0)
      fillpatch_force(time, vel_forces, nghost_force());
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
    EB_interp_CellCentroid_to_FaceCentroid (*ro_g_in[lev], ro_face[lev], 0, 0, 1, geom[lev], get_density_bcrec());
    EB_interp_CellCentroid_to_FaceCentroid (*ep_g_in[lev], ep_face[lev], 0, 0, 1, geom[lev], bcs_f);

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

    // Copy mac_phi MultiFabs into temporary variables
    Vector<MultiFab*> mac_phi_copy(finest_level+1);
    for (int lev(0); lev <= finest_level; ++lev) {
      mac_phi_copy[lev] = (MFHelpers::createFrom(*mac_phi[lev])).release();
    }

    macproj->getFluxes(amrex::GetVecOfArrOfPtrs(m_fluxes), mac_phi_copy, MLMG::Location::FaceCentroid);

    for (int lev(0); lev <= finest_level; ++lev) {
      delete mac_phi_copy[lev];
    }

    if( m_use_drag_in_godunov) {
      amrex::Print() << "WARNING: using mac_phi in Godunov MAC velocities.\n"
                     << "This should probably include the explicit drag force too.\n";
    }

  } else {
    for (int lev=0; lev <= finest_level; ++lev)
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        m_fluxes[lev][idim].setVal(0.);
  }



  // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
  //    arrays returned from this call are on face CENTROIDS.
  for (int lev = 0; lev < nlev; ++lev) {

    mac_phi[lev]->FillBoundary(geom[lev].periodicity());

    const EBFArrayBoxFactory* ebfact = &EBFactory(lev);

    // We need this to avoid FPE
    ep_u_mac[lev]->setVal(covered_val);
    ep_v_mac[lev]->setVal(covered_val);
    ep_w_mac[lev]->setVal(covered_val);

    std::string advection_string;
    if (advection_type() == AdvectionType::Godunov) 
        advection_string = "Godunov";
    else
        advection_string = "MOL";

    HydroUtils::ExtrapVelToFaces(*vel_in[lev], *vel_forces[lev],
                                 *ep_u_mac[lev], *ep_v_mac[lev], *ep_w_mac[lev],
                                  get_hydro_velocity_bcrec(), get_hydro_velocity_bcrec_device_ptr(),
                                  geom[lev], l_dt, *ebfact,
                                  m_godunov_ppm, m_godunov_use_forces_in_trans,
                                  advection_string);

    for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& ubx = mfi.nodaltilebox(0);
      const Box& vbx = mfi.nodaltilebox(1);
      const Box& wbx = mfi.nodaltilebox(2);

      // Face-centered velocity components
      Array4<Real      > const& umac_arr = ep_u_mac[lev]->array(mfi);
      Array4<Real      > const& vmac_arr = ep_v_mac[lev]->array(mfi);
      Array4<Real      > const& wmac_arr = ep_w_mac[lev]->array(mfi);

      // Face-centroid volume fractions
      Array4<Real const> const& epx_arr = ep_face[lev][0]->const_array(mfi);
      Array4<Real const> const& epy_arr = ep_face[lev][1]->const_array(mfi);
      Array4<Real const> const& epz_arr = ep_face[lev][2]->const_array(mfi);

      // Now we multiply the face velocities by the phasic volume fraction
      // so we have {ep * u_mac, ep * v_mac, ep * w_mac}.
      amrex::ParallelFor(ubx, vbx, wbx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { umac_arr(i,j,k) *= epx_arr(i,j,k); },
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { vmac_arr(i,j,k) *= epy_arr(i,j,k); },
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { wmac_arr(i,j,k) *= epz_arr(i,j,k); });
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
      tmp_div.emplace_back(grids[lev], dmap[lev], 1, nghost_state(), MFInfo(), EBFactory(lev));
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

    if (m_use_mac_phi_in_godunov) {
      for (int lev=0; lev <= finest_level; ++lev)
        mac_phi[lev]->mult(l_dt/2.,0,1,1);

      macproj->project(mac_phi,mac_mg_rtol,mac_mg_atol);

      for (int lev=0; lev <= finest_level; ++lev)
        mac_phi[lev]->mult(2./l_dt,0,1,1);

    } else {

      // Solve with initial guess of zero
      macproj->project(mac_mg_rtol, mac_mg_atol);
    }
  }

  // Get MAC velocities at face CENTER by dividing solution by ep at faces
  if (m_verbose)
    Print() << " >> After projection\n" ;

  for ( int lev=0; lev <= finest_level ; ++lev ) {

    if (m_verbose) {

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
