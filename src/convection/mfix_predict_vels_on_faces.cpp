#include <mfix.H>

void
mfix::mfix_predict_vels_on_faces (int lev, Real time,
                                  Vector< MultiFab* > const& vel_in,
                                  Vector< MultiFab* > const& ep_u_mac,
                                  Vector< MultiFab* > const& ep_v_mac,
                                  Vector< MultiFab* > const& ep_w_mac,
                                  Vector< MultiFab* > const& ep_in)

{
  BL_PROFILE("mfix::mfix_predict_vels_on_faces");

  Array< MultiFab*, 3> ep_face;

  iMultiFab cc_mask(grids[lev], dmap[lev], 1, 1);

  const int covered_value = 1;
  const int notcovered_value = 0;
  const int physical_boundaries_value = 0;
  const int interior_value = 1;

  cc_mask.BuildMask(geom[lev].Domain(), geom[lev].periodicity(),
                    covered_value, notcovered_value,
                    physical_boundaries_value, interior_value);

  // Get EB geometric info
  Array< const MultiCutFab*,3> areafrac;
  Array< const MultiCutFab*,3> facecent;

  areafrac = ebfactory[lev]->getAreaFrac();
  facecent = ebfactory[lev]->getFaceCent();
  const auto& cellcent = ebfactory[lev]->getCentroid();

  Real small_vel = 1.e-10;

  // ****************************************************************************
  // We will need ep on face centers to interpolate to face centroids below
  // ****************************************************************************

  ep_face[0] = new MultiFab(ep_u_mac[lev]->boxArray(),dmap[lev],1,1,MFInfo(),*ebfactory[lev]);
  ep_face[1] = new MultiFab(ep_v_mac[lev]->boxArray(),dmap[lev],1,1,MFInfo(),*ebfactory[lev]);
  ep_face[2] = new MultiFab(ep_w_mac[lev]->boxArray(),dmap[lev],1,1,MFInfo(),*ebfactory[lev]);

  // This is to make sure ep_face is defined everywhere
  ep_face[0]->setVal(covered_val);
  ep_face[1]->setVal(covered_val);
  ep_face[2]->setVal(covered_val);

  ep_in[lev]->FillBoundary(geom[lev].periodicity());
  // average_cellcenter_to_face(ep_face, *ep_in[lev], geom[lev]);
  // Vector<BCRec> bcs_s; // Just needed for this to compile
  EB_interp_CellCentroid_to_FaceCentroid (*ep_in[lev], ep_face, 0, 0, 1, geom[lev], bcs_s);

  ep_face[0]->FillBoundary();
  ep_face[1]->FillBoundary();
  ep_face[2]->FillBoundary();

  // ****************************************************************************
  // First compute the slopes
  // ****************************************************************************
  int slopes_comp = 0;

  mfix_compute_slopes(lev, time, *vel_in[lev],
                      get_xslopes_u(), get_yslopes_u(), get_zslopes_u(),
                      slopes_comp, m_vel_g_bc_types);

  // ****************************************************************************
  // Then predict to face centers
  // ****************************************************************************

  for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    // Tilebox
    const Box bx = mfi.tilebox();

    Box ubx = mfi.nodaltilebox(0);
    Box vbx = mfi.nodaltilebox(1);
    Box wbx = mfi.nodaltilebox(2);

    const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
    const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

    // Cell-centered velocity
    const auto& ccvel_fab = vel_in[lev]->array(mfi);

    // Cell-centered slopes
    const auto& xslopes_fab = m_leveldata[lev]->xslopes_u->array(mfi);
    const auto& yslopes_fab = m_leveldata[lev]->yslopes_u->array(mfi);
    const auto& zslopes_fab = m_leveldata[lev]->zslopes_u->array(mfi);

    // Face-centered velocity components
    const auto& umac_fab = (ep_u_mac[lev])->array(mfi);
    const auto& vmac_fab = (ep_v_mac[lev])->array(mfi);
    const auto& wmac_fab = (ep_w_mac[lev])->array(mfi);

    if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
    {
      Real val = 1.2345e300;

      (*ep_u_mac[lev])[mfi].setVal<RunOn::Device>(val, ubx, 0, 1);
      (*ep_v_mac[lev])[mfi].setVal<RunOn::Device>(val, vbx, 0, 1);
      (*ep_w_mac[lev])[mfi].setVal<RunOn::Device>(val, wbx, 0, 1);
    }

    // No cut cells in this FAB
    else if (flags.getType(amrex::grow(bx,2)) == FabType::regular )
    {
      // Face-centered ep
      const auto& epx_fab = (ep_face[0])->array(mfi);
      const auto& epy_fab = (ep_face[1])->array(mfi);
      const auto& epz_fab = (ep_face[2])->array(mfi);

      const int ubx_npoints = ubx.numPts();
      const auto ubx_lo = amrex::lbound(ubx);
      const auto ubx_len = amrex::length(ubx);

      const int vbx_npoints = vbx.numPts();
      const auto vbx_lo = amrex::lbound(vbx);
      const auto vbx_len = amrex::length(vbx);

      const int wbx_npoints = wbx.numPts();
      const auto wbx_lo = amrex::lbound(wbx);
      const auto wbx_len = amrex::length(wbx);

      const int npoints = amrex::max(ubx_npoints,vbx_npoints,wbx_npoints);

      ParallelFor(npoints, [=] AMREX_GPU_DEVICE (int idx) noexcept
      {
        if(idx < ubx_npoints)
        {
          int k = idx / (ubx_len.x*ubx_len.y);
          int j = (idx - k*(ubx_len.x*ubx_len.y)) / (ubx_len.x);
          int i = (idx - k*(ubx_len.x*ubx_len.y)) - j*ubx_len.x;

          i += ubx_lo.x;
          j += ubx_lo.y;
          k += ubx_lo.z;

          // X-faces
          const Real upls = ccvel_fab(i  ,j,k,0) - .5 * xslopes_fab(i  ,j,k,0);
          const Real umns = ccvel_fab(i-1,j,k,0) + .5 * xslopes_fab(i-1,j,k,0);
          Real umac(0);

          if (umns >= 0 or upls <= 0) {
            Real avg = .5 * ( upls + umns );

            if (avg >= small_vel) {
              umac = umns;
            }
            else if(avg <= -small_vel) {
              umac = upls;
            }

            umac *= epx_fab(i,j,k);
          }

          umac_fab(i,j,k) = umac;
        }

        if(idx < vbx_npoints)
        {
          int k = idx / (vbx_len.x*vbx_len.y);
          int j = (idx - k*(vbx_len.x*vbx_len.y)) / (vbx_len.x);
          int i = (idx - k*(vbx_len.x*vbx_len.y)) - j*vbx_len.x;

          i += vbx_lo.x;
          j += vbx_lo.y;
          k += vbx_lo.z;

          // Y-faces
          const Real vpls = ccvel_fab(i,j  ,k,1) - .5 * yslopes_fab(i,j  ,k,1);
          const Real vmns = ccvel_fab(i,j-1,k,1) + .5 * yslopes_fab(i,j-1,k,1);
          Real vmac(0);

          if (vmns >= 0 or vpls <= 0) {
            Real avg = .5 * (vpls + vmns);

            if (avg >= small_vel) {
              vmac = vmns;
            }
            else if (avg <= -small_vel) {
              vmac = vpls;
            }

            vmac *= epy_fab(i,j,k);
          }
          vmac_fab(i,j,k) = vmac;
        }

        if(idx < wbx_npoints)
        {
          int k = idx / (wbx_len.x*wbx_len.y);
          int j = (idx - k*(wbx_len.x*wbx_len.y)) / (wbx_len.x);
          int i = (idx - k*(wbx_len.x*wbx_len.y)) - j*wbx_len.x;

          i += wbx_lo.x;
          j += wbx_lo.y;
          k += wbx_lo.z;

          // Z-faces
          const Real wpls = ccvel_fab(i,j,k  ,2) - .5 * zslopes_fab(i,j,k  ,2);
          const Real wmns = ccvel_fab(i,j,k-1,2) + .5 * zslopes_fab(i,j,k-1,2);
          Real wmac(0);

          if (wmns >= 0 or wpls <= 0) {
            Real avg = .5 * (wpls + wmns);

            if (avg >= small_vel) {
              wmac = wmns;
            }
            else if (avg <= -small_vel) {
              wmac = wpls;
            }

            wmac *= epz_fab(i,j,k);
          }

          wmac_fab(i,j,k) = wmac;
        }
      });
    }
    // Cut cells in this FAB
    else
    {
      // Face centroids
      const auto& fcx_fab = facecent[0]->array(mfi);
      const auto& fcy_fab = facecent[1]->array(mfi);
      const auto& fcz_fab = facecent[2]->array(mfi);

      // Cell centroids
      const auto& ccc_fab = cellcent.array(mfi);

      // Face-centered ep
      const auto& epx_fab = (ep_face[0])->array(mfi);
      const auto& epy_fab = (ep_face[1])->array(mfi);
      const auto& epz_fab = (ep_face[2])->array(mfi);

      // Face-centered areas
      const auto& apx_fab = areafrac[0]->array(mfi);
      const auto& apy_fab = areafrac[1]->array(mfi);
      const auto& apz_fab = areafrac[2]->array(mfi);

      const int ubx_npoints = ubx.numPts();
      const auto ubx_lo = amrex::lbound(ubx);
      const auto ubx_len = amrex::length(ubx);

      const int vbx_npoints = vbx.numPts();
      const auto vbx_lo = amrex::lbound(vbx);
      const auto vbx_len = amrex::length(vbx);

      const int wbx_npoints = wbx.numPts();
      const auto wbx_lo = amrex::lbound(wbx);
      const auto wbx_len = amrex::length(wbx);

      const int npoints = amrex::max(ubx_npoints,vbx_npoints,wbx_npoints);

      // This FAB has cut cells -- we predict from cell centroids to face
      // centroids
      ParallelFor(npoints, [=] AMREX_GPU_DEVICE (int idx) noexcept
      {
        if(idx < ubx_npoints)
        {
          int k = idx / (ubx_len.x*ubx_len.y);
          int j = (idx - k*(ubx_len.x*ubx_len.y)) / (ubx_len.x);
          int i = (idx - k*(ubx_len.x*ubx_len.y)) - j*ubx_len.x;

          i += ubx_lo.x;
          j += ubx_lo.y;
          k += ubx_lo.z;

          // X-faces
          if (apx_fab(i,j,k) > 0.0)
          {
            // local (y,z) of centroid of x-face we are extrapolating to
            Real yf = fcx_fab(i,j,k,0);
            Real zf = fcx_fab(i,j,k,1);

            Real delta_x = .5 + ccc_fab(i,j,k,0);
            Real delta_y = yf - ccc_fab(i,j,k,1);
            Real delta_z = zf - ccc_fab(i,j,k,2);

            const Real ccvel_pls = ccvel_fab(i,j,k,0);
            const Real ccvel_mns = ccvel_fab(i-1,j,k,0);

            Real cc_umax = amrex::max(ccvel_pls, ccvel_mns);
            Real cc_umin = amrex::min(ccvel_pls, ccvel_mns);

            Real upls = ccvel_pls - delta_x * xslopes_fab(i,j,k,0)
                                  + delta_y * yslopes_fab(i,j,k,0)
                                  + delta_z * zslopes_fab(i,j,k,0);

            upls = amrex::min(upls, cc_umax);
            upls = amrex::max(upls, cc_umin);

            delta_x = .5 - ccc_fab(i-1,j,k,0);
            delta_y = yf - ccc_fab(i-1,j,k,1);
            delta_z = zf - ccc_fab(i-1,j,k,2);

            Real umns = ccvel_mns + delta_x * xslopes_fab(i-1,j,k,0)
                                  + delta_y * yslopes_fab(i-1,j,k,0)
                                  + delta_z * zslopes_fab(i-1,j,k,0);

            umns = amrex::min(umns, cc_umax);
            umns = amrex::max(umns, cc_umin);

            Real umac(0);

            if (umns >= 0 or upls <= 0) {
              Real avg = .5 * (upls + umns);

              if (avg >= small_vel) {
                umac = umns;
              }
              else if (avg <= -small_vel) {
                umac = upls;
              }

              umac *= epx_fab(i,j,k);
            }

            umac_fab(i,j,k) = umac;
          }
        }

        if(idx < vbx_npoints)
        {
          int k = idx / (vbx_len.x*vbx_len.y);
          int j = (idx - k*(vbx_len.x*vbx_len.y)) / (vbx_len.x);
          int i = (idx - k*(vbx_len.x*vbx_len.y)) - j*vbx_len.x;

          i += vbx_lo.x;
          j += vbx_lo.y;
          k += vbx_lo.z;

          // Y-faces
          if (apy_fab(i,j,k) > 0.0)
          {
            Real xf = fcy_fab(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
            Real zf = fcy_fab(i,j,k,1);

            Real delta_x = xf - ccc_fab(i,j,k,0);
            Real delta_y = .5 + ccc_fab(i,j,k,1);
            Real delta_z = zf - ccc_fab(i,j,k,2);

            const Real ccvel_pls = ccvel_fab(i,j,k,1);
            const Real ccvel_mns = ccvel_fab(i,j-1,k,1);

            Real cc_vmax = amrex::max(ccvel_pls, ccvel_mns);
            Real cc_vmin = amrex::min(ccvel_pls, ccvel_mns);

            Real vpls = ccvel_pls - delta_y * yslopes_fab(i,j,k,1)
                                  + delta_x * xslopes_fab(i,j,k,1)
                                  + delta_z * zslopes_fab(i,j,k,1);

            vpls = amrex::min(vpls, cc_vmax);
            vpls = amrex::max(vpls, cc_vmin);

            delta_x = xf - ccc_fab(i,j-1,k,0);
            delta_y = .5 - ccc_fab(i,j-1,k,1);
            delta_z = zf - ccc_fab(i,j-1,k,2);

            Real vmns = ccvel_mns + delta_y * yslopes_fab(i,j-1,k,1)
                                  + delta_x * xslopes_fab(i,j-1,k,1)
                                  + delta_z * zslopes_fab(i,j-1,k,1);

            vmns = amrex::min(vmns, cc_vmax);
            vmns = amrex::max(vmns, cc_vmin);

            Real vmac(0);

            if (vmns >= 0 or vpls <= 0) {
              Real avg = .5 * (vpls + vmns);

              if (avg >= small_vel) {
                vmac = vmns;
              }
              else if (avg <= -small_vel) {
                vmac = vpls;
              }

              vmac *= epy_fab(i,j,k);
            }

            vmac_fab(i,j,k) = vmac;
          }
        }

        if(idx < wbx_npoints)
        {
          int k = idx / (wbx_len.x*wbx_len.y);
          int j = (idx - k*(wbx_len.x*wbx_len.y)) / (wbx_len.x);
          int i = (idx - k*(wbx_len.x*wbx_len.y)) - j*wbx_len.x;

          i += wbx_lo.x;
          j += wbx_lo.y;
          k += wbx_lo.z;

          // Z-faces
          if (apz_fab(i,j,k) > 0.0)
          {
            Real xf = fcz_fab(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
            Real yf = fcz_fab(i,j,k,1);

            Real delta_x = xf - ccc_fab(i,j,k,0);
            Real delta_y = yf - ccc_fab(i,j,k,1);
            Real delta_z = .5 + ccc_fab(i,j,k,2);

            const Real ccvel_pls = ccvel_fab(i,j,k,2);
            const Real ccvel_mns = ccvel_fab(i,j,k-1,2);

            Real cc_wmax = amrex::max(ccvel_pls, ccvel_mns);
            Real cc_wmin = amrex::min(ccvel_pls, ccvel_mns);

            Real wpls = ccvel_pls - delta_z * zslopes_fab(i,j,k,2)
                                  + delta_x * xslopes_fab(i,j,k,2)
                                  + delta_y * yslopes_fab(i,j,k,2);

            wpls = amrex::min(wpls, cc_wmax);
            wpls = amrex::max(wpls, cc_wmin);

            delta_x = xf - ccc_fab(i,j,k-1,0);
            delta_y = yf - ccc_fab(i,j,k-1,1);
            delta_z = .5 - ccc_fab(i,j,k-1,2);

            Real wmns = ccvel_mns + delta_z * zslopes_fab(i,j,k-1,2)
                                  + delta_x * xslopes_fab(i,j,k-1,2)
                                  + delta_y * yslopes_fab(i,j,k-1,2);

            wmns = amrex::min(wmns, cc_wmax);
            wmns = amrex::max(wmns, cc_wmin);

            Real wmac(0);

            if ( wmns >= 0 or wpls <= 0) {
              Real avg = .5 * (wpls + wmns);

              if (avg >= small_vel) {
                wmac = wmns;
              }
              else if (avg <= -small_vel) {
                wmac = wpls;
              }

              wmac *= epz_fab(i,j,k);
            }

            wmac_fab(i,j,k) = wmac;
          }
        }
      });
    } // Cut cells
  } // MFIter

  delete ep_face[0];
  delete ep_face[1];
  delete ep_face[2];
}
