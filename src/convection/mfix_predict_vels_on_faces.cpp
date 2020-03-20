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

    Box domain(geom[lev].Domain());

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
    average_cellcenter_to_face(ep_face, *ep_in[lev], geom[lev]);

    ep_face[0]->FillBoundary();
    ep_face[1]->FillBoundary();
    ep_face[2]->FillBoundary();

    // ****************************************************************************
    // We will store the left and right states in arrays for interpolation
    // ****************************************************************************

    MultiFab upls(ep_face[0]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
    MultiFab umns(ep_face[0]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);

    MultiFab vpls(ep_face[1]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
    MultiFab vmns(ep_face[1]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);

    MultiFab wpls(ep_face[2]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
    MultiFab wmns(ep_face[2]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);

    // We need this just to avoid FPE (eg for outflow faces)
    upls.setVal(covered_val);
    umns.setVal(covered_val);
    vpls.setVal(covered_val);
    vmns.setVal(covered_val);
    wpls.setVal(covered_val);
    wmns.setVal(covered_val);

    // ****************************************************************************
    // First compute the slopes
    // ****************************************************************************
    int slopes_comp = 0;

    mfix_compute_slopes(lev, time, *vel_in[lev],
                        get_xslopes_u(), get_yslopes_u(), get_zslopes_u(),
                        slopes_comp);

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
       else if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
       {

          // Face-centered left and right states
          const auto& upls_fab = upls.array(mfi);
          const auto& vpls_fab = vpls.array(mfi);
          const auto& wpls_fab = wpls.array(mfi);
          const auto& umns_fab = umns.array(mfi);
          const auto& vmns_fab = vmns.array(mfi);
          const auto& wmns_fab = wmns.array(mfi);

          // Face-centered ep
          const auto& epx_fab = (ep_face[0])->array(mfi);
          const auto& epy_fab = (ep_face[1])->array(mfi);
          const auto& epz_fab = (ep_face[2])->array(mfi);

          amrex::ParallelFor(ubx,
            [small_vel,ccvel_fab,epx_fab,xslopes_fab,upls_fab,umns_fab,umac_fab]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              // X-faces
              upls_fab(i,j,k) = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
              umns_fab(i,j,k) = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
              if ( umns_fab(i,j,k) < 0.0 && upls_fab(i,j,k) > 0.0 ) {
                 umac_fab(i,j,k) = 0.0;
              } else {
                 Real avg = 0.5 * ( upls_fab(i,j,k) + umns_fab(i,j,k) );
                 if ( std::abs(avg) <  small_vel) { umac_fab(i,j,k) = 0.0;
                 } else if (avg >= 0)             { umac_fab(i,j,k) = umns_fab(i,j,k);
                 } else                           { umac_fab(i,j,k) = upls_fab(i,j,k);
                 }
                 umac_fab(i,j,k) *= epx_fab(i,j,k);
              }
          });

          amrex::ParallelFor(vbx,
            [small_vel,ccvel_fab,epy_fab,yslopes_fab,vpls_fab,vmns_fab,vmac_fab]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              // Y-faces
              vpls_fab(i,j,k) = ccvel_fab(i,j  ,k,1) - 0.5 * yslopes_fab(i,j  ,k,1);
              vmns_fab(i,j,k) = ccvel_fab(i,j-1,k,1) + 0.5 * yslopes_fab(i,j-1,k,1);
              if ( vmns_fab(i,j,k) < 0.0 && vpls_fab(i,j,k) > 0.0 ) {
                 vmac_fab(i,j,k) = 0.0;
              } else {
                 Real avg = 0.5 * ( vpls_fab(i,j,k) + vmns_fab(i,j,k) );
                 if ( std::abs(avg) <  small_vel) { vmac_fab(i,j,k) = 0.0;
                 } else if (avg >= 0)             { vmac_fab(i,j,k) = vmns_fab(i,j,k);
                 } else                           { vmac_fab(i,j,k) = vpls_fab(i,j,k);
                 }
                 vmac_fab(i,j,k) *= epy_fab(i,j,k);
              }
          });

          amrex::ParallelFor(wbx,
            [small_vel,ccvel_fab,epz_fab,zslopes_fab,wpls_fab,wmns_fab,wmac_fab]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              // Z-faces
              wpls_fab(i,j,k) = ccvel_fab(i,j,k  ,2) - 0.5 * zslopes_fab(i,j,k  ,2);
              wmns_fab(i,j,k) = ccvel_fab(i,j,k-1,2) + 0.5 * zslopes_fab(i,j,k-1,2);
              if ( wmns_fab(i,j,k) < 0.0 && wpls_fab(i,j,k) > 0.0 ) {
                 wmac_fab(i,j,k) = 0.0;
              } else {
                 Real avg = 0.5 * ( wpls_fab(i,j,k) + wmns_fab(i,j,k) );
                 if ( std::abs(avg) <  small_vel) { wmac_fab(i,j,k) = 0.0;
                 } else if (avg >= 0)             { wmac_fab(i,j,k) = wmns_fab(i,j,k);
                 } else                           { wmac_fab(i,j,k) = wpls_fab(i,j,k);
                 }
                 wmac_fab(i,j,k) *= epz_fab(i,j,k);
              }
          });

       // Cut cells in this FAB
       }
       else
       {
          const auto& ccvel_fab = vel_in[lev]->array(mfi);

          // Face centroids
          const auto& fcx_fab = facecent[0]->array(mfi);
          const auto& fcy_fab = facecent[1]->array(mfi);
          const auto& fcz_fab = facecent[2]->array(mfi);

          // Cell centroids
          const auto& ccc_fab = cellcent.array(mfi);

          // Cell-based slopes
          const auto& xslopes_fab = m_leveldata[lev]->xslopes_u->array(mfi);
          const auto& yslopes_fab = m_leveldata[lev]->yslopes_u->array(mfi);
          const auto& zslopes_fab = m_leveldata[lev]->zslopes_u->array(mfi);

          // Face-centered ep
          const auto& epx_fab = (ep_face[0])->array(mfi);
          const auto& epy_fab = (ep_face[1])->array(mfi);
          const auto& epz_fab = (ep_face[2])->array(mfi);

          // Face-centered left and right states
          const auto& upls_fab = upls.array(mfi);
          const auto& vpls_fab = vpls.array(mfi);
          const auto& wpls_fab = wpls.array(mfi);
          const auto& umns_fab = umns.array(mfi);
          const auto& vmns_fab = vmns.array(mfi);
          const auto& wmns_fab = wmns.array(mfi);

          // Face-centered areas
          const auto& apx_fab = areafrac[0]->array(mfi);
          const auto& apy_fab = areafrac[1]->array(mfi);
          const auto& apz_fab = areafrac[2]->array(mfi);

          // This FAB has cut cells -- we predict from cell centroids to face centroids
          amrex::ParallelFor(ubx,
            [apx_fab,fcx_fab,epx_fab,ccc_fab,upls_fab,umns_fab,ccvel_fab,xslopes_fab,yslopes_fab,zslopes_fab,umac_fab,small_vel]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              // X-faces
              if (apx_fab(i,j,k) > 0.0)
              {
                 Real yf = fcx_fab(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
                 Real zf = fcx_fab(i,j,k,1);

                 Real xc = ccc_fab(i,j,k,0); // centroid of cell (i,j,k)
                 Real yc = ccc_fab(i,j,k,1);
                 Real zc = ccc_fab(i,j,k,2);

                 Real delta_x = 0.5 + xc;
                 Real delta_y = yf  - yc;
                 Real delta_z = zf  - zc;

                 Real cc_umax = std::max(ccvel_fab(i,j,k,0), ccvel_fab(i-1,j,k,0));
                 Real cc_umin = std::min(ccvel_fab(i,j,k,0), ccvel_fab(i-1,j,k,0));

                 upls_fab(i,j,k) = ccvel_fab(i  ,j,k,0) - delta_x * xslopes_fab(i,j,k,0) 
                                                        + delta_y * yslopes_fab(i,j,k,0) 
                                                        + delta_z * zslopes_fab(i,j,k,0);

                 upls_fab(i,j,k) = std::min(upls_fab(i,j,k), cc_umax);
                 upls_fab(i,j,k) = std::max(upls_fab(i,j,k), cc_umin);

                 xc = ccc_fab(i-1,j,k,0); // centroid of cell (i-1,j,k)
                 yc = ccc_fab(i-1,j,k,1);
                 zc = ccc_fab(i-1,j,k,2);

                 delta_x = 0.5 - xc;
                 delta_y = yf  - yc;
                 delta_z = zf  - zc;

                 umns_fab(i,j,k) = ccvel_fab(i-1,j,k,0) + delta_x * xslopes_fab(i-1,j,k,0) 
                                                        + delta_y * yslopes_fab(i-1,j,k,0) 
                                                        + delta_z * zslopes_fab(i-1,j,k,0);

                 umns_fab(i,j,k) = std::min(umns_fab(i,j,k), cc_umax);
                 umns_fab(i,j,k) = std::max(umns_fab(i,j,k), cc_umin);

                 if ( umns_fab(i,j,k) < 0.0 && upls_fab(i,j,k) > 0.0 ) {
                    umac_fab(i,j,k) = 0.0;
                 } else {
                    Real avg = 0.5 * ( upls_fab(i,j,k) + umns_fab(i,j,k) );
                    if ( std::abs(avg) <  small_vel) { umac_fab(i,j,k) = 0.0;
                    } else if (avg >= 0)             { umac_fab(i,j,k) = umns_fab(i,j,k);
                    } else                           { umac_fab(i,j,k) = upls_fab(i,j,k);
                    }
                    umac_fab(i,j,k) *= epx_fab(i,j,k);
                 }
              }
          });

          amrex::ParallelFor(vbx,
            [apy_fab,fcy_fab,epy_fab,ccc_fab,vpls_fab,vmns_fab,ccvel_fab,xslopes_fab,yslopes_fab,zslopes_fab,vmac_fab,small_vel]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              // Y-faces
              if (apy_fab(i,j,k) > 0.0)
              {
                 Real xf = fcy_fab(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
                 Real zf = fcy_fab(i,j,k,1);

                 Real xc = ccc_fab(i,j,k,0); // centroid of cell (i,j,k)
                 Real yc = ccc_fab(i,j,k,1);
                 Real zc = ccc_fab(i,j,k,2);

                 Real delta_x = xf  - xc;
                 Real delta_y = 0.5 + yc;
                 Real delta_z = zf  - zc;

                 Real cc_umax = std::max(ccvel_fab(i,j,k,1), ccvel_fab(i,j-1,k,1));
                 Real cc_umin = std::min(ccvel_fab(i,j,k,1), ccvel_fab(i,j-1,k,1));

                 vpls_fab(i,j,k) = ccvel_fab(i,j  ,k,1) - delta_y * yslopes_fab(i,j,k,1) 
                                                        + delta_x * xslopes_fab(i,j,k,1) 
                                                        + delta_z * zslopes_fab(i,j,k,1);

                 vpls_fab(i,j,k) = std::min(vpls_fab(i,j,k), cc_umax);
                 vpls_fab(i,j,k) = std::max(vpls_fab(i,j,k), cc_umin);

                 xc = ccc_fab(i,j-1,k,0); // centroid of cell (i,j-1,k)
                 yc = ccc_fab(i,j-1,k,1);
                 zc = ccc_fab(i,j-1,k,2);

                 delta_x = xf  - xc;
                 delta_y = 0.5 - yc;
                 delta_z = zf  - zc;

                 vmns_fab(i,j,k) = ccvel_fab(i,j-1,k,1) + delta_y * yslopes_fab(i,j-1,k,1) 
                                                        + delta_x * xslopes_fab(i,j-1,k,1) 
                                                        + delta_z * zslopes_fab(i,j-1,k,1);

                 vmns_fab(i,j,k) = std::min(vmns_fab(i,j,k), cc_umax);
                 vmns_fab(i,j,k) = std::max(vmns_fab(i,j,k), cc_umin);

                 if ( vmns_fab(i,j,k) < 0.0 && vpls_fab(i,j,k) > 0.0 ) {
                    vmac_fab(i,j,k) = 0.0;
                 } else {
                    Real avg = 0.5 * ( vpls_fab(i,j,k) + vmns_fab(i,j,k) );
                    if ( std::abs(avg) <  small_vel) { vmac_fab(i,j,k) = 0.0;
                    } else if (avg >= 0)             { vmac_fab(i,j,k) = vmns_fab(i,j,k);
                    } else                           { vmac_fab(i,j,k) = vpls_fab(i,j,k);
                    }
                    vmac_fab(i,j,k) *= epy_fab(i,j,k);
                 }
              }
          });

          amrex::ParallelFor(wbx,
            [apz_fab,fcz_fab,epz_fab,ccc_fab,wpls_fab,wmns_fab,ccvel_fab,xslopes_fab,yslopes_fab,zslopes_fab,wmac_fab,small_vel]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              // Z-faces
              if (apz_fab(i,j,k) > 0.0) 
              {
                 Real xf = fcz_fab(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
                 Real yf = fcz_fab(i,j,k,1);

                 Real xc = ccc_fab(i,j,k,0); // centroid of cell (i,j,k)
                 Real yc = ccc_fab(i,j,k,1);
                 Real zc = ccc_fab(i,j,k,2);

                 Real delta_x = xf  - xc;
                 Real delta_y = yf  - yc;
                 Real delta_z = 0.5 + zc;

                 Real cc_umax = std::max(ccvel_fab(i,j,k,2), ccvel_fab(i,j,k-1,2));
                 Real cc_umin = std::min(ccvel_fab(i,j,k,2), ccvel_fab(i,j,k-1,2));

                 wpls_fab(i,j,k) = ccvel_fab(i,j,k  ,2) - delta_z * zslopes_fab(i,j,k,2) 
                                                        + delta_x * xslopes_fab(i,j,k,2) 
                                                        + delta_y * yslopes_fab(i,j,k,2);

                 wpls_fab(i,j,k) = std::min(wpls_fab(i,j,k), cc_umax);
                 wpls_fab(i,j,k) = std::max(wpls_fab(i,j,k), cc_umin);

                 xc = ccc_fab(i,j,k-1,0); // centroid of cell (i,j,k-1)
                 yc = ccc_fab(i,j,k-1,1);
                 zc = ccc_fab(i,j,k-1,2);

                 delta_x = xf  - xc;
                 delta_y = yf  - yc;
                 delta_z = 0.5 - zc;

                 wmns_fab(i,j,k) = ccvel_fab(i,j,k-1,2) + delta_z * zslopes_fab(i,j,k-1,2) 
                                                        + delta_x * xslopes_fab(i,j,k-1,2) 
                                                        + delta_y * yslopes_fab(i,j,k-1,2);

                 wmns_fab(i,j,k) = std::min(wmns_fab(i,j,k), cc_umax);
                 wmns_fab(i,j,k) = std::max(wmns_fab(i,j,k), cc_umin);

                 if ( wmns_fab(i,j,k) < 0.0 && wpls_fab(i,j,k) > 0.0 ) {
                    wmac_fab(i,j,k) = 0.0;
                 } else {
                    Real avg = 0.5 * ( wpls_fab(i,j,k) + wmns_fab(i,j,k) );
                    if ( std::abs(avg) <  small_vel) { wmac_fab(i,j,k) = 0.0;
                    } else if (avg >= 0)             { wmac_fab(i,j,k) = wmns_fab(i,j,k);
                    } else                           { wmac_fab(i,j,k) = wpls_fab(i,j,k);
                    }
                    wmac_fab(i,j,k) *= epz_fab(i,j,k);
                 }
              }
          });
       } // Cut cells
    } // MFIter

    delete ep_face[0];
    delete ep_face[1];
    delete ep_face[2];
}
