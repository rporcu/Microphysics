#include <mfix.H>

void
mfix::mfix_predict_vels_on_faces (int lev, Real time,
                                  Vector< std::unique_ptr<MultiFab> >& vel_in,
                                  Vector< std::unique_ptr<MultiFab> >& ep_u_mac,
                                  Vector< std::unique_ptr<MultiFab> >& ep_v_mac,
                                  Vector< std::unique_ptr<MultiFab> >& ep_w_mac,
                                  Vector< std::unique_ptr<MultiFab> >& ep_in)

{
    BL_PROFILE("mfix::mfix_predict_vels_on_faces");

    Array< std::unique_ptr<MultiFab>, AMREX_SPACEDIM> ep_face;

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
       Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
       Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;

       areafrac = ebfactory[lev]->getAreaFrac();
       facecent = ebfactory[lev]->getFaceCent();
       const auto& cellcent = ebfactory[lev]->getCentroid();

       Real small_vel = 1.e-10;
       Real  huge_vel = 1.e100;


       // ****************************************************************************
       // We will need ep on face centers to interpolate to face centroids below
       // ****************************************************************************

       ep_face[0].reset(new MultiFab(ep_u_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]));
       ep_face[1].reset(new MultiFab(ep_v_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]));
       ep_face[2].reset(new MultiFab(ep_w_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]));

       // This is to make sure ep_face is defined everywhere
       ep_face[0]->setVal(covered_val);
       ep_face[1]->setVal(covered_val);
       ep_face[2]->setVal(covered_val);

       ep_in[lev]->FillBoundary(geom[lev].periodicity());
       average_cellcenter_to_face(GetArrOfPtrs(ep_face), *ep_in[lev], geom[lev]);

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
       mfix_compute_slopes(lev, time, *vel_in[lev], xslopes_u, yslopes_u, zslopes_u, slopes_comp);

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

          Box ubx_grown = mfi.growntilebox(IntVect::TheDimensionVector(0));
          Box vbx_grown = mfi.growntilebox(IntVect::TheDimensionVector(1));
          Box wbx_grown = mfi.growntilebox(IntVect::TheDimensionVector(2));

          const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
          const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

          // Cell-centered velocity
          const auto& ccvel_fab = vel_in[lev]->array(mfi);

          // Cell-centered slopes
          const auto& xslopes_fab = (xslopes_u[lev])->array(mfi);
          const auto& yslopes_fab = (yslopes_u[lev])->array(mfi);
          const auto& zslopes_fab = (zslopes_u[lev])->array(mfi);

          // Face-centered velocity components
          const auto& umac_fab = (ep_u_mac[lev])->array(mfi);
          const auto& vmac_fab = (ep_v_mac[lev])->array(mfi);
          const auto& wmac_fab = (ep_w_mac[lev])->array(mfi);

          if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
          {
            Real val = 1.2345e300;

            (*ep_u_mac[lev])[mfi].setVal(val, ubx, 0, 1);
            (*ep_v_mac[lev])[mfi].setVal(val, vbx, 0, 1);
            (*ep_w_mac[lev])[mfi].setVal(val, wbx, 0, 1);
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
             const auto& xslopes_fab = (xslopes_u[lev])->array(mfi);
             const auto& yslopes_fab = (yslopes_u[lev])->array(mfi);
             const auto& zslopes_fab = (zslopes_u[lev])->array(mfi);

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

                    if (std::abs(delta_x) > 1.0) amrex::Print() << "HI:BOGUS DX IN X-EXTRAP " << delta_x << " " << i << " " << j << " " << k << std::endl;
                    if (std::abs(delta_y) > 0.5) amrex::Print() << "HI:BOGUS DY IN X-EXTRAP " << delta_y << " " << i << " " << j << " " << k << std::endl;
                    if (std::abs(delta_z) > 0.5) amrex::Print() << "HI:BOGUS DZ IN X-EXTRAP " << delta_z << " " << i << " " << j << " " << k << std::endl;

                    upls_fab(i,j,k) = ccvel_fab(i  ,j,k,0) - delta_x * xslopes_fab(i,j,k,0);
                                                           + delta_y * yslopes_fab(i,j,k,0);
                                                           + delta_z * zslopes_fab(i,j,k,0);

                         xc = ccc_fab(i-1,j,k,0); // centroid of cell (i-1,j,k)
                         yc = ccc_fab(i-1,j,k,1);
                         zc = ccc_fab(i-1,j,k,2);

                         delta_x = 0.5 - xc;
                         delta_y = yf  - yc;
                         delta_z = zf  - zc;

                    if (std::abs(delta_x) > 1.0) amrex::Print() << "LO:BOGUS DX IN X-EXTRAP " << delta_x << " " << i << " " << j << " " << k << std::endl;
                    if (std::abs(delta_y) > 0.5) amrex::Print() << "LO:BOGUS DY IN X-EXTRAP " << delta_y << " " << i << " " << j << " " << k << std::endl;
                    if (std::abs(delta_z) > 0.5) amrex::Print() << "LO:BOGUS DZ IN X-EXTRAP " << delta_z << " " << i << " " << j << " " << k << std::endl;

                    umns_fab(i,j,k) = ccvel_fab(i-1,j,k,0) + delta_x * xslopes_fab(i-1,j,k,0);
                                                           + delta_y * yslopes_fab(i-1,j,k,0);
                                                           + delta_z * zslopes_fab(i-1,j,k,0);

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

                    if (std::abs(delta_x) > 0.5) amrex::Print() << "HI:BOGUS DX IN X-EXTRAP " << delta_x << std::endl;
                    if (std::abs(delta_y) > 1.0) amrex::Print() << "HI:BOGUS DY IN X-EXTRAP " << delta_y << std::endl;
                    if (std::abs(delta_z) > 0.5) amrex::Print() << "HI:BOGUS DZ IN X-EXTRAP " << delta_z << std::endl;

                    vpls_fab(i,j,k) = ccvel_fab(i,j  ,k,1) - delta_y * yslopes_fab(i,j,k,1);
                                                           + delta_x * xslopes_fab(i,j,k,1);
                                                           + delta_z * zslopes_fab(i,j,k,1);

                         xc = ccc_fab(i,j-1,k,0); // centroid of cell (i,j-1,k)
                         yc = ccc_fab(i,j-1,k,1);
                         zc = ccc_fab(i,j-1,k,2);

                         delta_x = xf  - xc;
                         delta_y = 0.5 - yc;
                         delta_z = zf  - zc;

                    if (std::abs(delta_x) > 0.5) amrex::Print() << "LO:BOGUS DX IN X-EXTRAP " << delta_x << std::endl;
                    if (std::abs(delta_y) > 1.0) amrex::Print() << "LO:BOGUS DY IN X-EXTRAP " << delta_y << std::endl;
                    if (std::abs(delta_z) > 0.5) amrex::Print() << "LO:BOGUS DZ IN X-EXTRAP " << delta_z << std::endl;

                    vmns_fab(i,j,k) = ccvel_fab(i,j-1,k,1) + delta_y * yslopes_fab(i,j-1,k,1);
                                                           + delta_x * xslopes_fab(i,j-1,k,1);
                                                           + delta_z * zslopes_fab(i,j-1,k,1);

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

                    if (std::abs(delta_x) > 0.5) amrex::Print() << "HI:BOGUS DX IN X-EXTRAP " << delta_x << std::endl;
                    if (std::abs(delta_y) > 0.5) amrex::Print() << "HI:BOGUS DY IN X-EXTRAP " << delta_y << std::endl;
                    if (std::abs(delta_z) > 1.0) amrex::Print() << "HI:BOGUS DZ IN X-EXTRAP " << delta_z << std::endl;

                    wpls_fab(i,j,k) = ccvel_fab(i,j,k  ,2) - delta_z * zslopes_fab(i,j,k,2);
                                                           + delta_x * xslopes_fab(i,j,k,2);
                                                           + delta_y * yslopes_fab(i,j,k,2);

                         xc = ccc_fab(i,j,k-1,0); // centroid of cell (i,j,k-1)
                         yc = ccc_fab(i,j,k-1,1);
                         zc = ccc_fab(i,j,k-1,2);

                         delta_x = xf  - xc;
                         delta_y = yf  - yc;
                         delta_z = 0.5 - zc;

                    wmns_fab(i,j,k) = ccvel_fab(i,j,k-1,2) + delta_z * zslopes_fab(i,j,k-1,2);
                                                           + delta_x * xslopes_fab(i,j,k-1,2);
                                                           + delta_y * yslopes_fab(i,j,k-1,2);

                    if (std::abs(delta_x) > 0.5) amrex::Print() << "LO:BOGUS DX IN X-EXTRAP " << delta_x << std::endl;
                    if (std::abs(delta_y) > 0.5) amrex::Print() << "LO:BOGUS DY IN X-EXTRAP " << delta_y << std::endl;
                    if (std::abs(delta_z) > 1.0) amrex::Print() << "LO:BOGUS DZ IN X-EXTRAP " << delta_z << std::endl;

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

#if 0
       // ****************************************************************************
       // Make sure to fill face-centered values outside our grid before interpolating
       // ****************************************************************************
       upls.FillBoundary(geom[lev].periodicity());
       umns.FillBoundary(geom[lev].periodicity());
       vpls.FillBoundary(geom[lev].periodicity());
       vmns.FillBoundary(geom[lev].periodicity());
       wpls.FillBoundary(geom[lev].periodicity());
       wmns.FillBoundary(geom[lev].periodicity());

       // ****************************************************************************
       // Do interpolation to centroids -- only for cut cells
       // ****************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
          // Tilebox
          const Box  bx = mfi.tilebox();

          Box ubx = mfi.nodaltilebox(0);
          Box vbx = mfi.nodaltilebox(1);
          Box wbx = mfi.nodaltilebox(2);
      
          // Check efficiently if this tile contains any eb stuff

          const EBFArrayBox& vel_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
          const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

          if ( !(flags.getType(amrex::grow(bx,0)) == FabType::covered or
                 flags.getType(amrex::grow(bx,1)) == FabType::regular ) )
          {
             // Face-centered velocity components
             const auto& umac_fab = (ep_u_mac[lev])->array(mfi);
             const auto& vmac_fab = (ep_v_mac[lev])->array(mfi);
             const auto& wmac_fab = (ep_w_mac[lev])->array(mfi);

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

             // Face-centered ep
             const auto& epx_fab = (ep_face[0])->array(mfi);
             const auto& epy_fab = (ep_face[1])->array(mfi);
             const auto& epz_fab = (ep_face[2])->array(mfi);

             // Face centroids
             const auto& fcx_fab = facecent[0]->array(mfi);
             const auto& fcy_fab = facecent[1]->array(mfi);
             const auto& fcz_fab = facecent[2]->array(mfi);

             const auto& ccm_fab = cc_mask.const_array(mfi);

             amrex::ParallelFor(ubx,
               [small_vel,huge_vel,apx_fab,epx_fab,ccm_fab,fcx_fab,upls_fab,umns_fab,umac_fab]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                if (apx_fab(i,j,k) == 0.0)

                    umac_fab(i,j,k) = huge_vel;

                else if (apx_fab(i,j,k) < 1.0) {

                    int jj = j + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,0)));
                    int kk = k + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,1)));

                    Real fracy = (ccm_fab(i-1,jj,k) or ccm_fab(i,jj,k)) ? std::abs(fcx_fab(i,j,k,0)) : 0.0;
                    Real fracz = (ccm_fab(i-1,j,kk) or ccm_fab(i,j,kk)) ? std::abs(fcx_fab(i,j,k,1)) : 0.0;

                    Real upls_on_centroid = (1.0-fracy)*(1.0-fracz)*upls_fab(i, j,k ) +
                                                 fracy *(1.0-fracz)*upls_fab(i,jj,k ) +
                                                 fracz *(1.0-fracy)*upls_fab(i, j,kk) +
                                                 fracy *     fracz *upls_fab(i,jj,kk);
                    Real umns_on_centroid = (1.0-fracy)*(1.0-fracz)*umns_fab(i, j,k ) +
                                                 fracy *(1.0-fracz)*umns_fab(i,jj,k ) +
                                                 fracz *(1.0-fracy)*umns_fab(i, j,kk) +
                                                 fracy *     fracz *umns_fab(i,jj,kk);

                    Real   ep_on_centroid = (1.0-fracy)*(1.0-fracz)*epx_fab(i, j,k ) +
                                                 fracy *(1.0-fracz)*epx_fab(i,jj,k ) +
                                                 fracz *(1.0-fracy)*epx_fab(i, j,kk) +
                                                 fracy *     fracz *epx_fab(i,jj,kk);

                    if ( umns_on_centroid < 0.0 and upls_on_centroid > 0.0 ) {
                       umac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( upls_on_centroid + umns_on_centroid );
                       if ( std::abs(avg) < small_vel) { umac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)            { umac_fab(i,j,k) = umns_on_centroid;
                       } else                          { umac_fab(i,j,k) = upls_on_centroid;
                       }
                       umac_fab(i,j,k) *= ep_on_centroid;
                    }

                 } else {

                    if ( umns_fab(i,j,k) < 0.0 and upls_fab(i,j,k) > 0.0 ) {
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
               [small_vel,huge_vel,apy_fab,epy_fab,ccm_fab,fcy_fab,vpls_fab,vmns_fab,vmac_fab]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                if (apy_fab(i,j,k) == 0.0) {

                    vmac_fab(i,j,k) = huge_vel;

                } else if (apy_fab(i,j,k) < 1.0) {
                    int ii = i + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,0)));
                    int kk = k + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,1)));

                    Real fracx = (ccm_fab(ii,j-1,k ) or ccm_fab(ii,j,k)) ? std::abs(fcy_fab(i,j,k,0)) : 0.0;
                    Real fracz = (ccm_fab( i,j-1,kk) or ccm_fab(i,j,kk)) ? std::abs(fcy_fab(i,j,k,1)) : 0.0;

                    Real vpls_on_centroid = (1.0-fracx)*(1.0-fracz)*vpls_fab(i ,j,k ) +
                                                 fracx *(1.0-fracz)*vpls_fab(ii,j,k ) +
                                                 fracz *(1.0-fracx)*vpls_fab(i ,j,kk) +
                                                 fracx *     fracz *vpls_fab(ii,j,kk);
                    Real vmns_on_centroid = (1.0-fracx)*(1.0-fracz)*vmns_fab(i ,j,k ) +
                                                 fracx *(1.0-fracz)*vmns_fab(ii,j,k ) +
                                                 fracz *(1.0-fracx)*vmns_fab(i ,j,kk) +
                                                 fracx *     fracz *vmns_fab(ii,j,kk);
                    Real   ep_on_centroid = (1.0-fracx)*(1.0-fracz)* epy_fab(i ,j,k ) +
                                                 fracx *(1.0-fracz)* epy_fab(ii,j,k ) +
                                                 fracz *(1.0-fracx)* epy_fab(i ,j,kk) +
                                                 fracx *     fracz * epy_fab(ii,j,kk);

                    if ( vmns_on_centroid < 0.0 and vpls_on_centroid > 0.0 ) {
                       vmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( vpls_on_centroid + vmns_on_centroid );
                       if ( std::abs(avg) < small_vel) { vmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)            { vmac_fab(i,j,k) = vmns_on_centroid;
                       } else                          { vmac_fab(i,j,k) = vpls_on_centroid;
                       }
                       vmac_fab(i,j,k) *= ep_on_centroid;
                    }

                 } else {

                    if ( vmns_fab(i,j,k) < 0.0 and vpls_fab(i,j,k) > 0.0 ) {
                       vmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( vpls_fab(i,j,k) + vmns_fab(i,j,k) );
                       if ( std::abs(avg) < small_vel) { vmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)            { vmac_fab(i,j,k) = vmns_fab(i,j,k);
                       } else                          { vmac_fab(i,j,k) = vpls_fab(i,j,k);
                       }
                       vmac_fab(i,j,k) *= epy_fab(i,j,k);
                    }
                 }
             });

             amrex::ParallelFor(wbx,
               [small_vel,huge_vel,apz_fab,epz_fab,ccm_fab,fcz_fab,wpls_fab,wmns_fab,wmac_fab]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                if (apz_fab(i,j,k) == 0.0)

                    wmac_fab(i,j,k) = huge_vel;

                else if (apz_fab(i,j,k) < 1.0) {

                    int ii = i + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,0)));
                    int jj = j + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,1)));

                    Real fracx = (ccm_fab(ii,j,k-1) or ccm_fab(ii,j,k)) ? std::abs(fcz_fab(i,j,k,0)) : 0.0;
                    Real fracy = (ccm_fab(i,jj,k-1) or ccm_fab(i,jj,k)) ? std::abs(fcz_fab(i,j,k,1)) : 0.0;

                    Real wpls_on_centroid = (1.0-fracx)*(1.0-fracy)*wpls_fab(i ,j ,k) +
                                                 fracx *(1.0-fracy)*wpls_fab(ii,j ,k) +
                                                 fracy *(1.0-fracx)*wpls_fab(i ,jj,k) +
                                                 fracx *     fracy *wpls_fab(ii,jj,k);
                    Real wmns_on_centroid = (1.0-fracx)*(1.0-fracy)*wmns_fab(i ,j ,k) +
                                                 fracx *(1.0-fracy)*wmns_fab(ii,j ,k) +
                                                 fracy *(1.0-fracx)*wmns_fab(i ,jj,k) +
                                                 fracx *     fracy *wmns_fab(ii,jj,k);
                    Real   ep_on_centroid = (1.0-fracx)*(1.0-fracy)* epz_fab(i ,j ,k) +
                                                 fracx *(1.0-fracy)* epz_fab(ii,j ,k) +
                                                 fracy *(1.0-fracx)* epz_fab(i ,jj,k) +
                                                 fracx *     fracy * epz_fab(ii,jj,k);

                    if ( wmns_on_centroid < 0.0 and wpls_on_centroid > 0.0 ) {
                       wmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( wpls_on_centroid + wmns_on_centroid );
                       if ( std::abs(avg) < small_vel) { wmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)            { wmac_fab(i,j,k) = wmns_on_centroid;
                       } else                          { wmac_fab(i,j,k) = wpls_on_centroid;
                       }
                       wmac_fab(i,j,k) *= ep_on_centroid;
                    }

                 } else {

                    if ( wmns_fab(i,j,k) < 0.0 and wpls_fab(i,j,k) > 0.0 ) {
                       wmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( wpls_fab(i,j,k) + wmns_fab(i,j,k) );
                       if ( std::abs(avg) < small_vel) { wmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)            { wmac_fab(i,j,k) = wmns_fab(i,j,k);
                       } else                          { wmac_fab(i,j,k) = wpls_fab(i,j,k);
                       }
                       wmac_fab(i,j,k) *= epz_fab(i,j,k);
                    }
                 }
             });

          } // Cut cells
       } // MFIter
#endif
}
