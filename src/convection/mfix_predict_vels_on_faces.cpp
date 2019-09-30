#include <mfix.H>

void
mfix::mfix_predict_vels_on_faces ( Real time,
                                   Vector< std::unique_ptr<MultiFab> >& vel_in,
                                   Vector< std::unique_ptr<MultiFab> >& ep_u_mac,
                                   Vector< std::unique_ptr<MultiFab> >& ep_v_mac,
                                   Vector< std::unique_ptr<MultiFab> >& ep_w_mac,
                                   Vector< std::unique_ptr<MultiFab> >& ep_in)

{
    BL_PROFILE("mfix::mfix_predict_vels_on_faces");

    Vector< Array< std::unique_ptr<MultiFab>, AMREX_SPACEDIM> > ep_face;
    ep_face.resize(finest_level+1);

    for (int lev=0; lev < nlev; ++lev)
    {
       Box domain(geom[lev].Domain());   

       iMultiFab cc_mask(grids[lev], dmap[lev], 1, 1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       {
           std::vector< std::pair<int,Box> > isects;
           const std::vector<IntVect>& pshifts = geom[lev].periodicity().shiftIntVect();
           const BoxArray& ba = cc_mask.boxArray();
           for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
           {
               Array4<int> const& fab = cc_mask.array(mfi);
               const Box& bx = mfi.fabbox();
               for (const auto& iv : pshifts)
               {
                   ba.intersections(bx+iv, isects);
                   for (const auto& is : isects)
                   {
                       const Box& b = is.second-iv;
                       AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( b, i, j, k,
                       {
                           fab(i,j,k) = 1;
                       });
                   }
               }
           }
       }

       // We will need these ep on face centers to interpolate to face centroids below
       ep_in[lev]->FillBoundary(geom[lev].periodicity());
       ep_face[lev][0].reset(new MultiFab(ep_u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
       ep_face[lev][1].reset(new MultiFab(ep_v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
       ep_face[lev][2].reset(new MultiFab(ep_w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]));

       // Define ep on face centers
       average_cellcenter_to_face( GetArrOfPtrs(ep_face[lev]), *ep_in[lev], geom[lev] );

       MultiFab upls(ep_face[lev][0]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
       MultiFab umns(ep_face[lev][0]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
       MultiFab vpls(ep_face[lev][1]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
       MultiFab vmns(ep_face[lev][1]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
       MultiFab wpls(ep_face[lev][2]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
       MultiFab wmns(ep_face[lev][2]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);

       // First compute the slopes
       int slopes_comp = 0;
       mfix_compute_slopes(lev, time, *vel_in[lev], xslopes_u, yslopes_u, zslopes_u, slopes_comp);

       // Get EB geometric info
       Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
       Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;

       areafrac  =   ebfactory[lev] -> getAreaFrac();
       facecent  =   ebfactory[lev] -> getFaceCent();

       Real small_vel = 1.e-10;
       Real  huge_vel = 1.e100;

       // Then compute velocity at faces
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

       // Predict to centers
       for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
           // Tilebox
          Box  bx = mfi.tilebox();
          Box ubx = mfi.tilebox(e_x);
          Box vbx = mfi.tilebox(e_y);
          Box wbx = mfi.tilebox(e_z);
      
          // Check efficiently if this tile contains any eb stuff

          const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
          const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

          if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
          {
             (*ep_u_mac[lev])[mfi].setVal( 1.2345e300, ubx, 0, 1);
             (*ep_v_mac[lev])[mfi].setVal( 1.2345e300, vbx, 0, 1);
             (*ep_w_mac[lev])[mfi].setVal( 1.2345e300, wbx, 0, 1);
          }
          else if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
          {

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

             // Face-centered left and right states
             const auto& upls_fab = upls.array(mfi);
             const auto& vpls_fab = vpls.array(mfi);
             const auto& wpls_fab = wpls.array(mfi);
             const auto& umns_fab = umns.array(mfi);
             const auto& vmns_fab = vmns.array(mfi);
             const auto& wmns_fab = wmns.array(mfi);

             // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
             AMREX_HOST_DEVICE_FOR_3D(ubx, i, j, k, 
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
                 }
             });

             AMREX_HOST_DEVICE_FOR_3D(vbx, i, j, k,
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
                 }
             });

             AMREX_HOST_DEVICE_FOR_3D(wbx, i, j, k,
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
                 }
             });

          } else {

             // Cell-centered velocity
             const auto& ccvel_fab = vel_in[lev]->array(mfi);

             // Cell-centered slopes
             const auto& xslopes_fab = (xslopes_u[lev])->array(mfi);
             const auto& yslopes_fab = (yslopes_u[lev])->array(mfi);
             const auto& zslopes_fab = (zslopes_u[lev])->array(mfi);

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

             // This FAB has cut cells
             AMREX_HOST_DEVICE_FOR_3D(ubx, i, j, k, 
             {
                 // X-faces
                 if (apx_fab(i,j,k) > 0.0)
                 {
                    upls_fab(i,j,k) = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                    umns_fab(i,j,k) = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                 }
             });


             AMREX_HOST_DEVICE_FOR_3D(vbx, i, j, k,
             {
                 // Y-faces
                 if (apy_fab(i,j,k) > 0.0)
                 {
                    vpls_fab(i,j,k) = ccvel_fab(i,j  ,k,1) - 0.5 * yslopes_fab(i,j  ,k,1);
                    vmns_fab(i,j,k) = ccvel_fab(i,j-1,k,1) + 0.5 * yslopes_fab(i,j-1,k,1);
                 }
             });

             AMREX_HOST_DEVICE_FOR_3D(wbx, i, j, k,
             {
                 // Z-faces
                 if (apz_fab(i,j,k) > 0.0) {
                    wpls_fab(i,j,k) = ccvel_fab(i,j,k  ,2) - 0.5 * zslopes_fab(i,j,k  ,2);
                    wmns_fab(i,j,k) = ccvel_fab(i,j,k-1,2) + 0.5 * zslopes_fab(i,j,k-1,2);
                 }
             });
          } // Cut cells
       } // MFIter

       // Make sure to fill face-centered values outside our grid before interpolating
       upls.FillBoundary();
       umns.FillBoundary();
       vpls.FillBoundary();
       vmns.FillBoundary();
       wpls.FillBoundary();
       wmns.FillBoundary();

       // Do interpolation to centroids -- only for cut cells
       for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
           // Tilebox
          Box  bx = mfi.tilebox();
          Box ubx = mfi.tilebox(e_x);
          Box vbx = mfi.tilebox(e_y);
          Box wbx = mfi.tilebox(e_z);
      
          // Check efficiently if this tile contains any eb stuff

          const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
          const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

          if ( !(flags.getType(amrex::grow(bx,0)) == FabType::covered ||
                 flags.getType(amrex::grow(bx,1)) == FabType::regular ) )
          {
             // Cell-centered velocity
             const auto& ccvel_fab = vel_in[lev]->array(mfi);

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
             const auto& epx_fab = (ep_face[lev][0])->array(mfi);
             const auto& epy_fab = (ep_face[lev][1])->array(mfi);
             const auto& epz_fab = (ep_face[lev][2])->array(mfi);

             // Face centroids
             const auto& fcx_fab = facecent[0]->array(mfi);
             const auto& fcy_fab = facecent[1]->array(mfi);
             const auto& fcz_fab = facecent[2]->array(mfi);

             const auto& ccm_fab = cc_mask.const_array(mfi);

             AMREX_HOST_DEVICE_FOR_3D(ubx, i, j, k, 
             {
                if (apx_fab(i,j,k) == 0.0)

                    umac_fab(i,j,k) = huge_vel;

                else if (apx_fab(i,j,k) < 1.0) {

                    int jj = j + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,0)));
                    int kk = k + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,1)));

                    Real fracy = (ccm_fab(i-1,jj,k) || ccm_fab(i,jj,k)) ? std::abs(fcx_fab(i,j,k,0)) : 0.0;
                    Real fracz = (ccm_fab(i-1,j,kk) || ccm_fab(i,j,kk)) ? std::abs(fcx_fab(i,j,k,1)) : 0.0;

                    Real upls_on_centroid = (1.0-fracy)*(1.0-fracz)*upls_fab(i, j,k )+
                                                 fracy *(1.0-fracz)*upls_fab(i,jj,k )+
                                                 fracz *(1.0-fracy)*upls_fab(i, j,kk)+
                                                 fracy *     fracz *upls_fab(i,jj,kk);
                    Real umns_on_centroid = (1.0-fracy)*(1.0-fracz)*umns_fab(i, j,k )+
                                                 fracy *(1.0-fracz)*umns_fab(i,jj,k )+
                                                 fracz *(1.0-fracy)*umns_fab(i, j,kk)+
                                                 fracy *     fracz *umns_fab(i,jj,kk);

                    Real   ep_on_centroid = (1.0-fracy)*(1.0-fracz)*epx_fab(i, j,k )+
                                                 fracy *(1.0-fracz)*epx_fab(i,jj,k )+
                                                 fracz *(1.0-fracy)*epx_fab(i, j,kk)+
                                                 fracy *     fracz *epx_fab(i,jj,kk);

                    if ( umns_on_centroid < 0.0 && upls_on_centroid > 0.0 ) {
                       umac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( upls_on_centroid + umns_on_centroid );
                       if ( std::abs(avg) <  small_vel) { umac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { umac_fab(i,j,k) = umns_on_centroid;
                       } else                           { umac_fab(i,j,k) = upls_on_centroid;
                       }
                       umac_fab(i,j,k) *= ep_on_centroid;
                    }

                 } else {

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

             AMREX_HOST_DEVICE_FOR_3D(vbx, i, j, k,
             {
                if (apy_fab(i,j,k) == 0.0)

                    vmac_fab(i,j,k) = huge_vel;

                else if (apy_fab(i,j,k) < 1.0) {

                    int ii = i + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,0)));
                    int kk = k + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,1)));

                    Real fracx = (ccm_fab(ii,j-1,k) || ccm_fab(ii,j,k)) ? std::abs(fcy_fab(i,j,k,0)) : 0.0;
                    Real fracz = (ccm_fab(i,j-1,kk) || ccm_fab(i,j,kk)) ? std::abs(fcy_fab(i,j,k,1)) : 0.0;

                    Real vpls_on_centroid = (1.0-fracx)*(1.0-fracz)*vpls_fab(i ,j,k )+
                                                 fracx *(1.0-fracz)*vpls_fab(ii,j,k )+
                                                 fracz *(1.0-fracx)*vpls_fab(i ,j,kk)+
                                                 fracx *     fracz *vpls_fab(ii,j,kk);
                    Real vmns_on_centroid = (1.0-fracx)*(1.0-fracz)*vmns_fab(i ,j,k )+
                                                 fracx *(1.0-fracz)*vmns_fab(ii,j,k )+
                                                 fracz *(1.0-fracx)*vmns_fab(i ,j,kk)+
                                                 fracx *     fracz *vmns_fab(ii,j,kk);
                    Real   ep_on_centroid = (1.0-fracx)*(1.0-fracz)* epy_fab(i ,j,k )+
                                                 fracx *(1.0-fracz)* epy_fab(ii,j,k )+
                                                 fracz *(1.0-fracx)* epy_fab(i ,j,kk)+
                                                 fracx *     fracz * epy_fab(ii,j,kk);

                    if ( vmns_on_centroid < 0.0 && vpls_on_centroid > 0.0 ) {
                       vmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( vpls_on_centroid + vmns_on_centroid );
                       if ( std::abs(avg) <  small_vel) { vmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { vmac_fab(i,j,k) = vmns_on_centroid;
                       } else                           { vmac_fab(i,j,k) = vpls_on_centroid;
                       }
                       vmac_fab(i,j,k) *= ep_on_centroid;
                    }

                 } else {

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

             AMREX_HOST_DEVICE_FOR_3D(wbx, i, j, k,
             {
                if (apz_fab(i,j,k) == 0.0)

                    wmac_fab(i,j,k) = huge_vel;

                else if (apz_fab(i,j,k) < 1.0) {

                    int ii = i + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,0)));
                    int jj = j + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,1)));

                    Real fracx = (ccm_fab(ii,j,k-1) || ccm_fab(ii,j,k)) ? std::abs(fcz_fab(i,j,k,0)) : 0.0;
                    Real fracy = (ccm_fab(i,jj,k-1) || ccm_fab(i,jj,k)) ? std::abs(fcz_fab(i,j,k,1)) : 0.0;

                    Real wpls_on_centroid = (1.0-fracx)*(1.0-fracy)*wpls_fab(i ,j ,k)+
                                                 fracx *(1.0-fracy)*wpls_fab(ii,j ,k)+
                                                 fracy *(1.0-fracx)*wpls_fab(i ,jj,k)+
                                                 fracx *     fracy *wpls_fab(ii,jj,k);
                    Real wmns_on_centroid = (1.0-fracx)*(1.0-fracy)*wmns_fab(i ,j ,k)+
                                                 fracx *(1.0-fracy)*wmns_fab(ii,j ,k)+
                                                 fracy *(1.0-fracx)*wmns_fab(i ,jj,k)+
                                                 fracx *     fracy *wmns_fab(ii,jj,k);
                    Real   ep_on_centroid = (1.0-fracx)*(1.0-fracy)* epz_fab(i ,j ,k)+
                                                 fracx *(1.0-fracy)* epz_fab(ii,j ,k)+
                                                 fracy *(1.0-fracx)* epz_fab(i ,jj,k)+
                                                 fracx *     fracy * epz_fab(ii,jj,k);

                    if ( wmns_on_centroid < 0.0 && wpls_on_centroid > 0.0 ) {
                       wmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( wpls_on_centroid + wmns_on_centroid );
                       if ( std::abs(avg) <  small_vel) { wmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { wmac_fab(i,j,k) = wmns_on_centroid;
                       } else                           { wmac_fab(i,j,k) = wpls_on_centroid;
                       }
                       wmac_fab(i,j,k) *= ep_on_centroid;
                    }

                 } else {

                    if ( wmns_fab(i,j,k) < 0.0 && wpls_fab(i,j,k) > 0.0 ) {
                       wmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( wpls_fab(i,j,k) + wmns_fab(i,j,k) );
                       if ( std::abs(avg) <  small_vel) { wmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { wmac_fab(i,j,k) = wmns_fab(i,j,k);
                       } else                           { wmac_fab(i,j,k) = wpls_fab(i,j,k);
                       }
                    }
                       wmac_fab(i,j,k) *= epz_fab(i,j,k);
                 }
             });

          } // Cut cells
       } // MFIter

    }   // end loop over levels

#ifdef AMREX_USE_CUDA
    Gpu::Device::synchronize();
#endif
}
