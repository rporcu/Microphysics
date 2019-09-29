#include <mfix.H>

void
mfix::mfix_compute_MAC_velocity_at_faces ( Real time,
                                           Vector< std::unique_ptr<MultiFab> >& vel_in,
                                           Vector< std::unique_ptr<MultiFab> >& ep_u_mac,
                                           Vector< std::unique_ptr<MultiFab> >& ep_v_mac,
                                           Vector< std::unique_ptr<MultiFab> >& ep_w_mac,
                                           Vector< std::unique_ptr<MultiFab> >& ep_in)

{
    BL_PROFILE("mfix::mfix_compute_MAC_velocity_at_faces");

    // ep_face is a temporary, no need to keep them outside this routine
    Vector< Array< std::unique_ptr<MultiFab>, AMREX_SPACEDIM> > ep_face;

    ep_face.resize(finest_level+1);

    for (int lev=0; lev < nlev; ++lev)
    {
       Box domain(geom[lev].Domain());   

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

             // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
             AMREX_HOST_DEVICE_FOR_3D(ubx, i, j, k, 
             {
                 // X-faces
                 Real upls     = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                 Real umns     = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                 if ( umns < 0.0 && upls > 0.0 ) {
                    umac_fab(i,j,k) = 0.0;
                 } else {
                    Real avg = 0.5 * ( upls + umns );
                    if ( std::abs(avg) <  small_vel) { umac_fab(i,j,k) = 0.0;
                    } else if (avg >= 0)             { umac_fab(i,j,k) = umns;
                    } else                           { umac_fab(i,j,k) = upls;
                    }
                 }
             });

             AMREX_HOST_DEVICE_FOR_3D(vbx, i, j, k,
             {
                 // Y-faces
                 Real upls     = ccvel_fab(i,j  ,k,1) - 0.5 * yslopes_fab(i,j  ,k,1);
                 Real umns     = ccvel_fab(i,j-1,k,1) + 0.5 * yslopes_fab(i,j-1,k,1);
                 if ( umns < 0.0 && upls > 0.0 ) {
                    vmac_fab(i,j,k) = 0.0;
                 } else {
                    Real avg = 0.5 * ( upls + umns );
                    if ( std::abs(avg) <  small_vel) { vmac_fab(i,j,k) = 0.0;
                    } else if (avg >= 0)             { vmac_fab(i,j,k) = umns;
                    } else                           { vmac_fab(i,j,k) = upls;
                    }
                 }
             });

             AMREX_HOST_DEVICE_FOR_3D(wbx, i, j, k,
             {
                 // Z-faces
                 Real upls     = ccvel_fab(i,j,k  ,2) - 0.5 * zslopes_fab(i,j,k  ,2);
                 Real umns     = ccvel_fab(i,j,k-1,2) + 0.5 * zslopes_fab(i,j,k-1,2);
                 if ( umns < 0.0 && upls > 0.0 ) {
                    wmac_fab(i,j,k) = 0.0;
                 } else {
                    Real avg = 0.5 * ( upls + umns );
                    if ( std::abs(avg) <  small_vel) { wmac_fab(i,j,k) = 0.0;
                    } else if (avg >= 0)             { wmac_fab(i,j,k) = umns;
                    } else                           { wmac_fab(i,j,k) = upls;
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

             // Face-centered velocity components
             const auto& umac_fab = (ep_u_mac[lev])->array(mfi);
             const auto& vmac_fab = (ep_v_mac[lev])->array(mfi);
             const auto& wmac_fab = (ep_w_mac[lev])->array(mfi);

             // Face-centered areas
             const auto& ax_fab = areafrac[0]->array(mfi);
             const auto& ay_fab = areafrac[1]->array(mfi);
             const auto& az_fab = areafrac[2]->array(mfi);

             // This FAB has cut cells
             AMREX_HOST_DEVICE_FOR_3D(ubx, i, j, k, 
             {
                 // X-faces
                 if (ax_fab(i,j,k) > 0.0)
                 {
                    Real upls     = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                    Real umns     = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                    if ( umns < 0.0 && upls > 0.0 ) {
                       umac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( upls + umns );
                       if ( std::abs(avg) <  small_vel) { umac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { umac_fab(i,j,k) = umns;
                       } else                           { umac_fab(i,j,k) = upls;
                       }
                    }
                 } else {
                       umac_fab(i,j,k) = huge_vel;
                 }
             });

             AMREX_HOST_DEVICE_FOR_3D(vbx, i, j, k,
             {
                 // Y-faces
                 if (ay_fab(i,j,k) > 0.0)
                 {
                    Real upls     = ccvel_fab(i,j  ,k,1) - 0.5 * yslopes_fab(i,j  ,k,1);
                    Real umns     = ccvel_fab(i,j-1,k,1) + 0.5 * yslopes_fab(i,j-1,k,1);
                    if ( umns < 0.0 && upls > 0.0 ) {
                       vmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( upls + umns );
                       if ( std::abs(avg) <  small_vel) { vmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { vmac_fab(i,j,k) = umns;
                       } else                           { vmac_fab(i,j,k) = upls;
                       }
                    }
                 } else {
                       vmac_fab(i,j,k) = huge_vel;
                 }
             });

             AMREX_HOST_DEVICE_FOR_3D(wbx, i, j, k,
             {
                 // Z-faces
                 if (az_fab(i,j,k) > 0.0)
                 {
                    Real upls     = ccvel_fab(i,j,k  ,2) - 0.5 * zslopes_fab(i,j,k  ,2);
                    Real umns     = ccvel_fab(i,j,k-1,2) + 0.5 * zslopes_fab(i,j,k-1,2);
                    if ( umns < 0.0 && upls > 0.0 ) {
                       wmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( upls + umns );
                       if ( std::abs(avg) <  small_vel) { wmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { wmac_fab(i,j,k) = umns;
                       } else                           { wmac_fab(i,j,k) = upls;
                       }
                    }
                 } else {
                       wmac_fab(i,j,k) = huge_vel;
                 }
             });

          } // Cut cells
       } // MFIter

      ep_in[lev]->FillBoundary(geom[lev].periodicity());

      BoxArray x_edge_ba = grids[lev];
      x_edge_ba.surroundingNodes(0);
      BoxArray y_edge_ba = grids[lev];
      y_edge_ba.surroundingNodes(1);
      BoxArray z_edge_ba = grids[lev];
      z_edge_ba.surroundingNodes(2);

      ep_face[lev][0].reset(new MultiFab(x_edge_ba,dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      ep_face[lev][1].reset(new MultiFab(y_edge_ba,dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      ep_face[lev][2].reset(new MultiFab(z_edge_ba,dmap[lev],1,0,MFInfo(),*ebfactory[lev]));

      average_cellcenter_to_face( GetArrOfPtrs(ep_face[lev]), *ep_in[lev], geom[lev] );

      // Compute ep*u at faces and store it in ep_u_mac, ep_v_mac, ep_w_mac
      MultiFab::Multiply( *ep_u_mac[lev], *(ep_face[lev][0]), 0, 0, 1, 0 );
      MultiFab::Multiply( *ep_v_mac[lev], *(ep_face[lev][1]), 0, 0, 1, 0 );
      MultiFab::Multiply( *ep_w_mac[lev], *(ep_face[lev][2]), 0, 0, 1, 0 );

      // Set velocity bcs -- after we multiply by ep
      set_MAC_velocity_bcs( lev, ep_u_mac, ep_v_mac, ep_w_mac, time );
    }

#ifdef AMREX_USE_CUDA
    Gpu::Device::synchronize();
#endif

    // Do projection on all AMR levels in one shot 
    apply_MAC_projection (ep_u_mac, ep_v_mac, ep_w_mac, ep_g, ro_g, time, steady_state );
}
