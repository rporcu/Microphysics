#include <mfix.H>

void
mfix::mfix_compute_MAC_velocity_at_faces ( Real time,
                                           Vector< std::unique_ptr<MultiFab> >& vel)

{
    BL_PROFILE("mfix::mfix_compute_MAC_velocity_at_faces");

    for (int lev=0; lev < nlev; ++lev)
    {
       Box domain(geom[lev].Domain());   

       // State with ghost cells
       MultiFab Sborder(grids[lev], dmap[lev], vel_g[lev]->nComp(), nghost,  MFInfo(), *ebfactory[lev]);
       FillPatchVel(lev, time, Sborder, 0, Sborder.nComp(), bcs_u);
    
       // First compute the slopes
       mfix_compute_velocity_slopes( lev, time, Sborder );

       // Copy each FAB back from Sborder into the vel array, complete with filled ghost cells
       MultiFab::Copy (*vel[lev],  Sborder,  0, 0,  vel[lev]->nComp(),  vel[lev]->nGrow());

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
       for (MFIter mfi(Sborder,TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
           // Tilebox
          Box  bx = mfi.tilebox();
          Box ubx = mfi.tilebox(e_x);
          Box vbx = mfi.tilebox(e_y);
          Box wbx = mfi.tilebox(e_z);
      
          // Check efficiently if this tile contains any eb stuff

          const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_g[lev])[mfi]);
          const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

          if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
          {
             m_u_mac[lev] -> setVal( 1.2345e300, ubx, 0, 1);
             m_v_mac[lev] -> setVal( 1.2345e300, vbx, 0, 1);
             m_w_mac[lev] -> setVal( 1.2345e300, wbx, 0, 1);
          }
          else if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
          {

             // Cell-centered velocity
             const auto& ccvel_fab = vel_g[lev]->array(mfi);

             // Cell-centered slopes
             const auto& xslopes_fab = (xslopes[lev])->array(mfi);
             const auto& yslopes_fab = (yslopes[lev])->array(mfi);
             const auto& zslopes_fab = (zslopes[lev])->array(mfi);

             // Face-centered velocity components
             const auto& umac_fab = (m_u_mac[lev])->array(mfi);
             const auto& vmac_fab = (m_v_mac[lev])->array(mfi);
             const auto& wmac_fab = (m_w_mac[lev])->array(mfi);

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
             const auto& ccvel_fab = vel_g[lev]->array(mfi);

             // Cell-centered slopes
             const auto& xslopes_fab = (xslopes[lev])->array(mfi);
             const auto& yslopes_fab = (yslopes[lev])->array(mfi);
             const auto& zslopes_fab = (zslopes[lev])->array(mfi);

             // Face-centered velocity components
             const auto& umac_fab = (m_u_mac[lev])->array(mfi);
             const auto& vmac_fab = (m_v_mac[lev])->array(mfi);
             const auto& wmac_fab = (m_w_mac[lev])->array(mfi);

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
    }   // end loop over levels

#ifdef AMREX_USE_CUDA
    Gpu::Device::synchronize();
#endif

    // Note that we will call set_velocity_bcs in mac_projection so we don't need to call it here
    // set_velocity_bcs( lev, m_u_mac, m_v_mac, m_w_mac, time );

    // Do projection on all AMR levels in one shot 
    apply_MAC_projection (m_u_mac, m_v_mac, m_w_mac, ep_g, ro_g, time, steady_state );
}
