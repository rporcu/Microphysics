#include <mfix.H>
#include <mfix_F.H>
#include <mfix_proj_F.H>
#include <mfix_mac_F.H>

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

          if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
          {
             m_u_mac[lev] -> setVal( 1.2345e300, ubx, 0, 1);
             m_v_mac[lev] -> setVal( 1.2345e300, vbx, 0, 1);
             m_w_mac[lev] -> setVal( 1.2345e300, wbx, 0, 1);
          }
          else
          {
             // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
             if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
             {
                AMREX_CUDA_HOST_DEVICE_FOR_3D(ubx, i, j, k, 
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

                AMREX_CUDA_HOST_DEVICE_FOR_3D(vbx, i, j, k,
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

                AMREX_CUDA_HOST_DEVICE_FOR_3D(wbx, i, j, k,
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

#if 0
                compute_velocity_at_faces(
                   BL_TO_FORTRAN_BOX(bx),  
                   BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((    *vel[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((*xslopes[lev])[mfi]),
                   (*yslopes[lev])[mfi].dataPtr(),
                   (*zslopes[lev])[mfi].dataPtr(),
                   bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                   bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                   bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                   &nghost, domain.loVect(), domain.hiVect() );
#endif
             }
             else
             {
                compute_velocity_at_x_faces_eb(
                   BL_TO_FORTRAN_BOX(ubx),
                   BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((    *vel[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((*xslopes[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
                   BL_TO_FORTRAN_ANYD(flags),
                   bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                   &nghost, domain.loVect (), domain.hiVect () );           

                compute_velocity_at_y_faces_eb(
                   BL_TO_FORTRAN_BOX(vbx),
                   BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((    *vel[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((*yslopes[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
                   BL_TO_FORTRAN_ANYD(flags),
                   bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                   &nghost, domain.loVect (), domain.hiVect () );

                compute_velocity_at_z_faces_eb(
                   BL_TO_FORTRAN_BOX(wbx),
                   BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((    *vel[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
                   BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
                   BL_TO_FORTRAN_ANYD(flags),
                   bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                   &nghost, domain.loVect (), domain.hiVect () );
             }
          }
       }
    }   // end loop over levels

    // Note that we will call set_velocity_bcs in mac_projection so we don't need to call it here
    // set_velocity_bcs( lev, m_u_mac, m_v_mac, m_w_mac, time );

    // Do projection on all AMR levels in one shot 
    mac_projection -> apply_projection (m_u_mac, m_v_mac, m_w_mac, ep_g, ro_g, time, steady_state );
}
