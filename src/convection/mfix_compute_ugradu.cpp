#include <mfix.H>
#include <mfix_F.H>
#include <mfix_proj_F.H>
#include <mfix_mac_F.H>

//
// Compute acc using the vel passed in
//
void
mfix::mfix_compute_ugradu_predictor( Vector< std::unique_ptr<MultiFab> >& conv, 
                                     Vector< std::unique_ptr<MultiFab> >& vel,
	          		     Real time)
{
    BL_PROFILE("mfix::mfix_compute_ugradu_predictor");

    amrex::Print() << "In predictor at time " << time << std::endl;

    mfix_compute_MAC_velocity_at_faces( time, vel );

    for (int lev=0; lev < nlev; ++lev)
    {
        Box domain(geom[lev].Domain());
    
        // Get EB geometric info
        Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
        Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;
        const amrex::MultiFab*                    volfrac;
        const amrex::MultiCutFab*                 bndrycent;

        areafrac  =   ebfactory[lev] -> getAreaFrac();
        facecent  =   ebfactory[lev] -> getFaceCent();
        volfrac   = &(ebfactory[lev] -> getVolFrac());
        bndrycent = &(ebfactory[lev] -> getBndryCent());
       
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*vel[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox ();
            
            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
            const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

            if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
            {
                // If tile is completely covered by EB geometry, set slopes
                // value to some very large number so we know if
                // we accidentaly use these covered slopes later in calculations
                conv[lev] -> setVal( 1.2345e300, bx, 0, 3);
            }
            else
            {
                // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
                if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
                {
                    compute_ugradu(
                        BL_TO_FORTRAN_BOX(bx),  
                        BL_TO_FORTRAN_ANYD((*conv[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((    *vel[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((   *ep_g[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
                        (*xslopes[lev])[mfi].dataPtr(),
                        (*yslopes[lev])[mfi].dataPtr(),
                        BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
                        domain.loVect (), domain.hiVect (),
                        bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                        bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                        bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                        geom[lev].CellSize(), &nghost);
                }
                else
                {
                    compute_ugradu_eb(
                        BL_TO_FORTRAN_BOX(bx),  
                        BL_TO_FORTRAN_ANYD((*conv[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((    *vel[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((   *ep_g[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
                        BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
                        BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
                        BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
                        BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
                        BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
                        BL_TO_FORTRAN_ANYD(flags),
                        BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                        BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
                        (*xslopes[lev])[mfi].dataPtr(),
                        (*yslopes[lev])[mfi].dataPtr(),
                        BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
                        domain.loVect (), domain.hiVect (),
                        bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                        bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                        bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                        geom[lev].CellSize(), &nghost);
                }
            }
        }

    }
}

//
// Compute acc using the vel passed in
//
void
mfix::mfix_compute_ugradu_corrector( Vector< std::unique_ptr<MultiFab> >& conv, 
                                     Vector< std::unique_ptr<MultiFab> >& vel,
	          		   Real time)
{
    BL_PROFILE("mfix::mfix_compute_ugradu_corrector");

    amrex::Print() << "In corrector at time " << time << std::endl;

    mfix_compute_MAC_velocity_at_faces( time, vel );

    for (int lev=0; lev < nlev; ++lev)
    {
        Box domain(geom[lev].Domain());
    
        // Get EB geometric info
        Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
        Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;
        const amrex::MultiFab*                    volfrac;
        const amrex::MultiCutFab*                 bndrycent;

        areafrac  =   ebfactory[lev] -> getAreaFrac();
        facecent  =   ebfactory[lev] -> getFaceCent();
        volfrac   = &(ebfactory[lev] -> getVolFrac());
        bndrycent = &(ebfactory[lev] -> getBndryCent());
       
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*vel[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox ();
            
            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
            const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

            if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
            {
                // If tile is completely covered by EB geometry, set slopes
                // value to some very large number so we know if
                // we accidentaly use these covered slopes later in calculations
                conv[lev] -> setVal( 1.2345e300, bx, 0, 3);
            }
            else
            {
                // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
                if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
                {
                    compute_ugradu(
                        BL_TO_FORTRAN_BOX(bx),  
                        BL_TO_FORTRAN_ANYD((*conv[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((    *vel[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((   *ep_g[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
                        (*xslopes[lev])[mfi].dataPtr(),
                        (*yslopes[lev])[mfi].dataPtr(),
                        BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
                        domain.loVect (), domain.hiVect (),
                        bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                        bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                        bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                        geom[lev].CellSize(), &nghost);
                }
                else
                {
                    compute_ugradu_eb(
                        BL_TO_FORTRAN_BOX(bx),  
                        BL_TO_FORTRAN_ANYD((*conv[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((    *vel[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((   *ep_g[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
                        BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
                        BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
                        BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
                        BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
                        BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
                        BL_TO_FORTRAN_ANYD(flags),
                        BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                        BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
                        (*xslopes[lev])[mfi].dataPtr(),
                        (*yslopes[lev])[mfi].dataPtr(),
                        BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
                        domain.loVect (), domain.hiVect (),
                        bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                        bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                        bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                        geom[lev].CellSize(), &nghost);

                }
            }
        }

    }
}
