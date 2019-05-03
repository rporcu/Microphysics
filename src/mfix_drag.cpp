#include <AMReX_ParmParse.H>
#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix.H>
#include <mfix_des_F.H>
#include <mfix_drag_K.H>
#include <mfix_util_F.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>

void mfix::mfix_calc_drag_fluid(Real time)
{
    BL_PROFILE("mfix::mfix_calc_drag_fluid()");

    mfix_calc_particle_beta();

    // ******************************************************************************
    // Now use the beta of individual particles to create the drag terms on the fluid
    // ******************************************************************************
    for (int lev = 0; lev < nlev; lev++)
       drag[lev] ->setVal(0.0L);

    pc -> CalcDragOnFluid(drag, particle_ebfactory,
                          bc_ilo,bc_ihi,bc_jlo,bc_jhi,bc_klo,bc_khi,
                          nghost);

    // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
    for (int lev = 0; lev < nlev; lev++)
        drag[lev] -> FillBoundary(geom[lev].periodicity());
}

void
mfix::mfix_calc_drag_particle(Real time)
{
    BL_PROFILE("mfix::mfix_calc_drag_particle()");

    for (int lev = 0; lev < nlev; lev++)
    {

        Box      domain(geom[lev].Domain());
        MultiFab gp_tmp;

        gp_tmp.define(grids[lev],dmap[lev],3,1);

        MultiFab::Copy(gp_tmp, *gp[lev], 0, 0, 3, 1);
        gp_tmp.FillBoundary(geom[lev].periodicity());

        //
        // NOTE -- it is essential that we call set_gradp_bcs after calling FillBoundary
        //         because the set_gradp_bcs call hopefully sets the ghost cells exterior
        //         to the domain from ghost cells interior to the domain
        //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(gp_tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            set_gradp_bcs ( bx.loVect(), bx.hiVect(),
                            BL_TO_FORTRAN_ANYD(gp_tmp[mfi]),
                            bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                            bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                            bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                            domain.loVect(), domain.hiVect(),
                            &nghost);
        }

        // Extrapolate velocity Dirichlet bc's to ghost cells
        // HACK -- NOTE WE ARE CALLING THIS ON ALL LEVELS BUT ONLY NEED IT ON ONE LEVEL
        int extrap_dir_bcs = 1;
        mfix_set_velocity_bcs(time, extrap_dir_bcs);
        gp_tmp.FillBoundary(geom[lev].periodicity());

        bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                             (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );


        // Pointers to Multifabs of interest
        MultiFab*  gp_ptr;
        MultiFab* vel_ptr;

        // Temporararies for dual grid case
        std::unique_ptr<MultiFab>   gp_pba;
        std::unique_ptr<MultiFab>  vel_pba;

        if (OnSameGrids)
        {
            gp_ptr  = &gp_tmp;
            vel_ptr = vel_g[lev].get();
        }
        else
        {
            BoxArray            pba = pc->ParticleBoxArray(lev);
            DistributionMapping pdm = pc->ParticleDistributionMap(lev);

            int ng = gp_tmp.nGrow();
            gp_pba.reset(new MultiFab(pba,pdm,gp_tmp.nComp(),ng));
            gp_pba->copy(gp_tmp,0,0,gp_tmp.nComp(),ng,ng);
            gp_pba->FillBoundary(geom[lev].periodicity());

            EBFArrayBoxFactory ebfactory_loc( * eb_levels[lev], geom[lev], pba, pdm,
                                              {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                               m_eb_full_grow_cells}, EBSupport::basic);

            ng = vel_g[lev]->nGrow();
            vel_pba.reset(new MultiFab(pba,pdm,vel_g[lev]->nComp(),ng,MFInfo(), ebfactory_loc));
            vel_pba->copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),ng,ng);
            vel_pba->FillBoundary(geom[lev].periodicity());

            gp_ptr  = gp_pba.get();
            vel_ptr = vel_pba.get();
        }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            Real velfp[3];
            Real gradp[3];

            const auto dxi = geom[lev].InvCellSizeArray();
            const auto dx  = geom[lev].CellSizeArray();
            const auto plo = geom[lev].ProbLoArray();

            for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
            {
                auto& particles = pti.GetArrayOfStructs();
                const int np = particles.size();
                Box bx = pti.tilebox ();

                // this is to check efficiently if this tile contains any eb stuff
                const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_ptr)[pti]);
                const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

                if (flags.getType(amrex::grow(bx,0)) != FabType::covered)
                {
                    const auto& vel_array = vel_ptr->array(pti);
                    const auto&  gp_array =  gp_ptr->array(pti);

                    const auto& flags_array = flags.array();

                    if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
                    {

                        for (int ip = 0; ip < np; ++ip)
                        {
                           MFIXParticleContainer::ParticleType& particle = particles[ip];
                           Real pbeta = particle.rdata(realData::dragx);

                           trilinear_interp(particle, &velfp[0], vel_array, plo, dxi);
                           trilinear_interp(particle, &gradp[0],  gp_array, plo, dxi);

                           // Particle drag calculation
                           particle.rdata(realData::dragx) = pbeta * ( velfp[0] - particle.rdata(realData::velx) ) -
                                                            (gradp[0] + gp0[0]) * particle.rdata(realData::volume);

                           particle.rdata(realData::dragy) = pbeta * ( velfp[1] - particle.rdata(realData::vely) ) -
                                                            (gradp[1] + gp0[1]) * particle.rdata(realData::volume);

                           particle.rdata(realData::dragz) = pbeta * ( velfp[2] - particle.rdata(realData::velz) ) -
                                                            (gradp[2] + gp0[2]) * particle.rdata(realData::volume);
                        }
                    }
                    else // FAB not all regular
                    {
                        // Phi is always on the particles grid -- we do not use the refined level here
                        const MultiFab & phi = * level_sets[lev];
                        const auto& phi_array = phi.array(pti);

                        const MultiCutFab* bndrycent = &(ebfactory[lev] -> getBndryCent());
                        const auto& bct_array = bndrycent->array(pti);

                        const auto&  gp_array   =  gp_ptr->array(pti);

                        for (int ip = 0; ip < np; ++ip)
                        {
                            MFIXParticleContainer::ParticleType& particle = particles[ip];
                            Real pbeta = particle.rdata(realData::dragx);

                            // This identifies which cell the particle is in
                            int iloc = std::floor((particle.pos(0) - plo[0])*dxi[0]);
                            int jloc = std::floor((particle.pos(1) - plo[1])*dxi[1]);
                            int kloc = std::floor((particle.pos(2) - plo[2])*dxi[2]);

                            // Pick upper cell in the stencil
                            Real lx = (particle.pos(0) - plo[0])*dxi[0] + 0.5;
                            Real ly = (particle.pos(1) - plo[1])*dxi[1] + 0.5;
                            Real lz = (particle.pos(2) - plo[2])*dxi[2] + 0.5;

                            int i = std::floor(lx);
                            int j = std::floor(ly);
                            int k = std::floor(lz);

                            // If the particle has gone outside the fluid region then we will
                            //    let the wall collision term bring it back and not use any
                            //    extrapolated fluid quantities in a covered region
                            if (flags_array(iloc,jloc,kloc).isCovered())
                            {
                                particle.rdata(realData::dragx) = 0.0;
                                particle.rdata(realData::dragy) = 0.0;
                                particle.rdata(realData::dragz) = 0.0;
                            }

                            // Regular cell
                            else if (flags_array(iloc,jloc,kloc).isRegular()) 
                            {

                                // This cell is regular which means its not touching a covered cell
                                trilinear_interp(particle, &velfp[0], vel_array, plo, dxi);
                                trilinear_interp(particle, &gradp[0],  gp_array, plo, dxi);

                                particle.rdata(realData::dragx) = pbeta * ( velfp[0] - particle.rdata(realData::velx) ) -
                                                                 (gradp[0] + gp0[0]) * particle.rdata(realData::volume);
 
                                particle.rdata(realData::dragy) = pbeta * ( velfp[1] - particle.rdata(realData::vely) ) -
                                                                 (gradp[1] + gp0[1]) * particle.rdata(realData::volume);
 
                                particle.rdata(realData::dragz) = pbeta * ( velfp[2] - particle.rdata(realData::velz) ) -
                                                                 (gradp[2] + gp0[2]) * particle.rdata(realData::volume);
                            }

                            // Cut cell
                            else 
                            {
                                if (!flags_array(i-1,j-1,k-1).isCovered() &&
                                    !flags_array(i  ,j-1,k-1).isCovered() &&
                                    !flags_array(i-1,j  ,k-1).isCovered() &&
                                    !flags_array(i  ,j  ,k-1).isCovered() &&
                                    !flags_array(i-1,j-1,k  ).isCovered() &&
                                    !flags_array(i  ,j-1,k  ).isCovered() &&
                                    !flags_array(i-1,j  ,k  ).isCovered() &&
                                    !flags_array(i  ,j  ,k  ).isCovered()) 
                                {
                                    // None of the cells in the stencil is covered; we can use the regular formula
                                    trilinear_interp(particle, &velfp[0], vel_array, plo, dxi);
                                    trilinear_interp(particle, &gradp[0],  gp_array, plo, dxi);
    
                               } else {

                                    // One of the cells in the stencil is covered

                                    Real centroid_pos[3];
                                    centroid_pos[0] = bct_array(iloc,jloc,kloc,0);
                                    centroid_pos[1] = bct_array(iloc,jloc,kloc,1);
                                    centroid_pos[2] = bct_array(iloc,jloc,kloc,2);

                                    Real anrm[3];
                                    Real  pos[3];
                                    pos[0] = plo[0] + (iloc + 0.5 + centroid_pos[0])*dx[0];
                                    pos[1] = plo[1] + (jloc + 0.5 + centroid_pos[1])*dx[1];
                                    pos[2] = plo[2] + (kloc + 0.5 + centroid_pos[2])*dx[2];

                                    // We find the normal at the face centroid
                                    normal_from_ls(anrm, pos, phi_array, plo, dxi);

                                    // Scaled particle position relative to cell center (same as bndryCentroid scaling)
                                    Real px = particle.pos(0)*dxi[0]-(iloc+.5);
                                    Real py = particle.pos(1)*dxi[0]-(jloc+.5);
                                    Real pz = particle.pos(2)*dxi[0]-(kloc+.5);
    
                                    // Distance from plane:  (particle pos - centroid pos) dot (normal)
                                    Real dist = (centroid_pos[0] - px) * anrm[0] + 
                                                (centroid_pos[1] - py) * anrm[1] + 
                                                (centroid_pos[2] - pz) * anrm[2];
        
                                    // NOTE THIS ALGORITHM INTERPOLATES TO POINT BETWEEN CELL CENTER AND BNDRY CENTROID
                                    //      THAT IS SAME DISTANCE FROM WALL AS PARTICLE IS -- IT DOES NOT INTERPOLATE
                                    //      TO CORRECT PARTICLE LOCATION
                                    // Distance from cell center (iloc,jloc,kloc)
                                    Real gx = centroid_pos[0] - dist*anrm[0];
                                    Real gy = centroid_pos[1] - dist*anrm[1];
                                    Real gz = centroid_pos[2] - dist*anrm[2];
    
                                    int ii,jj,kk;

                                    if (anrm[0] < 0) {
                                        ii = iloc - 1;
                                    } else {
                                        ii = iloc + 1; 
                                        gx = -gx;
                                    }
                                    if (anrm[1] < 0) {
                                        jj = jloc - 1;
                                    } else {
                                        jj = jloc + 1; 
                                        gy = -gy;
                                    }
                                    if (anrm[2] < 0) {
                                        kk = kloc - 1;
                                    } else {
                                        kk = kloc + 1; 
                                        gz = -gz;
                                    }
    
                                    Real gxy = gx*gy;
                                    Real gxz = gx*gz;
                                    Real gyz = gy*gz;
                                    Real gxyz = gx*gy*gz;
    
                                    for (int n = 0; n < 3; n++)
                                    {
                                       velfp[n] = (1.0+gx+gy+gz+gxy+gxz+gyz+gxyz) * vel_array(iloc,jloc,kloc,n)
                                                + (-gz - gxz - gyz - gxyz)        * vel_array(iloc,jloc,kk  ,n)
                                                + (-gy - gxy - gyz - gxyz)        * vel_array(iloc,jj  ,kloc,n)
                                                + (gyz + gxyz)                    * vel_array(iloc,jj  ,kk ,n)
                                                + (-gx - gxy - gxz - gxyz)        * vel_array(ii  ,jloc,kloc,n)
                                                + (gxz + gxyz)                    * vel_array(ii  ,jloc,kk ,n)
                                                + (gxy + gxyz)                    * vel_array(ii  ,jj  ,kloc,n)
                                                + (-gxyz)                         * vel_array(ii  ,jj  ,kk ,n);
                                       gradp[n] = (1.0+gx+gy+gz+gxy+gxz+gyz+gxyz) *  gp_array(iloc,jloc,kloc,n)
                                                + (-gz - gxz - gyz - gxyz)        *  gp_array(iloc,jloc,kk  ,n)
                                                + (-gy - gxy - gyz - gxyz)        *  gp_array(iloc,jj  ,kloc,n)
                                                + (gyz + gxyz)                    *  gp_array(iloc,jj  ,kk ,n)
                                                + (-gx - gxy - gxz - gxyz)        *  gp_array(ii  ,jloc,kloc,n)
                                                + (gxz + gxyz)                    *  gp_array(ii  ,jloc,kk ,n)
                                                + (gxy + gxyz)                    *  gp_array(ii  ,jj  ,kloc,n)
                                                + (-gxyz)                         *  gp_array(ii  ,jj  ,kk ,n);
                                    }
                               }

                               particle.rdata(realData::dragx) = pbeta * ( velfp[0] - particle.rdata(realData::velx) ) -
                                                                (gradp[0] + gp0[0]) * particle.rdata(realData::volume);

                               particle.rdata(realData::dragy) = pbeta * ( velfp[1] - particle.rdata(realData::vely) ) -
                                                                (gradp[1] + gp0[1]) * particle.rdata(realData::volume);

                               particle.rdata(realData::dragz) = pbeta * ( velfp[2] - particle.rdata(realData::velz) ) -
                                                                (gradp[2] + gp0[2]) * particle.rdata(realData::volume);

                            } // Test on type of cell
                        } // particle loop
                } // if box not all regular

            } // FAB not covered
        } // pti

        // Reset velocity Dirichlet bc's to face values
        // HACK -- NOTE WE ARE CALLING THIS ON ALL LEVELS BUT ONLY NEED IT ON ONE LEVEL
        extrap_dir_bcs = 0;
        mfix_set_velocity_bcs(time, extrap_dir_bcs);

        } // omp region

    } // lev
}
