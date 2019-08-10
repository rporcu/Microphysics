#include <mfix_F.H>
#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_drag_K.H>
#include <AMReX_EBMultiFabUtil.H>

void mfix::mfix_calc_drag_fluid(Real time)
{
    mfix_calc_particle_beta(time);

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


    // Extrapolate velocity Dirichlet bc's to ghost cells
    int extrap_dir_bcs = 1;
    mfix_set_velocity_bcs(time, vel_g, extrap_dir_bcs);

    for (int lev = 0; lev < nlev; lev++)
    {

        Box      domain(geom[lev].Domain());
        MultiFab gp_tmp;

        gp_tmp.define(grids[lev],dmap[lev],3,1,MFInfo(),*ebfactory[lev]);

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

            set_gradp_bcs(bx, lev, gp_tmp[mfi], domain);
        }

        gp_tmp.FillBoundary(geom[lev].periodicity());

        bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                             (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

        // Pointers to Multifabs of interest
        MultiFab*  gp_ptr;
        MultiFab* vel_ptr;

        // Temporararies for dual grid case
        std::unique_ptr<MultiFab>   gp_pba;
        std::unique_ptr<MultiFab>  vel_pba;

        // This is just a sanity check to make sure we're not using covered values
        // We can remove these lines once we're confident in the algoirthm
        EB_set_covered(*vel_g[0], 0, 3, 1, covered_val);
        EB_set_covered( gp_tmp  , 0, 3, 1, covered_val);

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
            //const auto dx  = geom[lev].CellSizeArray(); // SET_BUT_NOT_USED
            const auto plo = geom[lev].ProbLoArray();

            for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
            {
                auto& particles = pti.GetArrayOfStructs();
                const int np = particles.size();

                Box bx = pti.tilebox ();

                // This is to check efficiently if this tile contains any eb stuff
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

                           trilinear_interp(particle, &velfp[0], vel_array, plo, dxi);
                           trilinear_interp(particle, &gradp[0],  gp_array, plo, dxi);

                           Real pbeta = particle.rdata(realData::dragx);

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
                        // phi is always on the particles grid -- we do not use the refined level here
                        const MultiFab & phi = * level_sets[lev];
                        const auto& phi_array = phi.array(pti);

                        for (int ip = 0; ip < np; ++ip)
                        {
                            MFIXParticleContainer::ParticleType& particle = particles[ip];
                            Real pbeta = particle.rdata(realData::dragx);

                            // This identifies which cell the particle is in
                            int iloc = floor((particle.pos(0) - plo[0])*dxi[0]);
                            int jloc = floor((particle.pos(1) - plo[1])*dxi[1]);
                            int kloc = floor((particle.pos(2) - plo[2])*dxi[2]);

                            // Pick upper cell in the stencil
                            Real lx = (particle.pos(0) - plo[0])*dxi[0] + 0.5;
                            Real ly = (particle.pos(1) - plo[1])*dxi[1] + 0.5;
                            Real lz = (particle.pos(2) - plo[2])*dxi[2] + 0.5;

                            int i = std::floor(lx);
                            int j = std::floor(ly);
                            int k = std::floor(lz);

                            // Covered cell
                            if (flags_array(iloc,jloc,kloc).isCovered())
                            {
                                particle.rdata(realData::dragx) = 0.0;
                                particle.rdata(realData::dragy) = 0.0;
                                particle.rdata(realData::dragz) = 0.0;

                            } else {

                                // Cut or regular cell and none of the cells in the stencil is covered
                                // (Note we can't assume regular cell has no covered cells in the stencil
                                //      because of the diagonal case)
                                if (!flags_array(i-1,j-1,k-1).isCovered() &&
                                    !flags_array(i  ,j-1,k-1).isCovered() &&
                                    !flags_array(i-1,j  ,k-1).isCovered() &&
                                    !flags_array(i  ,j  ,k-1).isCovered() &&
                                    !flags_array(i-1,j-1,k  ).isCovered() &&
                                    !flags_array(i  ,j-1,k  ).isCovered() &&
                                    !flags_array(i-1,j  ,k  ).isCovered() &&
                                    !flags_array(i  ,j  ,k  ).isCovered())
                                {
                                    trilinear_interp(particle, &velfp[0], vel_array, plo, dxi);
                                    trilinear_interp(particle, &gradp[0],  gp_array, plo, dxi);

                                // At least one of the cells in the stencil is covered
                                } else {

                                    Real anrm[3];

                                    // Compute distance of the particle from the wall.
                                    // (This is the same function we call when computing the particle-wall collisions)
                                    int ls_refinement = 1;
                                    //Real dist = interp_level_set(particle, ls_refinement, phi_array, plo, dxi); // UNUSED_VARIABLE

                                    // Compute the normal to the wall in this cell -- it doesn't matter
                                    // whether we compute it "at the particle location" or "at the centroid location"
                                    // because it interpolates from the same values of phi.
                                    level_set_normal(particle, ls_refinement, &anrm[0], phi_array, plo, dxi);

                                    // Particle position must be in [-.5:.5] is relative to cell center and scaled by dx
                                    Real gx = particle.pos(0)*dxi[0] - (iloc + 0.5);
                                    Real gy = particle.pos(1)*dxi[1] - (jloc + 0.5);
                                    Real gz = particle.pos(2)*dxi[2] - (kloc + 0.5);

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
                                                + (gyz + gxyz)                    * vel_array(iloc,jj  ,kk  ,n)
                                                + (-gx - gxy - gxz - gxyz)        * vel_array(ii  ,jloc,kloc,n)
                                                + (gxz + gxyz)                    * vel_array(ii  ,jloc,kk  ,n)
                                                + (gxy + gxyz)                    * vel_array(ii  ,jj  ,kloc,n)
                                                + (-gxyz)                         * vel_array(ii  ,jj  ,kk  ,n);
                                       gradp[n] = (1.0+gx+gy+gz+gxy+gxz+gyz+gxyz) *  gp_array(iloc,jloc,kloc,n)
                                                + (-gz - gxz - gyz - gxyz)        *  gp_array(iloc,jloc,kk  ,n)
                                                + (-gy - gxy - gyz - gxyz)        *  gp_array(iloc,jj  ,kloc,n)
                                                + (gyz + gxyz)                    *  gp_array(iloc,jj  ,kk ,n)
                                                + (-gx - gxy - gxz - gxyz)        *  gp_array(ii  ,jloc,kloc,n)
                                                + (gxz + gxyz)                    *  gp_array(ii  ,jloc,kk ,n)
                                                + (gxy + gxyz)                    *  gp_array(ii  ,jj  ,kloc,n)
                                                + (-gxyz)                         *  gp_array(ii  ,jj  ,kk ,n);

                                       // Keep the interpolated velocity between the cell value and the wall value (0)
                                       if ( (velfp[n] > 0.0 && velfp[n] > vel_array(iloc,jloc,kloc,n)) ||
                                            (velfp[n] < 0.0 && velfp[n] < vel_array(iloc,jloc,kloc,n)) )
                                           velfp[n] = vel_array(iloc,jloc,kloc,n);
                                       if ( (gradp[n] > 0.0 && gradp[n] > gp_array(iloc,jloc,kloc,n)) ||
                                            (gradp[n] < 0.0 && gradp[n] < gp_array(iloc,jloc,kloc,n)) )
                                           gradp[n] = gp_array(iloc,jloc,kloc,n);
                                    }
                               } // Cut cell

                               particle.rdata(realData::dragx) = pbeta * ( velfp[0] - particle.rdata(realData::velx) ) -
                                                                (gradp[0] + gp0[0]) * particle.rdata(realData::volume);

                               particle.rdata(realData::dragy) = pbeta * ( velfp[1] - particle.rdata(realData::vely) ) -
                                                                (gradp[1] + gp0[1]) * particle.rdata(realData::volume);

                               particle.rdata(realData::dragz) = pbeta * ( velfp[2] - particle.rdata(realData::velz) ) -
                                                                (gradp[2] + gp0[2]) * particle.rdata(realData::volume);

                            } // Not covered
                        } // particle loop
                } // if box not all regular

            } // FAB not covered
        } // pti


        } // omp region


    } // lev

    // Reset velocity Dirichlet bc's to face values
    extrap_dir_bcs = 0;
    mfix_set_velocity_bcs(time, vel_g, extrap_dir_bcs);

}
