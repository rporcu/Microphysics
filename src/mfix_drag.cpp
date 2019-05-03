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
        Box domain(geom[lev].Domain());
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

        // Phi is always on the particles grid
        const MultiFab & phi = * level_sets[lev];

        int band_width = 2;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            // Create fab to host reconstructed velocity field
            FArrayBox vel_r;

            Real x0 = geom[lev].ProbLo(0);
            Real y0 = geom[lev].ProbLo(1);
            Real z0 = geom[lev].ProbLo(2);

            Real velfp[3];
            Real gradp[3];

            const auto dxi = geom[lev].InvCellSizeArray();
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
                    if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
                    {
                        const auto& vel_array = vel_ptr->array(pti);
                        const auto&  gp_array =  gp_ptr->array(pti);

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
                    else
                    {
                        Box gbox = amrex::grow(bx,2);
                        vel_r.resize(gbox,3);

                        reconstruct_velocity( BL_TO_FORTRAN_ANYD(vel_r),
                                              BL_TO_FORTRAN_ANYD((*vel_ptr)[pti]),
                                              BL_TO_FORTRAN_ANYD((phi)[pti]),
                                              1, //level_set -> get_ls_ref(),
                                              BL_TO_FORTRAN_ANYD(flags),
                                              geom[lev].ProbLo(), geom[lev].CellSize(),
                                              &band_width);

                        const auto& vel_array     = vel_r.array();
                        const auto&  gp_array     =  gp_ptr->array(pti);
                        const auto& flags_array = flags.array();

                        for (int ip = 0; ip < np; ++ip)
                        {
                            MFIXParticleContainer::ParticleType& particle = particles[ip];
                            Real pbeta = particle.rdata(realData::dragx);

                            // This identifies which cell the particle is in
                            int iloc = std::floor((particle.pos(0) - x0)*dxi[0]);
                            int jloc = std::floor((particle.pos(1) - y0)*dxi[1]);
                            int kloc = std::floor((particle.pos(2) - z0)*dxi[2]);

                            // This part was in trilinear_interp before
                            // Pick upper cell in the stencil
                            Real lx = (particle.pos(0) - x0)*dxi[0] + 0.5;
                            Real ly = (particle.pos(1) - y0)*dxi[1] + 0.5;
                            Real lz = (particle.pos(2) - z0)*dxi[2] + 0.5;

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

                            } else {

                               trilinear_interp(particle, &velfp[0], vel_array, plo, dxi);

                               //
                               // gradp is interpolated only if there are not covered cells
                               // in the stencil. If any covered cell is present in the stencil,
                               // we use the cell gradp
                               //

                               if (flags_array(i-1,j,k).isCovered() ||
                                   flags_array(i  ,j,k).isCovered() ||
                                   flags_array(i,j-1,k).isCovered() ||
                                   flags_array(i,j  ,k).isCovered() ||
                                   flags_array(i,j,k-1).isCovered() ||
                                   flags_array(i,j,k  ).isCovered() )
                               {
                                   for (int n = 0; n < 3; n++)
                                      gradp[n] = gp_array(iloc,jloc,kloc,n);

                               } else {

                                   trilinear_interp(particle, &gradp[0], gp_array, plo, dxi);

                               }

                               particle.rdata(realData::dragx) = pbeta * ( velfp[0] - particle.rdata(realData::velx) ) -
                                                                (gradp[0] + gp0[0]) * particle.rdata(realData::volume);

                               particle.rdata(realData::dragy) = pbeta * ( velfp[1] - particle.rdata(realData::vely) ) -
                                                                (gradp[1] + gp0[1]) * particle.rdata(realData::volume);

                               particle.rdata(realData::dragz) = pbeta * ( velfp[2] - particle.rdata(realData::velz) ) -
                                                                (gradp[2] + gp0[2]) * particle.rdata(realData::volume);
                            } // cell not covered

                        } // particle loop
                } // if box not all regular

            } // if not covered
        } // pti

        // Reset velocity Dirichlet bc's to face values
        // HACK -- NOTE WE ARE CALLING THIS ON ALL LEVELS BUT ONLY NEED IT ON ONE LEVEL
        extrap_dir_bcs = 0;
        mfix_set_velocity_bcs(time, extrap_dir_bcs);

        } // omp region

    } // lev
}
