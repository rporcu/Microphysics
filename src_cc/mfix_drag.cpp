#include <AMReX_ParmParse.H>
#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix.H>
#include <mfix_des_F.H>
#include <mfix_util_F.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>

void mfix::mfix_calc_drag_fluid(Real time)
{
    BL_PROFILE("mfix::mfix_calc_drag_fluid()");

    for (int lev = 0; lev < nlev; lev++)
    {
        Real dx = geom[lev].CellSize(0);
        Real dy = geom[lev].CellSize(1);
        Real dz = geom[lev].CellSize(2);

        bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                             (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

        MultiFab* ep_ptr;
        MultiFab* ro_ptr;
        MultiFab* mu_ptr;
        MultiFab* vel_ptr;
        MultiFab* drag_ptr;
        MultiFab* f_gds_ptr;

        // This will be temporaries only for the dual grid case
        std::unique_ptr<MultiFab>  ep_g_pba;
        std::unique_ptr<MultiFab>  ro_g_pba;
        std::unique_ptr<MultiFab>  mu_g_pba;
        std::unique_ptr<MultiFab> vel_g_pba;
        std::unique_ptr<MultiFab>  drag_pba;
        std::unique_ptr<MultiFab> f_gds_pba;

        if (OnSameGrids)
        {
            ep_ptr    =  ep_g[lev].get();
            ro_ptr    =  ro_g[lev].get();
            mu_ptr    =  mu_g[lev].get();
            vel_ptr   = vel_g[lev].get();
            drag_ptr  =  drag[lev].get();
            f_gds_ptr = f_gds[lev].get();
        }
        else
        {
            BoxArray            pba = pc->ParticleBoxArray(lev);
            DistributionMapping pdm = pc->ParticleDistributionMap(lev);

            // Temporary arrays  -- CHECK COPY() ROUTINE
            int ng = ep_g[lev]->nGrow();
            ep_g_pba.reset(new MultiFab(pba,pdm,ep_g[lev]->nComp(),ng));
            ep_g_pba->copy(*ep_g[lev],0,0,1,ng,ng);
            ep_g_pba->FillBoundary(geom[lev].periodicity());

            ng = ro_g[lev]->nGrow();
            ro_g_pba.reset(new MultiFab(pba,pdm,ro_g[lev]->nComp(),ng));
            ro_g_pba->copy(*ro_g[lev],0,0,1,ng,ng);
            ro_g_pba->FillBoundary(geom[lev].periodicity());

            ng = mu_g[lev]->nGrow();
            mu_g_pba.reset(new MultiFab(pba,pdm,mu_g[lev]->nComp(),ng));
            mu_g_pba->copy(*mu_g[lev],0,0,1,ng,ng);
            mu_g_pba->FillBoundary(geom[lev].periodicity());


            EBFArrayBoxFactory ebfactory_loc( * eb_level_fluid, geom[lev], pba, pdm,
                                              {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                               m_eb_full_grow_cells}, EBSupport::basic);

            ng = vel_g[lev]->nGrow();
            vel_g_pba.reset(new MultiFab(pba,pdm,vel_g[lev]->nComp(),ng, MFInfo(), ebfactory_loc));
            vel_g_pba->copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),ng,ng);
            vel_g_pba->FillBoundary(geom[lev].periodicity());

            drag_pba.reset(new MultiFab(pba,pdm,drag[lev]->nComp(),drag[lev]->nGrow()));

            f_gds_pba.reset(new MultiFab(pba,pdm,f_gds[lev]->nComp(),f_gds[lev]->nGrow()));

            ep_ptr    =  ep_g_pba.get();
            ro_ptr    =  ro_g_pba.get();
            mu_ptr    =  mu_g_pba.get();
            vel_ptr   = vel_g_pba.get();
            drag_ptr  =  drag_pba.get();
            f_gds_ptr = f_gds_pba.get();

        }

        // Create fab to host reconstructed velocity field
        FArrayBox vel_r;

        // Phi is always on the particles grid
        const MultiFab & phi = * ls[lev];



#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
        {

            auto& particles = pti.GetArrayOfStructs();
            const int np = particles.size();

            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_ptr)[pti]);
            const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();
            Box bx = pti.tilebox ();

            if (flags.getType(amrex::grow(bx,0)) != FabType::covered)
            {
                if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
                {

                    calc_particle_beta( BL_TO_FORTRAN_ANYD((*ep_ptr)[pti]),
                                        (*ro_ptr)[pti].dataPtr(),
                                        (*mu_ptr)[pti].dataPtr(),
                                        BL_TO_FORTRAN_ANYD((*vel_ptr)[pti]),
                                        &np, particles.data(),
                                        geom[lev].ProbLo(), geom[lev].CellSize());
                }
                else
                {

                    Box gbox = amrex::grow(bx,2);
                    vel_r.resize(gbox,3);

                    reconstruct_velocity( BL_TO_FORTRAN_ANYD(vel_r),
                                          BL_TO_FORTRAN_ANYD((*vel_ptr)[pti]),
                                          BL_TO_FORTRAN_ANYD((phi)[pti]),
                                          1,
                                          BL_TO_FORTRAN_ANYD(flags),
                                          geom[lev].ProbLo(), geom[lev].CellSize());

                    calc_particle_beta( BL_TO_FORTRAN_ANYD((*ep_ptr)[pti]),
                                        (*ro_ptr)[pti].dataPtr(),
                                        (*mu_ptr)[pti].dataPtr(),
                                        BL_TO_FORTRAN_ANYD(vel_r),
                                        &np, particles.data(),
                                        geom[lev].ProbLo(), geom[lev].CellSize());
                }
            }

        }

        // ******************************************************************************
        // Now use the beta of individual particles to create the drag terms on the fluid
        // ******************************************************************************

        drag_ptr  -> setVal(0.0L);
        f_gds_ptr -> setVal(0.0L);

        pc -> CalcDragOnFluid(*f_gds_ptr, *drag_ptr, *particle_ebfactory[lev],
                              *bc_ilo[lev],*bc_ihi[lev],*bc_jlo[lev],*bc_jhi[lev],
                              *bc_klo[lev],*bc_khi[lev],nghost);


        if (!OnSameGrids)
        {
            // Copy back from the dual grids.
            drag[lev] ->copy(*drag_pba);
            f_gds[lev]->copy(*f_gds_pba);
        }

        // The projection method uses drag to update u, not (cell_vol * u), so we must divide by vol here
        //     and we will divide by density in the update.
        Real ovol = 1./(dx*dy*dz);
        drag[lev]->mult(ovol);
        f_gds[lev]->mult(ovol);

        // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
        drag[lev] -> FillBoundary(geom[lev].periodicity());
        f_gds[lev]-> FillBoundary(geom[lev].periodicity());

    } // lev
}

void
mfix::mfix_calc_drag_particle(Real time)
{
    BL_PROFILE("mfix::mfix_calc_drag_particle()");

    for (int lev = 0; lev < nlev; lev++)
    {

        Box      domain(geom[lev].Domain());
        MultiFab gp_tmp, gp0_tmp;

        gp_tmp.define(grids[lev],dmap[lev],3,1);

        MultiFab::Copy(gp_tmp, *gp[lev], 0, 0, 3, 1);
        gp_tmp.FillBoundary(geom[lev].periodicity());

        //
        // NOTE -- it is essential that we call set_gradp_bcs after calling FillBoundary
        //         because the set_gradp_bcs call hopefully sets the ghost cells exterior
        //         to the domain from ghost cells interior to the domain
        //
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(gp_tmp, true); mfi.isValid(); ++mfi)
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
        MultiFab* gp0_ptr;
        MultiFab* vel_ptr;

        // Temporararies for dual grid case
        std::unique_ptr<MultiFab>   gp_pba;
        std::unique_ptr<MultiFab>  gp0_pba;
        std::unique_ptr<MultiFab>  vel_pba;

        if (OnSameGrids)
        {
            gp_ptr  = &gp_tmp;
            gp0_ptr = gp0[lev].get();
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

            ng = gp0[lev]->nGrow();
            gp0_pba.reset(new MultiFab(pba,pdm,gp0[lev]->nComp(),ng));
            gp0_pba->copy(*gp0[lev],0,0,gp0[lev]->nComp(),ng,ng);
            gp0_pba->FillBoundary(geom[lev].periodicity());

            EBFArrayBoxFactory ebfactory_loc( * eb_level_fluid,
                                              geom[lev], pba, pdm,
                                              {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                               m_eb_full_grow_cells}, EBSupport::basic);

            ng = vel_g[lev]->nGrow();
            vel_pba.reset(new MultiFab(pba,pdm,vel_g[lev]->nComp(),ng,MFInfo(), ebfactory_loc));
            vel_pba->copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),ng,ng);
            vel_pba->FillBoundary(geom[lev].periodicity());

            gp_ptr  = gp_pba.get();
            gp0_ptr = gp0_pba.get();
            vel_ptr = vel_pba.get();

        }


        // Create fab to host reconstructed velocity field
        FArrayBox vel_r;

        // Phi is always on the particles grid
        const MultiFab & phi = * ls[lev];


#ifdef _OPENMP
#pragma omp parallel
#endif
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
                        calc_drag_particle( BL_TO_FORTRAN_ANYD((*gp_ptr)[pti]),
                                            BL_TO_FORTRAN_ANYD((*gp0_ptr)[pti]),
                                            BL_TO_FORTRAN_ANYD((*vel_ptr)[pti]),
                                            &np, particles.data(),
                                            geom[lev].CellSize(), geom[lev].ProbLo());
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
                                              geom[lev].ProbLo(), geom[lev].CellSize()) ;

                        calc_drag_particle_eb( BL_TO_FORTRAN_ANYD((*gp_ptr)[pti]),
                                               BL_TO_FORTRAN_ANYD((*gp0_ptr)[pti]),
                                               BL_TO_FORTRAN_ANYD(vel_r),
                                               BL_TO_FORTRAN_ANYD(flags),
                                               &np, particles.data(),
                                               geom[lev].CellSize(), geom[lev].ProbLo());

                    }
                }

            }

            // Reset velocity Dirichlet bc's to face values
            // HACK -- NOTE WE ARE CALLING THIS ON ALL LEVELS BUT ONLY NEED IT ON ONE LEVEL
            extrap_dir_bcs = 0;
            mfix_set_velocity_bcs(time, extrap_dir_bcs);

    } // lev
}
