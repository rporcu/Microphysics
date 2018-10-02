#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>

void mfix::mfix_calc_drag_fluid()
{
    BL_PROFILE("mfix::mfix_calc_drag_fluid()");
    for (int lev = 0; lev < nlev; lev++)
    {
       Real dx = geom[lev].CellSize(0);
       Real dy = geom[lev].CellSize(1);
       Real dz = geom[lev].CellSize(2);

       bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                            (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

       if (OnSameGrids)
       {
       // ************************************************************
       // First create the beta of individual particles
       // ************************************************************
#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
       {
           const Box& sbx = (*ep_g[lev])[pti].box();
           auto& particles = pti.GetArrayOfStructs();
           const int np = particles.size();

           calc_particle_beta(
               sbx.loVect(), sbx.hiVect(),
               (*ep_g[lev])[pti].dataPtr() , (*ro_g[lev])[pti].dataPtr(),
               (*vel_g[lev])[pti].dataPtr(), (*mu_g[lev])[pti].dataPtr(),
               &np, particles.data(), &dx, &dy, &dz);
       }

       // ******************************************************************************
       // Now use the beta of individual particles to create the drag terms on the fluid
       // ******************************************************************************

       drag[lev]->setVal(0.0L);
       f_gds[lev]->setVal(0.0L);

       pc -> CalcDragOnFluid(*f_gds[lev], *drag[lev],
                             *bc_ilo[lev],*bc_ihi[lev],*bc_jlo[lev],*bc_jhi[lev],
                             *bc_klo[lev],*bc_khi[lev],nghost);
       }
       else
       {

       BoxArray            pba = pc->ParticleBoxArray(lev);
       DistributionMapping pdm = pc->ParticleDistributionMap(lev);

       // Temporary arrays
       int ng = ep_g[lev]->nGrow();
       std::unique_ptr<MultiFab> ep_g_pba(new MultiFab(pba,pdm,ep_g[lev]->nComp(),ng));
       ep_g_pba->copy(*ep_g[lev],0,0,1,ng,ng);
       ep_g_pba->FillBoundary(geom[lev].periodicity());

       ng = ro_g[lev]->nGrow();
       std::unique_ptr<MultiFab> ro_g_pba(new MultiFab(pba,pdm,ro_g[lev]->nComp(),ng));
       ro_g_pba->copy(*ro_g[lev],0,0,1,ng,ng);
       ro_g_pba->FillBoundary(geom[lev].periodicity());

       ng = mu_g[lev]->nGrow();
       std::unique_ptr<MultiFab> mu_g_pba(new MultiFab(pba,pdm,mu_g[lev]->nComp(),ng));
       mu_g_pba->copy(*mu_g[lev],0,0,1,ng,ng);
       mu_g_pba->FillBoundary(geom[lev].periodicity());

       ng = vel_g[lev]->nGrow();
       std::unique_ptr<MultiFab> vel_g_pba(new MultiFab(pba,pdm,vel_g[lev]->nComp(),ng));
       vel_g_pba->copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),ng,ng);
       vel_g_pba->FillBoundary(geom[lev].periodicity());

       // ************************************************************
       // First create the beta of individual particles
       // ************************************************************

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
       {
           const Box& sbx = (*ep_g_pba)[pti].box();
           auto& particles = pti.GetArrayOfStructs();
           const int np = particles.size();

           calc_particle_beta(
               sbx.loVect(), sbx.hiVect(),
               (*ep_g_pba)[pti].dataPtr(), (*ro_g_pba)[pti].dataPtr(),
               (*vel_g_pba)[pti].dataPtr(),  (*mu_g_pba)[pti].dataPtr(),
               &np, particles.data(), &dx, &dy, &dz);
       }

       // ******************************************************************************
       // Now use the beta of individual particles to create the drag terms on the fluid
       // ******************************************************************************

       std::unique_ptr<MultiFab>  drag_pba(new MultiFab(pba,pdm,drag[lev]->nComp(),drag[lev]->nGrow()));
       std::unique_ptr<MultiFab> f_gds_pba(new MultiFab(pba,pdm,f_gds[lev]->nComp(),f_gds[lev]->nGrow()));

       f_gds_pba->setVal(0.0L);
       drag_pba->setVal(0.0L);

       pc -> CalcDragOnFluid(*f_gds_pba,*drag_pba,
                             *bc_ilo[lev],*bc_ihi[lev],*bc_jlo[lev],*bc_jhi[lev],*bc_klo[lev],*bc_khi[lev],nghost);

       // Copy back from the dual grids.
        drag[lev] ->copy(*drag_pba);
       f_gds[lev] ->copy(*f_gds_pba);

       } // if not OnSameGrids

       // The projection method uses drag to update u, not (cell_vol * u), so we must divide by vol here
       //     and we will divide by density in the update.
       Real ovol = 1./(dx*dy*dz);
        drag[lev]->mult(ovol);
       f_gds[lev]->mult(ovol);
   
       // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
        drag[lev]->FillBoundary(geom[lev].periodicity());
       f_gds[lev]->FillBoundary(geom[lev].periodicity());

    } // lev
}

void
mfix::mfix_calc_drag_particle()
{
    BL_PROFILE("mfix::mfix_calc_drag_particle()");

    for (int lev = 0; lev < nlev; lev++)
    {

       Real dx = geom[lev].CellSize(0);
       Real dy = geom[lev].CellSize(1);
       Real dz = geom[lev].CellSize(2);

       bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                            (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

       Box domain(geom[lev].Domain());

       if (OnSameGrids)
       {
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
          mfix_set_velocity_bcs(extrap_dir_bcs);

          gp_tmp.FillBoundary(geom[lev].periodicity());

          int use_slopes = 0;
          if (use_slopes)
             mfix_compute_velocity_slopes( lev, vel_g );
#ifdef _OPENMP
#pragma omp parallel
#endif
          for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
          {
              auto& particles = pti.GetArrayOfStructs();
              const int np = particles.size();
   
              calc_drag_particle( BL_TO_FORTRAN_ANYD(       gp_tmp[pti]),
                                  BL_TO_FORTRAN_ANYD((  *gp0[lev])[pti]),
                                  BL_TO_FORTRAN_ANYD((*vel_g[lev])[pti]),
                                  (*xslopes[lev])[pti].dataPtr(),
                                  (*yslopes[lev])[pti].dataPtr(),
                                  BL_TO_FORTRAN_ANYD((*zslopes[lev])[pti]),
                                  &np, particles.data(), &dx, &dy, &dz, use_slopes);
          }

          // Reset velocity Dirichlet bc's to face values
          // HACK -- NOTE WE ARE CALLING THIS ON ALL LEVELS BUT ONLY NEED IT ON ONE LEVEL
          extrap_dir_bcs = 0;
          mfix_set_velocity_bcs(extrap_dir_bcs);
       }
#if 0
       else
       {

       BoxArray            pba = pc->ParticleBoxArray(lev);
       DistributionMapping pdm = pc->ParticleDistributionMap(lev);

       ng = gp[lev]->nGrow();
       std::unique_ptr<MultiFab> gp_tmp(new MultiFab(pba,pdm,gp[lev]->nComp(),ng));
       gp_tmp->copy(*gp[lev],0,0,gp[lev]->nComp(),ng,ng);
       gp_tmp->FillBoundary(geom[lev].periodicity());

       ng = gp0[lev]->nGrow();
       std::unique_ptr<MultiFab> gp0_tmp(new MultiFab(pba,pdm,gp0[lev]->nComp(),ng));
       gp0_tmp->copy(*gp0[lev],0,0,gp0[lev]->nComp(),ng,ng);
       gp0_tmp->FillBoundary(geom[lev].periodicity());

       ng = vel_g[lev]->nGrow();
       std::unique_ptr<MultiFab> vel_g_tmp(new MultiFab(tmp,pdm,vel_g[lev]->nComp(),ng));
       vel_g_tmp->copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),ng,ng);
       vel_g_tmp->FillBoundary(geom[lev].periodicity());

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
           set_gradp_bcs ( bx.loVect(), bx.hiVect(),
                           BL_TO_FORTRAN_ANYD(gp_tmp[mfi]),
                           bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                           bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                           bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                           domain.loVect(), domain.hiVect(),
                           &nghost);
       }

       int extrap_dir_bcs = 1;
       mfix_set_velocity_bcs(lev, extrap_dir_bcs);

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
       {
           auto& particles = pti.GetArrayOfStructs();
           const int np = particles.size();

           calc_drag_particle( BL_TO_FORTRAN_ANYD(     gp_tmp[pti]),
                               BL_TO_FORTRAN_ANYD((  *gp0_tmp[pti])),
                               BL_TO_FORTRAN_ANYD((*vel_g_tmp)[pti]),
                               &np, particles.data(), &dx, &dy, &dz);
       }
    }
#endif

    } // lev
}
