#include <mfix.H>
#include <mfix_util_F.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_EBMultiFabUtil.H>

// Compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
mfix::FillPatch (int lev, MultiFab& mf, MultiFab& cmf, MultiFab& fmf, int icomp, int ncomp)
{
#if 0
    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime);

        PhysBCFunct physbc(geom[lev],bcs,BndryFunctBase(phifill));
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                     geom[lev], physbc);
    }
    else
    {
        PhysBCFunct cphysbc(geom[lev-1],bcs,BndryFunctBase(phifill));
        PhysBCFunct fphysbc(geom[lev  ],bcs,BndryFunctBase(phifill));

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                   0, icomp, ncomp, geom[lev-1], geom[lev],
                                   cphysbc, fphysbc, refRatio(lev-1),
                                   mapper, bcs);
    }
#endif
}

//
// Set the BCs for all the variables EXCEPT pressure or velocity.
//
void
mfix::mfix_set_scalar_bcs ()
{
  BL_PROFILE("mfix::mfix_set_scalar_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
     for (MFIter mfi(*ep_g[lev], true); mfi.isValid(); ++mfi)
     {
        set_scalar_bcs ( BL_TO_FORTRAN_ANYD((*ep_g[lev])[mfi]),
                        (*ro_g[lev])[mfi].dataPtr (),
                        (*rop_g[lev])[mfi].dataPtr (),
                        (*mu_g[lev])[mfi].dataPtr (),
                        (*lambda_g[lev])[mfi].dataPtr (),
                        bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                        bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                        bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                        domain.loVect(), domain.hiVect(),
                        &nghost );
      }
        ep_g[lev] -> FillBoundary (geom[lev].periodicity());
        ro_g[lev] -> FillBoundary (geom[lev].periodicity());
       rop_g[lev] -> FillBoundary (geom[lev].periodicity());
        mu_g[lev] -> FillBoundary (geom[lev].periodicity());
    lambda_g[lev] -> FillBoundary (geom[lev].periodicity());
  }
}

//
// Set the BCs for velocity only
//
void
mfix::mfix_set_velocity_bcs (Real time, int extrap_dir_bcs)
{
  BL_PROFILE("mfix::mfix_set_velocity_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
     vel_g[lev] -> FillBoundary (geom[lev].periodicity());
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
     for (MFIter mfi(*vel_g[lev], true); mfi.isValid(); ++mfi)
     {
        set_velocity_bcs ( &time, 
                           BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
                           bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                           bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                           bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                           domain.loVect(), domain.hiVect(),
                           &nghost, &extrap_dir_bcs );
     }

     // Do this after as well as before to pick up terms that got updated in the call above
     vel_g[lev] -> FillBoundary (geom[lev].periodicity());
  }
}

//
// Fills ghost cell values of pressure appropriately for the BC type
//
void
mfix::mfix_extrap_pressure (int lev, std::unique_ptr<amrex::MultiFab>& p)
{
    BL_PROFILE("mfix::mfix_extrap_pressure()");
    if (nodal_pressure == 1) return;
 
    Box domain(geom[lev].Domain());
 
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for (MFIter mfi(*p, true); mfi.isValid(); ++mfi) {
 
        extrap_pressure_to_ghost_cells (
            BL_TO_FORTRAN_ANYD((*p)[mfi]),
            bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
            bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
            bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
            domain.loVect(), domain.hiVect(),
            &nghost);
    }
}

