#include <mfix.H>
#include <mfix_F.H>
#include <mfix_util_F.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_EBMultiFabUtil.H>

namespace
{
  mfix* mfix_for_fillpatching;
}

void set_ptr_to_mfix(mfix& mfix_for_fillpatching_in)
{
   mfix_for_fillpatching = &mfix_for_fillpatching_in;
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
// We can't get around this so instead we create an mfix object
//    and use that to access the quantities that aren't passed here.
inline
void VelFillBox (Box const& bx, FArrayBox& dest,
                 const int dcomp, const int numcomp,
                 GeometryData const& geom, const Real time_in,
                 const BCRec* bcr, const int bcomp,
                 const int orig_comp)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in VelFillBox");
    if (numcomp != 3)
         amrex::Abort("Must have numcomp = 3 in VelFillBox");

    const Box& domain = geom.Domain();

    AmrParGDB* my_gdb = mfix_for_fillpatching->GetParGDB();

    // This is a bit hack-y but does get us the right level 
    int lev = 0;
    for (int ilev = 0; ilev < 10; ilev++)
    {
       const Geometry& lev_geom = mfix_for_fillpatching->GetParGDB()->Geom(ilev);
       if (domain.length()[0] == (lev_geom.Domain()).length()[0]) 
       {
         lev = ilev;
         break;
       }
    }

    // We are hard-wiring this fillpatch routine to define the Dirichlet values
    //    at the faces (not the ghost cell center)
    int extrap_dir_bcs = 0;

    // We only do this to make it not const
    Real time = time_in;

    const int* bc_ilo_ptr = mfix_for_fillpatching->get_bc_ilo_ptr(lev);
    const int* bc_ihi_ptr = mfix_for_fillpatching->get_bc_ihi_ptr(lev);
    const int* bc_jlo_ptr = mfix_for_fillpatching->get_bc_jlo_ptr(lev);
    const int* bc_jhi_ptr = mfix_for_fillpatching->get_bc_jhi_ptr(lev);
    const int* bc_klo_ptr = mfix_for_fillpatching->get_bc_klo_ptr(lev);
    const int* bc_khi_ptr = mfix_for_fillpatching->get_bc_khi_ptr(lev);

    int nghost = mfix_for_fillpatching->get_nghost();

#if 0

    const Real* dx = geom.CellSize();
    Real xlo[AMREX_SPACEDIM];

    for (int icomp = 0; icomp < numcomp; ++icomp)
    {
        generic_fill(BL_TO_FORTRAN_N_ANYD(dest,dcomp+icomp), 
                     BL_TO_FORTRAN_BOX(domain),
                     &icomp, dx, xlo, &time, bcr[bcomp+icomp].vect());
    }

#else
    set_velocity_bcs ( &time, 
                       BL_TO_FORTRAN_ANYD(dest), 
                       bc_ilo_ptr, bc_ihi_ptr, 
                       bc_jlo_ptr, bc_jhi_ptr, 
                       bc_klo_ptr, bc_khi_ptr, 
                       domain.loVect(), domain.hiVect(),
                       &nghost, &extrap_dir_bcs );
#endif
}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
mfix::FillPatchVel (int lev, Real time, MultiFab& mf, int icomp, int ncomp, const Vector<BCRec>& bcs)
{
    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetDataVel(0, time, smf, stime);

        CpuBndryFuncFab bfunc(VelFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                    geom[lev], physbc, 0);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataVel(lev-1, time, cmf, ctime);
        GetDataVel(lev  , time, fmf, ftime);

        CpuBndryFuncFab bfunc(VelFillBox);
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, 0);

    }

    // if (lev == 1) std::cout << " MF IN FILLPATCH " <<  mf[0] << std::endl;
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void
mfix::GetDataVel (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (lev > 0)
    {
       std::cout << "GETTING DATA AT LEVEL " << lev << std::endl;
       std::cout << "GETTING DATA AT TIME " << t_old[lev] << " < " << time << " < " << t_new[lev] << std::endl;

       std::cout << "NORM OF OLD U " <<  mfix_norm0(vel_go,lev,0)  << std::endl;;
       std::cout << "NORM OF OLD V " <<  mfix_norm0(vel_go,lev,1)  << std::endl;;
       std::cout << "NORM OF OLD W " <<  mfix_norm0(vel_go,lev,2)  << std::endl;;

       std::cout << "NORM OF NEW U " <<  mfix_norm0(vel_g,lev,0)  << std::endl;;
       std::cout << "NORM OF NEW V " <<  mfix_norm0(vel_g,lev,1)  << std::endl;;
       std::cout << "NORM OF NEW W " <<  mfix_norm0(vel_g,lev,2)  << std::endl;;
    }

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        data.push_back(vel_g[lev].get());
        datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        data.push_back(vel_go[lev].get());
        datatime.push_back(t_old[lev]);
    }
    else
    {
        data.push_back(vel_go[lev].get());
        data.push_back(vel_g[lev].get());
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
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

