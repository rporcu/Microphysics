#include <mfix_diff_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>

//
// Explicit diffusion
//
void
mfix::mfix_compute_divtau ( int lev,
                            MultiFab& divtau,
                            Vector< std::unique_ptr<MultiFab> >& vel,
                            int explicit_diffusion)
{
   BL_PROFILE("mfix::mfix_compute_divtau");
   Box domain(geom[lev].Domain());

   EB_set_covered(*vel[lev], 0, vel[lev]->nComp(), vel[lev]->nGrow(), covered_val);

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

      if (flags.getType(bx) == FabType::covered)
      {
         divtau[mfi].setVal(1.2345e200, bx, 0, 3);
      }
      else
      {
         if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular)
         {
            compute_divtau(
               BL_TO_FORTRAN_BOX(bx),
               BL_TO_FORTRAN_ANYD(divtau[mfi]),
               BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
               (*mu_g[lev])[mfi].dataPtr(),
               (*ro_g[lev])[mfi].dataPtr(),
               BL_TO_FORTRAN_ANYD((*ep_g[lev])[mfi]),
               domain.loVect (), domain.hiVect (),
               bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
               bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
               bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
               geom[lev].CellSize(), &nghost, &explicit_diffusion);
         }
         else
         {
            compute_divtau_eb(
               BL_TO_FORTRAN_BOX(bx),
               BL_TO_FORTRAN_ANYD(divtau[mfi]),
               BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
               (*mu_g[lev])[mfi].dataPtr(),
               (*ro_g[lev])[mfi].dataPtr(),
               BL_TO_FORTRAN_ANYD((*ep_g[lev])[mfi]),
               BL_TO_FORTRAN_ANYD(flags),
               BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
               BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
               BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
               BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
               BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
               BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
               BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
               BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
               domain.loVect (), domain.hiVect (),
               bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
               bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
               bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
               geom[lev].CellSize(), &nghost, &explicit_diffusion);

         }
      }
   }
   // Divide by (ro_g ep_g)
   for (int n = 0; n < 3; n++) 
   {
       MultiFab::Divide( divtau, *ep_g[lev], 0, n, 1, 0 );
       MultiFab::Divide( divtau, *ro_g[lev], 0, n, 1, 0 );
   }
}
