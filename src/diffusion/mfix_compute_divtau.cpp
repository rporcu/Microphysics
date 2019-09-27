#include <mfix_redist_diff.hpp>
#include <mfix_diff_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLTensorOp.H>
#include <AMReX_MLEBTensorOp.H>

//
// Computation of divtau
//
void
mfix::mfix_compute_divtau ( Vector< std::unique_ptr<MultiFab> >& divtau,
                            Vector< std::unique_ptr<MultiFab> >& vel    )
{
   BL_PROFILE("mfix::mfix_diffuse_explicit");

   Vector<std::unique_ptr<MultiFab> > divtau_aux(nlev);
   for (int lev = 0; lev < nlev; lev++)
   {
     divtau_aux[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost,
           MFInfo(), *ebfactory[lev]));
     divtau_aux[lev]->setVal(0.0);
   }

   // Swap ghost cells and apply BCs to velocity
   Real time = 0;
   int extrap_dir_bcs = 0;
   mfix_set_velocity_bcs (time, vel, extrap_dir_bcs);

   // The boundary conditions need only be set once -- we do this at level 0
   int bc_lo[3], bc_hi[3];

   // Whole domain
   Box domain(geom[0].Domain());

   // Set BCs for Poisson's solver
   set_diff_bc (bc_lo, bc_hi,
                domain.loVect(), domain.hiVect(),
                &nghost,
                bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

   //
   // First define the operator "ebtensorop"
   //
   //       (alpha * a - beta * (del dot b grad)) sol
   //
   // LPInfo                       info;
   MLEBTensorOp ebtensorop(geom, grids, dmap, LPInfo().setMaxCoarseningLevel(0),
                           amrex::GetVecOfConstPtrs(ebfactory));

   // It is essential that we set MaxOrder of the solver to 2
   // if we want to use the standard sol(i)-sol(i-1) approximation
   // for the gradient at Dirichlet boundaries.
   // The solver's default order is 3 and this uses three points for the
   // gradient at a Dirichlet boundary.
   ebtensorop.setMaxOrder(2);

   // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
   ebtensorop.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                            {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );

   // Return div (mu grad)) phi 
   ebtensorop.setScalars(0.0, -1.0);

   // Compute the coefficients
   for (int lev = 0; lev < nlev; lev++)
   {
       average_cellcenter_to_face( GetArrOfPtrs(bcoeff_cc[lev]), *mu_g[lev], geom[lev] );

       bcoeff_cc[lev][0] -> FillBoundary(geom[lev].periodicity());
       bcoeff_cc[lev][1] -> FillBoundary(geom[lev].periodicity());
       bcoeff_cc[lev][2] -> FillBoundary(geom[lev].periodicity());

       ebtensorop.setShearViscosity  (lev, GetArrOfConstPtrs(bcoeff_cc[lev]));
       ebtensorop.setEBShearViscosity(lev, (*mu_g[lev]));

       ebtensorop.setLevelBC ( lev, GetVecOfConstPtrs(vel)[lev] );
   }

   MLMG solver(ebtensorop);

   solver.apply(GetVecOfPtrs(divtau_aux), GetVecOfPtrs(vel));

   for (int lev = 0; lev < nlev; lev++)
   {
      divtau_aux[lev]->FillBoundary(geom[lev].periodicity());
      MultiFab::Copy(*divtau[lev],*divtau_aux[lev],0,0,3,0);
   }

   const amrex::MultiFab* volfrac;

   const int cyclic_x = geom[0].isPeriodic(0) ? 1 : 0;
   const int cyclic_y = geom[0].isPeriodic(1) ? 1 : 0;
   const int cyclic_z = geom[0].isPeriodic(2) ? 1 : 0;

   for (int lev = 0; lev < nlev; lev++)
   {

      volfrac   = &(ebfactory[lev] -> getVolFrac());
      const Real* dx = geom[lev].CellSize();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*divtau[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) 
      {
         // Tilebox
         Box bx = mfi.tilebox ();

         // this is to check efficiently if this tile contains any eb stuff
         const EBFArrayBox&  vel_fab    = static_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
         const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();
   
         if (flags.getType(bx) != FabType::covered)
         {
            if (flags.getType(amrex::grow(bx,nghost)) != FabType::regular)
            {

                // Do redistribution at EB boundaries
                compute_redist_diff(
                   bx, *divtau[lev], *ep_g[lev], *divtau_aux[lev], &mfi, flags,
                   volfrac, cyclic_x, cyclic_y, cyclic_z, domain, dx);
            }
         }
      }

      // Divide by (ro_g ep_g)
      for (int n = 0; n < 3; n++)
      {
          MultiFab::Divide( *divtau[lev], *ep_g[lev], 0, n, 1, 0 );
          MultiFab::Divide( *divtau[lev], *ro_g[lev], 0, n, 1, 0 );
      }
   }
}
