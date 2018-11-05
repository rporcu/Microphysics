#include <mfix_F.H>
#include <mfix_proj_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLNodeLaplacian.H>

//
// Compute div(ep_g * u)
// 
void
mfix::mfix_compute_diveu (Real time)
{
    if (nodal_pressure == 1)
    {
        Vector<std::unique_ptr<MultiFab> > epu;
        epu.resize(nlev);

        int extrap_dir_bcs = 0;
        mfix_set_velocity_bcs (time, extrap_dir_bcs);

        for (int lev = 0; lev < nlev; lev++)
        {
            // We only need one ghost cell here -- so no need to make it bigger
            epu[lev].reset(new MultiFab(vel_g[lev]->boxArray(), vel_g[lev]->DistributionMap(),
                                        vel_g[lev]->nComp()   , 1 , MFInfo(),
                                        *ebfactory[lev]));

            epu[lev]->setVal(1.e200);

            Box domain(geom[lev].Domain());

            MultiFab::Copy(*epu[lev], *vel_g[lev], 0, 0, 3, epu[lev]->nGrow() );
   
            for (int n = 0; n < 3; n++)
                MultiFab::Multiply( *epu[lev], *ep_g[lev], 0, n, 1, epu[lev]->nGrow() );

            epu[lev]->FillBoundary (geom[lev].periodicity());     

#ifdef _OPENMP
#pragma omp parallel
#endif
           // Extrapolate Dirichlet values to ghost cells -- but do it differently in that 
           //  no-slip walls are treated exactly like slip walls -- this is only relevant
           //  when going into the projection
           for (MFIter mfi((*epu[lev]), true); mfi.isValid(); ++mfi)
           {
                set_vec_bcs ( BL_TO_FORTRAN_ANYD((*epu[lev])[mfi]),
                              bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                              bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                              bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                              domain.loVect(), domain.hiVect(),
                              &nghost);
           }

           epu[lev]->FillBoundary (geom[lev].periodicity());
       }

        // Define the operator in order to compute the multi-level divergence
        //
        //        (del dot b sigma grad)) phi
        //
        LPInfo                       info;
        MLNodeLaplacian              matrix(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));

        // Set domain BCs for Poisson's solver
        // The domain BCs refer to level 0 only
        int bc_lo[3], bc_hi[3];
        Box domain(geom[0].Domain());

        set_ppe_bc(bc_lo, bc_hi,
                   domain.loVect(), domain.hiVect(),
                   &nghost,
                   bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                   bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                   bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

        matrix.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                             {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );

        matrix.compDivergence(GetVecOfPtrs(diveu), GetVecOfPtrs(epu)); 

    }
    else
    {

       int extrap_dir_bcs = 1;
       mfix_set_velocity_bcs (time, extrap_dir_bcs);

       for (int lev = 0; lev < nlev; lev++)
       {

          Box domain(geom[lev].Domain());
          vel_g[lev]->FillBoundary (geom[lev].periodicity());     

          // Create face centered multifabs for ep_g and vel_g
          MultiFab epu( vel_g[lev]->boxArray(),  vel_g[lev]-> DistributionMap(),
                        vel_g[lev]->nComp(), nghost, MFInfo(), *ebfactory[lev]);

          MultiFab::Copy( epu, *vel_g[lev], 0, 0, 3, vel_g[lev]->nGrow() );
   
          for (int n = 0; n < 3; n++)
             MultiFab::Multiply( epu, *ep_g[lev], 0, n, 1, vel_g[lev]->nGrow() );

          Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> epu_fc;
          mfix_average_cc_to_fc( lev, epu, epu_fc );
   
          // This does not need to have correct ghost values in place
          EB_computeDivergence( *diveu[lev], GetArrOfConstPtrs(epu_fc), geom[lev] );
       }
    }

    // Restore velocities to carry Dirichlet values on faces
    int extrap_dir_bcs = 0;
    mfix_set_velocity_bcs (time, extrap_dir_bcs);
}
