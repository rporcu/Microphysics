#include <mfix.H>
#include <mfix_util_F.H>
#include <AMReX_EBMultiFabUtil.H>

void
mfix::check_for_nans (int lev)
{
    bool ug_has_nans = vel_g[lev] -> contains_nan (0);
    bool vg_has_nans = vel_g[lev] -> contains_nan (1);
    bool wg_has_nans = vel_g[lev] -> contains_nan (2);
    bool pg_has_nans =   p_g[lev] -> contains_nan (0);

    if (ug_has_nans)
  amrex::Print() << "WARNING: u_g contains NaNs!!!";

    if (vg_has_nans)
  amrex::Print() << "WARNING: v_g contains NaNs!!!";

    if (wg_has_nans)
  amrex::Print() << "WARNING: w_g contains NaNs!!!";

    if (pg_has_nans)
  amrex::Print() << "WARNING: p_g contains NaNs!!!";
}

//
// Print the maximum values of the velocity components
//
void
mfix::mfix_print_max_vel(int lev)
{
    amrex::Print() << "max(abs(u/v/w/p))  = " <<
       mfix_norm0(vel_g, lev, 0) << "  " <<
       mfix_norm0(vel_g, lev, 1) << "  " <<
       mfix_norm0(vel_g, lev, 2) << "  " <<
       mfix_norm0(p_g,   lev, 0) << "  " << std::endl;
}


//
// Print the maximum values of the pressure gradient components
//
void
mfix::mfix_print_max_gp (int lev)
{
    amrex::Print() << "max(abs(gpx/gpy/gpz))  = " <<
       mfix_norm0(gp, lev, 0) << "  " <<
       mfix_norm0(gp, lev, 1) << "  " <<
       mfix_norm0(gp, lev, 2) << "  " << std::endl;
}


void
mfix::mfix_compute_vort ()
{
    BL_PROFILE("mfix::mfix_compute_vort");

    for (int lev = 0; lev < nlev; lev++)
    {
       Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(*vel_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
          // Tilebox
          Box bx = mfi.tilebox ();

          // This is to check efficiently if this tile contains any eb stuff
          const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_g[lev])[mfi]);
          const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

          if (flags.getType(amrex::grow(bx,0)) == FabType::regular )
          {
            compute_vort (
                        BL_TO_FORTRAN_BOX(bx),
                        BL_TO_FORTRAN_ANYD((* vort[lev])[mfi]),
                        BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
                        geom[lev].CellSize());
          } else {
             vort[lev]->setVal( 0.0, bx, 0, 1);
          }
       }
    }
}

//
// Compute norm0 of EB multifab
//
Real
mfix::mfix_norm0 ( const Vector< std::unique_ptr<MultiFab>>& mf, int lev, int comp )
{
   int ncomp = 1;
   int ngrow = 0;
   MultiFab mf_tmp( mf[lev]->boxArray(), mf[lev]->DistributionMap(), ncomp, ngrow, 
                    MFInfo(), *ebfactory[lev]);

   MultiFab::Copy( mf_tmp, *mf[lev], comp, 0, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );

   return mf_tmp.norm0( 0 );
}

//
// Compute norm0 of EB multifab
//
Real
mfix::mfix_norm0 ( MultiFab& mf, int lev, int comp )
{
   int ncomp = 1;
   int ngrow = 0;
   MultiFab mf_tmp( mf.boxArray(), mf.DistributionMap(), ncomp, ngrow, 
                    MFInfo(), *ebfactory[lev]);

   MultiFab::Copy( mf_tmp, mf, comp, 0, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );

   return mf_tmp.norm0( 0 );
}

//
// Compute max of EB multifab
//
Real
mfix::mfix_max ( MultiFab& mf, int lev, int comp )
{
   int ncomp = 1;
   int ngrow = 0;
   MultiFab mf_tmp( mf.boxArray(), mf.DistributionMap(), ncomp, ngrow, 
                    MFInfo(), *ebfactory[lev]);

   MultiFab::Copy( mf_tmp, mf, comp, 0, 1, 0 );
   EB_set_covered( mf_tmp, -1.e100 );

   return mf_tmp.max( 0 );
}

//
// Compute min of EB multifab
//
Real
mfix::mfix_min ( MultiFab& mf, int lev, int comp )
{
   int ncomp = 1;
   int ngrow = 0;
   MultiFab mf_tmp( mf.boxArray(), mf.DistributionMap(), ncomp, ngrow, 
                    MFInfo(), *ebfactory[lev]);

   MultiFab::Copy( mf_tmp, mf, comp, 0, 1, 0 );
   EB_set_covered( mf_tmp, 1.e100 );

   return mf_tmp.min( 0 );
}

//
// Compute norm1 of EB multifab
//
Real
mfix::mfix_norm1 ( const Vector< std::unique_ptr<MultiFab>>& mf, int lev, int comp )
{
   int ncomp = 1;
   int ngrow = 0;
   MultiFab mf_tmp( mf[lev]->boxArray(), mf[lev]->DistributionMap(), ncomp, ngrow, 
                    MFInfo(), *ebfactory[lev]);

   MultiFab::Copy( mf_tmp, *mf[lev], comp, 0, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );

   return mf_tmp.norm1( 0, geom[lev].periodicity() );
}

//
// Compute norm1 of EB multifab
//
Real
mfix::mfix_norm1 ( MultiFab& mf, int lev, int comp )
{
   int ncomp = 1;
   int ngrow = 0;
   MultiFab mf_tmp( mf.boxArray(), mf.DistributionMap(), ncomp, ngrow, 
                    MFInfo(), *ebfactory[lev]);

   MultiFab::Copy( mf_tmp, mf, comp, 0, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );

   return mf_tmp.norm1( 0, geom[lev].periodicity() );
}

//
// Subroutine to compute norm0 of EB multifab
//
Real
mfix::mfix_norm0 ( const Vector< std::unique_ptr<MultiFab>>& mf1,
                   const Vector< std::unique_ptr<MultiFab>>& mf2,
                   int lev, int comp1, int comp2 )
{
   BL_ASSERT((*mf1[lev]).boxArray() == (*mf2[lev]).boxArray());

   int ncomp = 1;
   int ngrow = 0;
   MultiFab mf_tmp( mf1[lev]->boxArray(), mf1[lev]->DistributionMap(), ncomp, ngrow,
                    MFInfo(), *ebfactory[lev]);

   MultiFab::Copy    ( mf_tmp, *mf1[lev], comp1, 0, 1, 0 );
   MultiFab::Multiply( mf_tmp, *mf2[lev], comp2, 0, 1, 0 );

   EB_set_covered( mf_tmp, 0.0 );

   return mf_tmp.norm0( 0 );
}


//
// Compute norm0 (max norm / inf norm) of product of EB multifabs
//
Real
mfix::mfix_norm0 ( MultiFab& mf1, MultiFab& mf2, int lev, int comp1, int comp2)
{
   BL_ASSERT(mf1.boxArray() == mf2.boxArray());

   int ncomp = 1;
   int ngrow = 0;
   MultiFab mf_tmp( mf1.boxArray(), mf1.DistributionMap(), ncomp, ngrow,  
                    MFInfo(), *ebfactory[lev]);

   MultiFab::Copy    ( mf_tmp, mf1, comp1, 0, 1, 0 );
   MultiFab::Multiply( mf_tmp, mf2, comp2, 0, 1, 0 );

   EB_set_covered( mf_tmp, 0.0 );

   return mf_tmp.norm0( 0 );
}

Real
mfix::volWgtSum (int lev, const MultiFab& mf, int comp, bool local)
{
    BL_PROFILE("mfix::volWgtSum()");

    Real        sum     = 0.0;
    const Real* dx      = geom[lev].CellSize();

    const MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif
    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
        const FArrayBox& fab = mf[mfi];

        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        mfix_sum_mf(lo, hi, BL_TO_FORTRAN_N_ANYD(fab, comp),
                    dx, BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                    &sum);
    }

    if (!local)
        ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}
