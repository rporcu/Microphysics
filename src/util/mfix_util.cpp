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
    bool ropg_has_nans = rop_g[lev] -> contains_nan (0);

    if (ug_has_nans)
  amrex::Print() << "WARNING: u_g contains NaNs!!!";

    if (vg_has_nans)
  amrex::Print() << "WARNING: v_g contains NaNs!!!";

    if (wg_has_nans)
  amrex::Print() << "WARNING: w_g contains NaNs!!!";

    if (pg_has_nans)
  amrex::Print() << "WARNING: p_g contains NaNs!!!";

    if (ropg_has_nans)
  amrex::Print() << "WARNING: rop_g contains NaNs!!!";

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
#pragma omp parallel
#endif
       for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi)
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
// Subroutine to compute norm0 of EB multifab
//
Real
mfix::mfix_norm0 ( const Vector< std::unique_ptr<MultiFab>>& mf, int lev, int comp )
{
   MultiFab mf_tmp( mf[lev]->boxArray(), mf[lev]->DistributionMap(), mf[lev]->nComp(),
                    0,  MFInfo(), *ebfactory[lev]);

   MultiFab::Copy( mf_tmp, *mf[lev], comp, comp, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );

   return mf_tmp.norm0( comp );
}

Real
mfix::mfix_norm0 ( MultiFab& mf, int lev, int comp )
{
   MultiFab mf_tmp( mf.boxArray(), mf.DistributionMap(), mf.nComp(),
                    0,  MFInfo(), *ebfactory[lev]);

   MultiFab::Copy( mf_tmp, mf, comp, comp, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );

   return mf_tmp.norm0( comp );
}

//
// Subroutine to compute norm1 of EB multifab
//
Real
mfix::mfix_norm1 ( const Vector< std::unique_ptr<MultiFab>>& mf, int lev, int comp )
{
   MultiFab mf_tmp( mf[lev]->boxArray(), mf[lev]->DistributionMap(), mf[lev]->nComp(),
                    0,  MFInfo(), *ebfactory[lev]);

   MultiFab::Copy( mf_tmp, *mf[lev], comp, comp, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );

   return mf_tmp.norm1( comp, geom[lev].periodicity() );
}

Real
mfix::mfix_norm1 ( MultiFab& mf, int lev, int comp )
{
   MultiFab mf_tmp( mf.boxArray(), mf.DistributionMap(), mf.nComp(),
                    0,  MFInfo(), *ebfactory[lev]);

   MultiFab::Copy( mf_tmp, mf, comp, comp, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );

   return mf_tmp.norm1( comp, geom[lev].periodicity() );
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

#pragma gpu
  mfix_sum_mf(AMREX_INT_ANYD(lo),AMREX_INT_ANYD(hi),BL_TO_FORTRAN_N_ANYD(fab,comp),
          AMREX_REAL_ANYD(dx),BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                    AMREX_MFITER_REDUCE_SUM(&sum));
    }

    if (!local)
  ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}
