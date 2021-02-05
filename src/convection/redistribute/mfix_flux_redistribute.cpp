#include <redistribution.H>

using namespace amrex;




void redistribution::flux_redistribute_eb (Box const& bx, int ncomp, int scomp,
                                           Array4<Real> const& dUdt,
                                           Array4<Real const> const& dUdt_in,
                                           Array4<Real> const& scratch,
                                           Array4<EBCellFlag const> const& flag,
                                           Array4<Real const> const& vfrac,
                                           Geometry& lev_geom)
{
    const Box dbox = lev_geom.growPeriodicDomain(2);

    Array4<Real> tmp(scratch, 0);
    Array4<Real> delm(scratch, ncomp);
    Array4<Real> wgt(scratch, 2*ncomp);

    Box const& bxg1 = amrex::grow(bx,1);
    Box const& bxg2 = amrex::grow(bx,2);

    // xxxxx TODO: more weight options
    amrex::ParallelFor(bxg2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        wgt(i,j,k) = (dbox.contains(IntVect(AMREX_D_DECL(i,j,k)))) ? 1.0 : 0.0;
    });

    amrex::ParallelFor(bxg1, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isSingleValued()) {
#if 1
            Real vtot = 0.0;
            Real divnc = 0.0;
            for (int kk = -1; kk <= 1; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
                 if ((ii != 0 or jj != 0 or kk != 0) and
                     flag(i,j,k).isConnected(ii,jj,kk) and
                    dbox.contains(IntVect(AMREX_D_DECL(i+ii,j+jj,k+kk))))
                {
                    Real vf = vfrac(i+ii,j+jj,k+kk);
                    vtot += vf;
                    divnc += vf * dUdt_in(i+ii,j+jj,k+kk,n+scomp);
                }
            }}}
            divnc /= (vtot + 1.e-80);
#else
            Real divnc = vfrac(i,j,k)*dUdt_in(i,j,k,n+scomp);
#endif
            Real optmp = (1.0-vfrac(i,j,k))*(divnc-dUdt_in(i,j,k,n+scomp));
            tmp(i,j,k,n+scomp) = optmp;
            delm(i,j,k,n+scomp) = -vfrac(i,j,k)*optmp;
        } else {
            tmp(i,j,k,n+scomp) = 0.0;
        }
    });

    amrex::ParallelFor(bxg1 & dbox, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isSingleValued()) {
            Real wtot = 0.0;
            for (int kk = -1; kk <= 1; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
                if ((ii != 0 or jj != 0 or kk != 0) and
                    flag(i,j,k).isConnected(ii,jj,kk))
                {
                    wtot += vfrac(i+ii,j+jj,k+kk) * wgt(i+ii,j+jj,k+kk);
                }
            }}}
            wtot = 1.0/(wtot+1.e-80);

            Real dtmp = delm(i,j,k,n+scomp) * wtot;
            for (int kk = -1; kk <= 1; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
                if ((ii != 0 or jj != 0 or kk != 0) and
                    bx.contains(IntVect(AMREX_D_DECL(i+ii,j+jj,k+kk))) and
                    flag(i,j,k).isConnected(ii,jj,kk))
                {
                    Gpu::Atomic::Add(&tmp(i+ii,j+jj,k+kk,n+scomp), dtmp*wgt(i+ii,j+jj,k+kk));
                }
            }}}
        }
    });

    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {

      Real org_div = dUdt(i,j,k,n+scomp);
      dUdt(i,j,k,n+scomp) = dUdt_in(i,j,k,n+scomp) + tmp(i,j,k,n+scomp);
    });

}







void redistribution::apply_eb_redistribution ( const Box& bx,
                                               amrex::Array4<amrex::Real>       const& div,
                                               amrex::Array4<amrex::Real const> const& divc,
                                               amrex::Array4<amrex::Real const> const& wt,
                                               MFIter const& mfi,
                                               const int icomp,
                                               const int ncomp,
                                               amrex::Array4<amrex::EBCellFlag const> const& flags,
                                               amrex::Array4<amrex::Real const> const& vfrac,
                                               const Geometry & geom)
{
  //
  // Check that grid is uniform
  //
  const Real* dx = geom.CellSize();

  if( ! amrex::almostEqual(dx[0],dx[1]) ||
      ! amrex::almostEqual(dx[0],dx[2]) ||
      ! amrex::almostEqual(dx[1],dx[2]) )
    amrex::Abort("apply_eb_redistribution(): grid spacing must be uniform");

  const Box dbox = geom.growPeriodicDomain(2);

   const Box& grown1_bx = amrex::grow(bx,1);
   const Box& grown2_bx = amrex::grow(bx,2);

   //
   // Working arrays
   //
   FArrayBox  delm_fab(grown1_bx,ncomp);
   FArrayBox  optmp_fab(grown2_bx,ncomp);
   FArrayBox  mask_fab(grown2_bx);

   Array4<Real> const& optmp = optmp_fab.array();
   Array4<Real> const& mask  = mask_fab.array();
   Array4<Real> const& delm  = delm_fab.array();

   //
   // Array "mask" is used to sever the link to ghost cells when the BCs
   // are not periodic
   // It is set to 1 when a cell can be used in computations, 0 otherwise
   //
   AMREX_FOR_3D(grown2_bx, i, j, k,
   {
       mask(i,j,k) = (dbox.contains(IntVect(AMREX_D_DECL(i,j,k)))) ? 1.0 : 0.0;
   });

   //
   // Init to zero tmp array
   //
   AMREX_FOR_4D(grown2_bx, ncomp, i, j, k, n,
   {
       optmp(i,j,k,n) = 0;
   });

   //
   // Step 2: compute delta M (mass gain or loss) on (lo-1,lo+1)
   //
   AMREX_FOR_4D(grown1_bx, ncomp, i, j, k, n,
   {
   if(flags(i,j,k).isSingleValued())
   {
       Real divnc(0.0);
       Real vtot(0.0);
       Real wted_frac(0.0);
       int  ks = (AMREX_SPACEDIM == 3) ? -1 : 0;
       int  ke = (AMREX_SPACEDIM == 3) ?  1 : 0;

       for (int kk(ks); kk <= ke; ++kk) {
           for (int jj(-1); jj <= 1; ++jj) {
               for (int ii(-1); ii <= 1; ++ii) {
                   if( (ii != 0 || jj != 0 || kk != 0) &&
                       flags(i,j,k).isConnected(ii,jj,kk) &&
                       dbox.contains(IntVect(AMREX_D_DECL(i+ii,j+jj,k+kk))))
                   {
                       wted_frac = vfrac(i+ii,j+jj,k+kk) * wt(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
                       vtot   += wted_frac;
                       divnc  += wted_frac * divc(i+ii,j+jj,k+kk,n);
                   }
               }
           }
       }
       divnc /=  (vtot + 1.e-80);

       // We need to multiply divc by mask to make sure optmp is zero for cells
       // outside the domain for non-cyclic BCs
       optmp(i,j,k,n) =  (1 - vfrac(i,j,k)) * (divnc - divc(i,j,k,n) * mask(i,j,k));
       delm(i,j,k,n)  = -(    vfrac(i,j,k)) * optmp(i,j,k,n);

   }
   else
   {
       delm(i,j,k,n) = 0;
   }
   });


   //
   // Step 3: redistribute excess/loss of mass
   //
   AMREX_FOR_4D(grown1_bx, ncomp, i, j, k, n,
   {
   if(flags(i,j,k).isSingleValued())
   {
       Real wtot(0.0);
       int  ks = (AMREX_SPACEDIM == 3) ? -1 : 0;
       int  ke = (AMREX_SPACEDIM == 3) ?  1 : 0;

       for (int kk(ks); kk <= ke; ++kk) {
         for (int jj(-1); jj <= 1; ++jj) {
           for (int ii(-1); ii <= 1; ++ii) {

                   if( (ii != 0 || jj != 0 || kk != 0) &&
                       (flags(i,j,k).isConnected(ii,jj,kk)) )
                   {
                       wtot += wt(i+ii,j+jj,k+kk) * vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
                   }

       }}}

       wtot = 1.0/(wtot + 1.e-80);

       for (int kk(ks); kk <= ke; ++kk) {
         for (int jj(-1); jj <= 1; ++jj) {
           for (int ii(-1); ii <= 1; ++ii) {

                   if( (ii != 0 || jj != 0 || kk != 0) &&
                       (flags(i,j,k).isConnected(ii,jj,kk)) &&
                       bx.contains(IntVect(AMREX_D_DECL(i+ii,j+jj,k+kk))) )
                   {
                       Gpu::Atomic::AddNoRet(&optmp(i+ii,j+jj,k+kk,n),
                                        delm(i,j,k,n) * wtot * mask(i+ii,j+jj,k+kk) * wt(i+ii,j+jj,k+kk));
                   }
           }}}

   }
   });

   //
   // Resume the correct sign, AKA return the negative
   //
   AMREX_FOR_4D(bx, ncomp, i, j, k, n,
   {
     div(i,j,k,icomp+n) = divc(i,j,k,n) + optmp(i,j,k,n);
   });

   Gpu::synchronize();
}
