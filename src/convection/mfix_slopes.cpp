#include <mfix.H>

//
// Compute the slopes of Sborder (velocity, density or tracer)
//
void
mfix::mfix_compute_slopes (int lev, Real time, MultiFab& Sborder,
                           Vector<std::unique_ptr<MultiFab>>& xslopes_in,
                           Vector<std::unique_ptr<MultiFab>>& yslopes_in,
                           Vector<std::unique_ptr<MultiFab>>& zslopes_in,
                           int slopes_comp)
{
    BL_PROFILE("mfix::mfix_compute_slopes");

    EB_set_covered(Sborder, 0, Sborder.nComp(), 1, covered_val);

    Box domain(geom[lev].Domain());

    int ncomp = Sborder.nComp();

    // We initialize slopes to zero in the grown domain ... this is essential
    //    to handle the up-winding at outflow faces
    xslopes_in[lev]->setVal(0.0, slopes_comp, ncomp, xslopes_in[lev]->nGrow());
    yslopes_in[lev]->setVal(0.0, slopes_comp, ncomp, yslopes_in[lev]->nGrow());
    zslopes_in[lev]->setVal(0.0, slopes_comp, ncomp, zslopes_in[lev]->nGrow());

    // ... then set them to this large number in the interior in order to be sure
    //     that no "bad" values go unnoticed
    xslopes_in[lev]->setVal(1.2345e300, slopes_comp, ncomp, 0);
    yslopes_in[lev]->setVal(1.2345e300, slopes_comp, ncomp, 0);
    zslopes_in[lev]->setVal(1.2345e300, slopes_comp, ncomp, 0);

    const auto cellcent = &(ebfactory[lev] -> getCentroid());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Sborder,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       // Tilebox
       Box bx = mfi.tilebox();

       // This is to check efficiently if this tile contains any eb stuff
       const EBFArrayBox& Sborder_fab = static_cast<EBFArrayBox const&>(Sborder[mfi]);
       const EBCellFlagFab& flags = Sborder_fab.getEBCellFlagFab();

       if (flags.getType(amrex::grow(bx,0)) != FabType::covered)
       {
           const auto& state_fab = Sborder.array(mfi);
           const auto& xs_fab = xslopes_in[lev]->array(mfi);
           const auto& ys_fab = yslopes_in[lev]->array(mfi);
           const auto& zs_fab = zslopes_in[lev]->array(mfi);

           // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
           if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
           {
             amrex::ParallelFor(bx, ncomp,
               [state_fab,xs_fab,ys_fab,zs_fab,slopes_comp]
               AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
               {
                   // X direction
                   Real du_xl = 2.0*(state_fab(i  ,j,k,n) - state_fab(i-1,j,k,n));
                   Real du_xr = 2.0*(state_fab(i+1,j,k,n) - state_fab(i  ,j,k,n));
                   Real du_xc = 0.5*(state_fab(i+1,j,k,n) - state_fab(i-1,j,k,n));

                   Real xslope = amrex::min(std::abs(du_xl), std::abs(du_xc), std::abs(du_xr));
                   xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                   xs_fab(i,j,k,slopes_comp+n) = (du_xc > 0.0) ? xslope : -xslope;

                   // Y direction
                   Real du_yl = 2.0*(state_fab(i,j  ,k,n) - state_fab(i,j-1,k,n));
                   Real du_yr = 2.0*(state_fab(i,j+1,k,n) - state_fab(i,j  ,k,n));
                   Real du_yc = 0.5*(state_fab(i,j+1,k,n) - state_fab(i,j-1,k,n));

                   Real yslope = amrex::min(std::abs(du_yl), std::abs(du_yc), std::abs(du_yr));
                   yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                   ys_fab(i,j,k,slopes_comp+n) = (du_yc > 0.0) ? yslope : -yslope;

                   // Z direction
                   Real du_zl = 2.0*(state_fab(i,j,k  ,n) - state_fab(i,j,k-1,n));
                   Real du_zr = 2.0*(state_fab(i,j,k+1,n) - state_fab(i,j,k  ,n));
                   Real du_zc = 0.5*(state_fab(i,j,k+1,n) - state_fab(i,j,k-1,n));

                   Real zslope = amrex::min(std::abs(du_zl), std::abs(du_zc), std::abs(du_zr));
                   zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                   zs_fab(i,j,k,slopes_comp+n) = (du_zc > 0.0) ? zslope : -zslope;
               });
           }
           else
           {

               const auto& flag_fab = flags.array();
               const auto& ccent_fab = cellcent->array(mfi);

               amrex::ParallelFor(bx, ncomp,
               [state_fab,xs_fab,ys_fab,zs_fab,slopes_comp,flag_fab,ccent_fab]
               AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
               {
                   if (flag_fab(i,j,k).isCovered())
                   {
                       xs_fab(i,j,k,slopes_comp+n) = 0.0;
                       ys_fab(i,j,k,slopes_comp+n) = 0.0;
                       zs_fab(i,j,k,slopes_comp+n) = 0.0;
                   }
                   else
                   {

                     amrex::Real A[8][2];
                     amrex::Real du[8];

                     {
                       int lc=0;
                       for(int jj(-1); jj<=1; jj++){
                         for(int ii(-1); ii<=1; ii++){
                           if ( ii == 0 and jj==0 ) continue;

                           if( flag_fab(i,j,k).isConnected(ii,jj,0)){

                             // Not multplying by dx to be consistent with how the
                             // slope is stored. Also not including the global shift
                             // wrt plo or i,j,k. We only need relative distance.

                             A[lc][0] = ii + ccent_fab(i+ii,j+jj,k,0) - ccent_fab(i,j,k,0);
                             A[lc][1] = jj + ccent_fab(i+ii,j+jj,k,1) - ccent_fab(i,j,k,1);

                             du[lc] = state_fab(i+ii,j+jj,k,n) - state_fab(i,j,k,n);

                           } else {

                             A[lc][0] = 0.0;
                             A[lc][1] = 0.0;

                             du[lc] = 0.0;
                           }

                           lc++;
                         }
                       }
                     }

                     amrex::Real a11(0.0), a12(0.0), a22(0.0);
                     amrex::Real  b1(0.0),  b2(0.0);

                     for(int lc(0); lc<8; ++lc){
                       a11 += A[lc][0]* A[lc][0];
                       a12 += A[lc][0]* A[lc][1];
                       a22 += A[lc][1]* A[lc][1];
                       b1  += A[lc][0]*du[lc];
                       b2  += A[lc][1]*du[lc];
                     }

                     amrex::Real det = 1.0/(a11*a22 - a12*a12);

                     // Slope at centroid of (i,j,k)
                     amrex::Real dudx = (a22*b1 - a12*b2) * det;
                     amrex::Real dudy = (a11*b2 - a12*b1) * det;


                     // X direction
                     if(flag_fab(i  ,j,k).isSingleValued() or
                        (flag_fab(i-1,j,k).isSingleValued() or not flag_fab(i,j,k).isConnected(-1,0,0)) or
                        (flag_fab(i+1,j,k).isSingleValued() or not flag_fab(i,j,k).isConnected( 1,0,0))) {

                       xs_fab(i,j,k,slopes_comp+n) = dudx;

                     } else {

                       Real du_xl = 2.0*(state_fab(i  ,j,k,n) - state_fab(i-1,j,k,n));
                       Real du_xr = 2.0*(state_fab(i+1,j,k,n) - state_fab(i  ,j,k,n));
                       Real du_xc = 0.5*(state_fab(i+1,j,k,n) - state_fab(i-1,j,k,n));

                       Real xslope = amrex::min(std::abs(du_xl), std::abs(du_xc), std::abs(du_xr));
                       xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                       xs_fab(i,j,k,slopes_comp+n) = (du_xc > 0.0) ? xslope : -xslope;

                     }


                     // Y direction
                     if(flag_fab(i,j  ,k).isSingleValued() or
                       (flag_fab(i,j-1,k).isSingleValued() or not flag_fab(i,j,k).isConnected(0,-1,0)) or
                       (flag_fab(i,j+1,k).isSingleValued() or not flag_fab(i,j,k).isConnected(0, 1,0))) {

                       ys_fab(i,j,k,slopes_comp+n) = dudy;

                     } else {

                       Real du_yl = 2.0*(state_fab(i,j  ,k,n) - state_fab(i,j-1,k,n));
                       Real du_yr = 2.0*(state_fab(i,j+1,k,n) - state_fab(i,j  ,k,n));
                       Real du_yc = 0.5*(state_fab(i,j+1,k,n) - state_fab(i,j-1,k,n));

                       Real yslope = amrex::min(std::abs(du_yl), std::abs(du_yc), std::abs(du_yr));
                       yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                       ys_fab(i,j,k,slopes_comp+n) = (du_yc > 0.0) ? yslope : -yslope;
                     }


                     // Z direction
                     if(flag_fab(i,j,k  ).isSingleValued() or
                       (flag_fab(i,j,k-1).isSingleValued() or not flag_fab(i,j,k).isConnected(0,0,-1)) or
                       (flag_fab(i,j,k+1).isSingleValued() or not flag_fab(i,j,k).isConnected(0,0, 1))) {

                       zs_fab(i,j,k,slopes_comp+n) = 0.0;

                     } else {

                       Real du_zl = 2.0*(state_fab(i,j,k  ,n) - state_fab(i,j,k-1,n));
                       Real du_zr = 2.0*(state_fab(i,j,k+1,n) - state_fab(i,j,k  ,n));
                       Real du_zc = 0.5*(state_fab(i,j,k+1,n) - state_fab(i,j,k-1,n));

                       Real zslope = amrex::min(std::abs(du_zl), std::abs(du_zc), std::abs(du_zr));
                       zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                       zs_fab(i,j,k,slopes_comp+n) = (du_zc > 0.0) ? zslope : -zslope;

                     }
                   }
               });
           } // end of cut cell region

           const int minf = bc_list.get_minf();

           const auto& flag_fab = flags.array();

           const auto& ilo_ifab = bc_ilo[lev]->array();
           const auto& ihi_ifab = bc_ihi[lev]->array();
           const auto& jlo_ifab = bc_jlo[lev]->array();
           const auto& jhi_ifab = bc_jhi[lev]->array();
           const auto& klo_ifab = bc_klo[lev]->array();
           const auto& khi_ifab = bc_khi[lev]->array();

           amrex::ParallelFor(bx, ncomp,
             [domain,flag_fab,ilo_ifab,ihi_ifab,jlo_ifab,jhi_ifab,klo_ifab,khi_ifab,state_fab,
              xs_fab,ys_fab,zs_fab,minf,slopes_comp]
             AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
             {
               if ( (i == domain.smallEnd(0)) && !flag_fab(i,j,k).isCovered() && ilo_ifab(i-1,j,k,0) == minf)
               {
                   Real du_xl = 2.0*(state_fab(i  ,j,k,n) - state_fab(i-1,j,k,n));
                   Real du_xr = 2.0*(state_fab(i+1,j,k,n) - state_fab(i  ,j,k,n));
                   Real du_xc = (state_fab(i+1,j,k,n)+3.0*state_fab(i,j,k,n)-4.0*state_fab(i-1,j,k,n))/3.0;

                   Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                   xslope          = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                   xs_fab(i,j,k,slopes_comp+n) = (du_xc       > 0.0) ? xslope : -xslope;
               }
               if ( (i == domain.bigEnd(0)) && !flag_fab(i,j,k).isCovered() && ihi_ifab(i+1,j,k,0) == minf)
               {
                   Real du_xl = 2.0*(state_fab(i  ,j,k,n) - state_fab(i-1,j,k,n));
                   Real du_xr = 2.0*(state_fab(i+1,j,k,n) - state_fab(i  ,j,k,n));
                   Real du_xc = -(state_fab(i-1,j,k,n)+3.0*state_fab(i,j,k,n)-4.0*state_fab(i+1,j,k,n))/3.0;

                   Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                   xslope          = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                   xs_fab(i,j,k,slopes_comp+n) = (du_xc       > 0.0) ? xslope : -xslope;
               }

               if ( (j == domain.smallEnd(1)) && !flag_fab(i,j,k).isCovered() && jlo_ifab(i,j-1,k,0) == minf)
               {
                   Real du_yl = 2.0*(state_fab(i,j  ,k,n) - state_fab(i,j-1,k,n));
                   Real du_yr = 2.0*(state_fab(i,j+1,k,n) - state_fab(i,j  ,k,n));
                   Real du_yc = (state_fab(i,j+1,k,n)+3.0*state_fab(i,j,k,n)-4.0*state_fab(i,j-1,k,n))/3.0;

                   Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                   yslope          = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                   ys_fab(i,j,k,slopes_comp+n) = (du_yc       > 0.0) ? yslope : -yslope;
               }
               if ( (j == domain.bigEnd(1)) && !flag_fab(i,j,k).isCovered() && jhi_ifab(i,j+1,k,0) == minf)
               {
                   Real du_yl = 2.0*(state_fab(i,j  ,k,n) - state_fab(i,j-1,k,n));
                   Real du_yr = 2.0*(state_fab(i,j+1,k,n) - state_fab(i,j  ,k,n));
                   Real du_yc = -(state_fab(i,j-1,k,n)+3.0*state_fab(i,j,k,n)-4.0*state_fab(i,j+1,k,n))/3.0;

                   Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                   yslope          = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                   ys_fab(i,j,k,slopes_comp+n) = (du_yc       > 0.0) ? yslope : -yslope;
               }

               if ( (k == domain.smallEnd(2)) && !flag_fab(i,j,k).isCovered() && klo_ifab(i,j,k-1,0) == minf)
               {
                   Real du_zl = 2.0*(state_fab(i,j,k  ,n) - state_fab(i,j,k-1,n));
                   Real du_zr = 2.0*(state_fab(i,j,k+1,n) - state_fab(i,j,k  ,n));
                   Real du_zc = (state_fab(i,j,k+1,n)+3.0*state_fab(i,j,k,n)-4.0*state_fab(i,j,k-1,n))/3.0;

                   Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                   zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                   zs_fab(i,j,k,slopes_comp+n) = (du_zc       > 0.0) ? zslope : -zslope;
               }
               if ( (k == domain.bigEnd(2)) && !flag_fab(i,j,k).isCovered() && khi_ifab(i,j,k+1,0) == minf)
               {
                   Real du_zl = 2.0*(state_fab(i,j,k  ,n) - state_fab(i,j,k-1,n));
                   Real du_zr = 2.0*(state_fab(i,j,k+1,n) - state_fab(i,j,k  ,n));
                   Real du_zc = -(state_fab(i,j,k-1,n)+3.0*state_fab(i,j,k,n)-4.0*state_fab(i,j,k+1,n))/3.0;

                   Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                   zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                   zs_fab(i,j,k,slopes_comp+n) = (du_zc       > 0.0) ? zslope : -zslope;
               }
           });
           
        } // not covered
    } // MFIter

    xslopes_in[lev]->FillBoundary(geom[lev].periodicity());
    yslopes_in[lev]->FillBoundary(geom[lev].periodicity());
    zslopes_in[lev]->FillBoundary(geom[lev].periodicity());
}
