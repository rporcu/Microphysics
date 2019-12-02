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
#ifdef AMREX_USE_CUDA
#define BLOCK_SIZE 4
             dim3 threadsPerBlock(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
             dim3 numBlocks(std::ceil(bx.length()[0]/float(threadsPerBlock.x)),
                            std::ceil(bx.length()[1]/float(threadsPerBlock.y)),
                            std::ceil(bx.length()[2]/float(threadsPerBlock.z)));

             const unsigned shared_mem_size = sizeof(amrex::Real) * 3 * (
               (threadsPerBlock.x+2)*(threadsPerBlock.y)*(threadsPerBlock.z) +
               (threadsPerBlock.x)*(threadsPerBlock.y+2)*(threadsPerBlock.z) +
               (threadsPerBlock.x)*(threadsPerBlock.y)*(threadsPerBlock.z+2));

             IntVect bx_lo(bx.loVect());
             IntVect bx_hi(bx.hiVect());

             amrex::launch_global<<<numBlocks, threadsPerBlock, shared_mem_size,
               Gpu::gpuStream()>>>(
                 [state_fab,xs_fab,ys_fab,zs_fab,slopes_comp,bx_lo,bx_hi]
                 AMREX_GPU_DEVICE() noexcept
               {
                 int i = threadIdx.x;
                 int j = threadIdx.y;
                 int k = threadIdx.z;

                 int stride_i = blockIdx.x*blockDim.x + bx_lo[0];
                 int stride_j = blockIdx.y*blockDim.y + bx_lo[1];
                 int stride_k = blockIdx.z*blockDim.z + bx_lo[2];

                 int global_i = i + stride_i;
                 int global_j = j + stride_j;
                 int global_k = k + stride_k;

                 int local_limit_i = std::min(BLOCK_SIZE, bx_hi[0]-stride_i+1) + 1;
                 int local_limit_j = std::min(BLOCK_SIZE, bx_hi[1]-stride_j+1) + 1;
                 int local_limit_k = std::min(BLOCK_SIZE, bx_hi[2]-stride_k+1) + 1;

                 int global_limit_i = std::min(stride_i+BLOCK_SIZE, bx_hi[0]+1);
                 int global_limit_j = std::min(stride_j+BLOCK_SIZE, bx_hi[1]+1);
                 int global_limit_k = std::min(stride_k+BLOCK_SIZE, bx_hi[2]+1);

                 __shared__ Real state_fab_i[BLOCK_SIZE+2][BLOCK_SIZE][BLOCK_SIZE][3];
                 __shared__ Real state_fab_j[BLOCK_SIZE][BLOCK_SIZE+2][BLOCK_SIZE][3];
                 __shared__ Real state_fab_k[BLOCK_SIZE][BLOCK_SIZE][BLOCK_SIZE+2][3];

                 if((global_i <= bx_hi[0]) and (global_j <= bx_hi[1]) and (global_k <= bx_hi[2]))
                 {
                   i++; j++; k++;

                   int local_list_i[3] = {0, i, local_limit_i};
                   int local_list_j[3] = {0, j, local_limit_j};
                   int local_list_k[3] = {0, k, local_limit_k};

                   int global_list_i[3] = {stride_i-1, global_i, global_limit_i};
                   int global_list_j[3] = {stride_j-1, global_j, global_limit_j};
                   int global_list_k[3] = {stride_k-1, global_k, global_limit_k};

                   for(unsigned int ii(0); ii < 3; ii++) {
                     int t_i = local_list_i[ii];
                     int g_i = global_list_i[ii];

                     for(unsigned int n(0); n < 3; n++)
                       state_fab_i[t_i][j][k][n] = state_fab(g_i,global_j,global_k,n);
                   }
                   
                   for(unsigned int jj(0); jj < 3; jj++) {
                     int t_j = local_list_j[jj];
                     int g_j = global_list_j[jj];

                     for(unsigned int n(0); n < 3; n++)
                       state_fab_j[i][t_j][k][n] = state_fab(global_i,g_j,global_k,n);
                   }
                   
                   for(unsigned int kk(0); kk < 3; kk++) {
                     int t_k = local_list_k[kk];
                     int g_k = global_list_k[kk];

                     for(unsigned int n(0); n < 3; n++)
                       state_fab_k[i][j][t_k][n] = state_fab(global_i,global_j,g_k,n);
                   }

                   __syncthreads();

                   for(unsigned int n(0); n < 3; n++) {
                     // X direction
                     Real du_xl = 2.0*(state_fab_i[i  ][j][k][n]- state_fab_i[i-1][j][k][n]);
                     Real du_xr = 2.0*(state_fab_i[i+1][j][k][n]- state_fab_i[i  ][j][k][n]);
                     Real du_xc = 0.5*(state_fab_i[i+1][j][k][n]- state_fab_i[i-1][j][k][n]);

                     Real xslope = amrex::min(std::abs(du_xl), std::abs(du_xc), std::abs(du_xr));
                     xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                     xs_fab(global_i,global_j,global_k,slopes_comp+n) = (du_xc > 0.0) ? xslope : -xslope;

                     // Y direction
                     Real du_yl = 2.0*(state_fab_j[i][j  ][k][n] - state_fab_j[i][j-1][k][n]);
                     Real du_yr = 2.0*(state_fab_j[i][j+1][k][n] - state_fab_j[i][j  ][k][n]);
                     Real du_yc = 0.5*(state_fab_j[i][j+1][k][n] - state_fab_j[i][j-1][k][n]);

                     Real yslope = amrex::min(std::abs(du_yl), std::abs(du_yc), std::abs(du_yr));
                     yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                     ys_fab(global_i,global_j,global_k,slopes_comp+n) = (du_yc > 0.0) ? yslope : -yslope;

                     // Z direction
                     Real du_zl = 2.0*(state_fab_k[i][j][k  ][n] - state_fab_k[i][j][k-1][n]);
                     Real du_zr = 2.0*(state_fab_k[i][j][k+1][n] - state_fab_k[i][j][k  ][n]);
                     Real du_zc = 0.5*(state_fab_k[i][j][k+1][n] - state_fab_k[i][j][k-1][n]);

                     Real zslope = amrex::min(std::abs(du_zl), std::abs(du_zc), std::abs(du_zr));
                     zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                     zs_fab(global_i,global_j,global_k,slopes_comp+n) = (du_zc > 0.0) ? zslope : -zslope;
                   }
                 }
               });
#else
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
#endif
           }
           else
           {
               const auto& flag_fab = flags.array();

             amrex::ParallelFor(bx, ncomp,
               [state_fab,xs_fab,ys_fab,zs_fab,slopes_comp,flag_fab]
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
                       // X direction
                       Real du_xl = (flag_fab(i-1,j,k).isCovered()) ? 0.0 :
                           2.0*(state_fab(i,j,k,n) - state_fab(i-1,j,k,n));
                       Real du_xr = (flag_fab(i+1,j,k).isCovered()) ? 0.0 :
                           2.0*(state_fab(i+1,j,k,n) - state_fab(i,j,k,n));
                       Real du_xc = 0.5*(state_fab(i+1,j,k,n) - state_fab(i-1,j,k,n));

                       Real xslope = amrex::min(std::abs(du_xl), std::abs(du_xc), std::abs(du_xr));
                       xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                       xs_fab(i,j,k,slopes_comp+n) = (du_xc > 0.0) ? xslope : -xslope;

                       // Y direction
                       Real du_yl = (flag_fab(i,j-1,k).isCovered()) ? 0.0 :
                           2.0*(state_fab(i,j,k,n) - state_fab(i,j-1,k,n));
                       Real du_yr = (flag_fab(i,j+1,k).isCovered()) ? 0.0 :
                           2.0*(state_fab(i,j+1,k,n) - state_fab(i,j,k,n));
                       Real du_yc = 0.5*(state_fab(i,j+1,k,n) - state_fab(i,j-1,k,n));

                       Real yslope = amrex::min(std::abs(du_yl), std::abs(du_yc), std::abs(du_yr));
                       yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                       ys_fab(i,j,k,slopes_comp+n) = (du_yc > 0.0) ? yslope : -yslope;

                       // Z direction
                       Real du_zl = (flag_fab(i,j,k-1).isCovered()) ? 0.0 :
                           2.0*(state_fab(i,j,k,n) - state_fab(i,j,k-1,n));
                       Real du_zr = (flag_fab(i,j,k+1).isCovered()) ? 0.0 :
                           2.0*(state_fab(i,j,k+1,n) - state_fab(i,j,k,n));
                       Real du_zc = 0.5*(state_fab(i,j,k+1,n) - state_fab(i,j,k-1,n));

                       Real zslope = amrex::min(std::abs(du_zl), std::abs(du_zc), std::abs(du_zr));
                       zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                       zs_fab(i,j,k,slopes_comp+n) = (du_zc > 0.0) ? zslope : -zslope;
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
