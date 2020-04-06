#include <mfix.H>

//
// Compute the slopes of Sborder (velocity, density or tracer)
//
void
mfix::mfix_compute_slopes (int lev,
                           Real time,
                           MultiFab& Sborder,
                           Vector< MultiFab* > const& xslopes_in,
                           Vector< MultiFab* > const& yslopes_in,
                           Vector< MultiFab* > const& zslopes_in,
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

           // No cut cells in tile + 1-cell width halo -> use non-eb routine
           if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
           {
             amrex::ParallelFor(bx, ncomp,
               [state_fab,xs_fab,ys_fab,zs_fab,slopes_comp]
               AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
               {
                   const Real state = state_fab(i,j,k,n);
                   const Real state_x_mns = state_fab(i-1,j,k,n);
                   const Real state_x_pls = state_fab(i+1,j,k,n);
                   const Real state_y_mns = state_fab(i,j-1,k,n);
                   const Real state_y_pls = state_fab(i,j+1,k,n);
                   const Real state_z_mns = state_fab(i,j,k-1,n);
                   const Real state_z_pls = state_fab(i,j,k+1,n);

                   // X direction
                   Real du_xl = 2.0*(state - state_x_mns);
                   Real du_xr = 2.0*(state_x_pls - state);
                   Real du_xc = 0.5*(state_x_pls - state_x_mns);

                   Real xslope = amrex::min(std::abs(du_xl), std::abs(du_xc), std::abs(du_xr));
                   xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                   xs_fab(i,j,k,slopes_comp+n) = (du_xc > 0.0) ? xslope : -xslope;

                   // Y direction
                   Real du_yl = 2.0*(state - state_y_mns);
                   Real du_yr = 2.0*(state_y_pls - state);
                   Real du_yc = 0.5*(state_y_pls - state_y_mns);

                   Real yslope = amrex::min(std::abs(du_yl), std::abs(du_yc), std::abs(du_yr));
                   yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                   ys_fab(i,j,k,slopes_comp+n) = (du_yc > 0.0) ? yslope : -yslope;

                   // Z direction
                   Real du_zl = 2.0*(state - state_z_mns);
                   Real du_zr = 2.0*(state_z_pls - state);
                   Real du_zc = 0.5*(state_z_pls - state_z_mns);

                   Real zslope = amrex::min(std::abs(du_zl), std::abs(du_zc), std::abs(du_zr));
                   zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                   zs_fab(i,j,k,slopes_comp+n) = (du_zc > 0.0) ? zslope : -zslope;
               });
           }
           else
           {

               const auto& flag_fab = flags.array();
               const auto& ccent_fab = cellcent->array(mfi);

        bool kernel_run_with_launch_global(false);

#ifdef AMREX_USE_CUDA
#define BLOCK_SIZE 4
        dim3 threadsPerBlock(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
        dim3 numBlocks(std::ceil(bx.length()[0]/float(threadsPerBlock.x)),
                       std::ceil(bx.length()[1]/float(threadsPerBlock.y)),
                       std::ceil(bx.length()[2]/float(threadsPerBlock.z)));

        const unsigned sm_size =
          6*sizeof(amrex::Real)*(BLOCK_SIZE+2)*(BLOCK_SIZE+2)*(BLOCK_SIZE+2) +
          sizeof(amrex::EBCellFlag)*(BLOCK_SIZE+2)*(BLOCK_SIZE+2)*(BLOCK_SIZE+2);

        if(sm_size <= Gpu::Device::sharedMemPerBlock()) {
          amrex::IntVect bx_lo(bx.loVect());
          amrex::IntVect bx_hi(bx.hiVect());

          amrex::launch_global<<<numBlocks,threadsPerBlock,sm_size,Gpu::gpuStream()>>>(
            [state_fab,xs_fab,ys_fab,zs_fab,slopes_comp,flag_fab,ccent_fab,ncomp,
             bx_lo,bx_hi] AMREX_GPU_DEVICE () noexcept
          {
            int tid_x = threadIdx.x;
            int tid_y = threadIdx.y;
            int tid_z = threadIdx.z;

            const int stride_x = blockIdx.x*BLOCK_SIZE + bx_lo[0];
            const int stride_y = blockIdx.y*BLOCK_SIZE + bx_lo[1];
            const int stride_z = blockIdx.z*BLOCK_SIZE + bx_lo[2];

            int i = tid_x + stride_x;
            int j = tid_y + stride_y;
            int k = tid_z + stride_z;

            __shared__ Real ccent_sm[BLOCK_SIZE+2][BLOCK_SIZE+2][BLOCK_SIZE+2][3];
            __shared__ Real state_sm[BLOCK_SIZE+2][BLOCK_SIZE+2][BLOCK_SIZE+2][3];
            __shared__ EBCellFlag flags_sm[BLOCK_SIZE+2][BLOCK_SIZE+2][BLOCK_SIZE+2];

            if((i <= bx_hi[0]) and (j <= bx_hi[1]) and (k <= bx_hi[2]))
            {
              tid_x += 1; tid_y += 1; tid_z += 1;

              int local_list_x[3] = {0, tid_x, amrex::min(BLOCK_SIZE, bx_hi[0]-stride_x+1)+1};
              int local_list_y[3] = {0, tid_y, amrex::min(BLOCK_SIZE, bx_hi[1]-stride_y+1)+1};
              int local_list_z[3] = {0, tid_z, amrex::min(BLOCK_SIZE, bx_hi[2]-stride_z+1)+1};

              int global_list_x[3] = {stride_x-1, i, amrex::min(stride_x+BLOCK_SIZE, bx_hi[0]+1)};
              int global_list_y[3] = {stride_y-1, j, amrex::min(stride_y+BLOCK_SIZE, bx_hi[1]+1)};
              int global_list_z[3] = {stride_z-1, k, amrex::min(stride_z+BLOCK_SIZE, bx_hi[2]+1)};

              for (int ii(0); ii < 3; ii++)
              for (int jj(0); jj < 3; jj++)
              for (int kk(0); kk < 3; kk++) {
                int t_x(local_list_x[ii]);
                int t_y(local_list_y[jj]);
                int t_z(local_list_z[kk]);

                int g_x(global_list_x[ii]);
                int g_y(global_list_y[jj]);
                int g_z(global_list_z[kk]);

                for (int n(0); n < ncomp; n++) {
                  ccent_sm[t_x][t_y][t_z][n] = ccent_fab(g_x,g_y,g_z,n);
                  state_sm[t_x][t_y][t_z][n] = state_fab(g_x,g_y,g_z,n);
                }

                flags_sm[t_x][t_y][t_z] = flag_fab(g_x,g_y,g_z);
              }

              __syncthreads();

              for (int n(0); n < ncomp; n++)
              {
                if (flags_sm[tid_x][tid_y][tid_z].isCovered())
                {
                  xs_fab(i,j,k,slopes_comp+n) = 0.0;
                  ys_fab(i,j,k,slopes_comp+n) = 0.0;
                  zs_fab(i,j,k,slopes_comp+n) = 0.0;
                }
                else
                {
                  amrex::Real A[27][3];
                  amrex::Real du[27];

                  {
                    int lc = 0;

                    for(int kk(-1); kk<=1; kk++)
                    for(int jj(-1); jj<=1; jj++)
                    for(int ii(-1); ii<=1; ii++) {
                      if((ii!=0 or jj!=0 or kk!=0) and flags_sm[tid_x][tid_y][tid_z].isConnected(ii,jj,kk)) {
                        // Not multiplying by dx to be consistent with how the
                        // slope is stored. Also not including the global shift
                        // wrt plo or i,j,k. We only need relative distance.
                        A[lc][0] = ii+ccent_sm[tid_x+ii][tid_y+jj][tid_z+kk][0]-ccent_sm[tid_x][tid_y][tid_z][0];
                        A[lc][1] = jj+ccent_sm[tid_x+ii][tid_y+jj][tid_z+kk][1]-ccent_sm[tid_x][tid_y][tid_z][1];
                        A[lc][2] = kk+ccent_sm[tid_x+ii][tid_y+jj][tid_z+kk][2]-ccent_sm[tid_x][tid_y][tid_z][2];

                        du[lc] = state_sm[tid_x+ii][tid_y+jj][tid_z+kk][n] - state_sm[tid_x][tid_y][tid_z][n];
                      }
                      else {
                        A[lc][0] = 0.0;
                        A[lc][1] = 0.0;
                        A[lc][2] = 0.0;

                        du[lc] = 0.0;
                      }

                      lc++;
                    }
                  }

                  amrex::Real AtA[3][3];
                  amrex::Real Atb[3];

                  for(int jj(0); jj<3; ++jj){
                    for(int ii(0); ii<3; ++ii){
                      AtA[ii][jj] = 0.0;
                    }
                    Atb[jj] = 0.0;
                  }

                  {
                    for(int lc(0); lc<27; ++lc){
                      AtA[0][0] += A[lc][0]* A[lc][0];
                      AtA[0][1] += A[lc][0]* A[lc][1];
                      AtA[0][2] += A[lc][0]* A[lc][2];
                      AtA[1][1] += A[lc][1]* A[lc][1];
                      AtA[1][2] += A[lc][1]* A[lc][2];
                      AtA[2][2] += A[lc][2]* A[lc][2];

                      Atb[0] += A[lc][0]*du[lc];
                      Atb[1] += A[lc][1]*du[lc];
                      Atb[2] += A[lc][2]*du[lc];
                    }
                  }

                  // Fill in symmetric
                  AtA[1][0] = AtA[0][1];
                  AtA[2][0] = AtA[0][2];
                  AtA[2][1] = AtA[1][2];

                  amrex::Real detAtA =
                    AtA[0][0]*(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
                    AtA[0][1]*(AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
                    AtA[0][2]*(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

                  // X direction
                  if(flags_sm[tid_x][tid_y][tid_z].isSingleValued() or
                     (flags_sm[tid_x-1][tid_y][tid_z].isSingleValued() or
                      (not flags_sm[tid_x][tid_y][tid_z].isConnected(-1,0,0))) or
                     (flags_sm[tid_x+1][tid_y][tid_z].isSingleValued() or
                      (not flags_sm[tid_x][tid_y][tid_z].isConnected( 1,0,0)))) {

                    amrex::Real detAtA_x =
                      Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
                      AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
                      AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );

                    // Slope at centroid of (i,j,k)
                    xs_fab(i,j,k,slopes_comp+n) = detAtA_x / detAtA;
                  }
                  else {
                    Real du_xl = 2.0*(state_sm[tid_x  ][tid_y][tid_z][n] - state_sm[tid_x-1][tid_y][tid_z][n]);
                    Real du_xr = 2.0*(state_sm[tid_x+1][tid_y][tid_z][n] - state_sm[tid_x  ][tid_y][tid_z][n]);
                    Real du_xc = 0.5*(state_sm[tid_x+1][tid_y][tid_z][n] - state_sm[tid_x-1][tid_y][tid_z][n]);

                    Real xslope = amrex::min(std::abs(du_xl), std::abs(du_xc), std::abs(du_xr));
                    xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                    xs_fab(i,j,k,slopes_comp+n) = (du_xc > 0.0) ? xslope : -xslope;
                  }

                  // Y direction
                  if(flags_sm[tid_x][tid_y][tid_z].isSingleValued() or
                     (flags_sm[tid_x][tid_y-1][tid_z].isSingleValued() or
                      (not flags_sm[tid_x][tid_y][tid_z].isConnected(0,-1,0))) or
                     (flags_sm[tid_x][tid_y+1][tid_z].isSingleValued() or
                      (not flags_sm[tid_x][tid_y][tid_z].isConnected(0, 1,0)))) {

                    amrex::Real detAtA_y =
                      AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
                      Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
                      AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);

                    // Slope at centroid of (i,j,k)
                    ys_fab(i,j,k,slopes_comp+n) = detAtA_y / detAtA;
                  }
                  else {
                    Real du_yl = 2.0*(state_sm[tid_x][tid_y  ][tid_z][n] - state_sm[tid_x][tid_y-1][tid_z][n]);
                    Real du_yr = 2.0*(state_sm[tid_x][tid_y+1][tid_z][n] - state_sm[tid_x][tid_y  ][tid_z][n]);
                    Real du_yc = 0.5*(state_sm[tid_x][tid_y+1][tid_z][n] - state_sm[tid_x][tid_y-1][tid_z][n]);

                    Real yslope = amrex::min(std::abs(du_yl), std::abs(du_yc), std::abs(du_yr));
                    yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                    ys_fab(i,j,k,slopes_comp+n) = (du_yc > 0.0) ? yslope : -yslope;
                  }

                  // Z direction
                  if(flags_sm[tid_x][tid_y][tid_z].isSingleValued() or
                     (flags_sm[tid_x][tid_y][tid_z-1].isSingleValued() or
                      (not flags_sm[tid_x][tid_y][tid_z].isConnected(0,0,-1))) or
                     (flags_sm[tid_x][tid_y][tid_z+1].isSingleValued() or
                      (not flags_sm[tid_x][tid_y][tid_z].isConnected(0,0, 1)))) {

                    amrex::Real detAtA_z =
                      AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
                      AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
                      Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

                    zs_fab(i,j,k,slopes_comp+n) = detAtA_z / detAtA;
                  }
                  else {
                    Real du_zl = 2.0*(state_sm[tid_x][tid_y][tid_z  ][n] - state_sm[tid_x][tid_y][tid_z-1][n]);
                    Real du_zr = 2.0*(state_sm[tid_x][tid_y][tid_z+1][n] - state_sm[tid_x][tid_y][tid_z  ][n]);
                    Real du_zc = 0.5*(state_sm[tid_x][tid_y][tid_z+1][n] - state_sm[tid_x][tid_y][tid_z-1][n]);

                    Real zslope = amrex::min(std::abs(du_zl), std::abs(du_zc), std::abs(du_zr));
                    zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                    zs_fab(i,j,k,slopes_comp+n) = (du_zc > 0.0) ? zslope : -zslope;
                  }
                }
              }
            }
          });
        }
#undef BLOCK_SIZE
#endif
        if(not kernel_run_with_launch_global) {
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

                     amrex::Real A[27][3];
                     amrex::Real du[27];

                     {
                       int lc=0;
                       for(int kk(-1); kk<=1; kk++){
                         for(int jj(-1); jj<=1; jj++){
                           for(int ii(-1); ii<=1; ii++){
                             if( flag_fab(i,j,k).isConnected(ii,jj,kk) and
                                 not (ii==0 and jj==0 and kk==0)) {

                               // Not multiplying by dx to be consistent with how the
                               // slope is stored. Also not including the global shift
                               // wrt plo or i,j,k. We only need relative distance.

                               A[lc][0] = ii + ccent_fab(i+ii,j+jj,k+kk,0) - ccent_fab(i,j,k,0);
                               A[lc][1] = jj + ccent_fab(i+ii,j+jj,k+kk,1) - ccent_fab(i,j,k,1);
                               A[lc][2] = kk + ccent_fab(i+ii,j+jj,k+kk,2) - ccent_fab(i,j,k,2);

                               du[lc] = state_fab(i+ii,j+jj,k+kk,n) - state_fab(i,j,k,n);

                             } else {

                               A[lc][0] = 0.0;
                               A[lc][1] = 0.0;
                               A[lc][2] = 0.0;

                               du[lc] = 0.0;
                             }

                             lc++;
                           }
                         }
                       }
                     }

                     amrex::Real AtA[3][3];
                     amrex::Real Atb[3];

                     for(int jj(0); jj<3; ++jj){
                       for(int ii(0); ii<3; ++ii){
                         AtA[ii][jj] = 0.0;
                       }
                       Atb[jj] = 0.0;
                     }


                     {

                       for(int lc(0); lc<27; ++lc){
                         AtA[0][0] += A[lc][0]* A[lc][0];
                         AtA[0][1] += A[lc][0]* A[lc][1];
                         AtA[0][2] += A[lc][0]* A[lc][2];
                         AtA[1][1] += A[lc][1]* A[lc][1];
                         AtA[1][2] += A[lc][1]* A[lc][2];
                         AtA[2][2] += A[lc][2]* A[lc][2];

                         Atb[0] += A[lc][0]*du[lc];
                         Atb[1] += A[lc][1]*du[lc];
                         Atb[2] += A[lc][2]*du[lc];
                       }
                     }

                     // Fill in symmetric
                     AtA[1][0] = AtA[0][1];
                     AtA[2][0] = AtA[0][2];
                     AtA[2][1] = AtA[1][2];


                     amrex::Real detAtA =
                       AtA[0][0]*(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
                       AtA[0][1]*(AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
                       AtA[0][2]*(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);



                     // X direction
                     if( flag_fab(i  ,j,k).isSingleValued() or
                        (flag_fab(i-1,j,k).isSingleValued() or not flag_fab(i,j,k).isConnected(-1,0,0)) or
                        (flag_fab(i+1,j,k).isSingleValued() or not flag_fab(i,j,k).isConnected( 1,0,0))) {

                       amrex::Real detAtA_x =
                         Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
                         AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
                         AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );

                       // Slope at centroid of (i,j,k)
                       xs_fab(i,j,k,slopes_comp+n) = detAtA_x / detAtA;

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

                       amrex::Real detAtA_y =
                         AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
                         Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
                         AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);

                       // Slope at centroid of (i,j,k)
                       ys_fab(i,j,k,slopes_comp+n) = detAtA_y / detAtA;

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

                       amrex::Real detAtA_z =
                         AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
                         AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
                         Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

                       zs_fab(i,j,k,slopes_comp+n) = detAtA_z / detAtA;

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
             }
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
