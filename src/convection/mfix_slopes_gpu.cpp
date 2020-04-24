#ifdef AMREX_USE_CUDA

#include <mfix.H>

namespace {
  struct get_index
  {
    const unsigned stride_i, stride_j, stride_k;

    AMREX_GPU_HOST_DEVICE
    get_index(unsigned dim_x, unsigned dim_y, unsigned dim_z, unsigned dim_comp)
      : stride_i(dim_y*dim_z*dim_comp)
      , stride_j(dim_z*dim_comp)
      , stride_k(dim_comp)
    {}

    AMREX_GPU_HOST_DEVICE
    unsigned operator() (unsigned i, unsigned j, unsigned k, unsigned n=0) const
    {
      return i*stride_i + j*stride_j + k*stride_k + n;
    }
  };
}

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

        dim3 Block(4,4,4);
        dim3 Grid(std::ceil(bx.length()[0] / float(Block.x)),
                  std::ceil(bx.length()[1] / float(Block.y)),
                  std::ceil(bx.length()[2] / float(Block.z)));

        const unsigned sm_size =
          ncomp*sizeof(Real)*(Block.x+2)*(Block.y+2)*(Block.z+2) + // ccent_sm
          ncomp*sizeof(Real)*(Block.x+2)*(Block.y+2)*(Block.z+2) + // state_sm
          sizeof(EBCellFlag)*(Block.x+2)*(Block.y+2)*(Block.z+2);  // flags_sm

        if(sm_size > Gpu::Device::sharedMemPerBlock()) {
          amrex::Abort("Exceeding GPU shared memory in mfix_slopes_gpu.cpp");
        }

        amrex::IntVect bx_lo(bx.loVect());
        amrex::IntVect bx_hi(bx.hiVect());

        amrex::launch_global<<<Grid,Block,sm_size,Gpu::gpuStream()>>>(
          [state_fab,xs_fab,ys_fab,zs_fab,slopes_comp,flag_fab,ccent_fab,ncomp,
           bx_lo,bx_hi] AMREX_GPU_DEVICE () noexcept
        {
          int tid_x = threadIdx.x;
          int tid_y = threadIdx.y;
          int tid_z = threadIdx.z;

          const int stride_x = blockIdx.x*blockDim.x + bx_lo[0];
          const int stride_y = blockIdx.y*blockDim.y + bx_lo[1];
          const int stride_z = blockIdx.z*blockDim.z + bx_lo[2];

          int i = tid_x + stride_x;
          int j = tid_y + stride_y;
          int k = tid_z + stride_z;

          get_index idx_sca(blockDim.x+2,blockDim.y+2,blockDim.z+2,1);
          get_index idx_vec(blockDim.x+2,blockDim.y+2,blockDim.z+2,ncomp);

          extern __shared__ Real sm[];

          const unsigned BlocksSize = (blockDim.x+2)*(blockDim.y+2)*(blockDim.z+2);

          Real* ccent_sm = sm;
          Real* state_sm = (Real*)&ccent_sm[ncomp*BlocksSize];
          EBCellFlag* flags_sm = (EBCellFlag*)&state_sm[ncomp*BlocksSize];

          if((i <= bx_hi[0]) and (j <= bx_hi[1]) and (k <= bx_hi[2]))
          {
            tid_x += 1; tid_y += 1; tid_z += 1;

            int loc_x[3] = {0, tid_x, amrex::min((int)blockDim.x, bx_hi[0]-stride_x+1)+1};
            int loc_y[3] = {0, tid_y, amrex::min((int)blockDim.y, bx_hi[1]-stride_y+1)+1};
            int loc_z[3] = {0, tid_z, amrex::min((int)blockDim.z, bx_hi[2]-stride_z+1)+1};

            int glob_x[3] = {stride_x-1, i, amrex::min(stride_x+(int)blockDim.x, bx_hi[0]+1)};
            int glob_y[3] = {stride_y-1, j, amrex::min(stride_y+(int)blockDim.y, bx_hi[1]+1)};
            int glob_z[3] = {stride_z-1, k, amrex::min(stride_z+(int)blockDim.z, bx_hi[2]+1)};

            for (int ii(0); ii < 3; ii++)
            for (int jj(0); jj < 3; jj++)
            for (int kk(0); kk < 3; kk++) {
              int t_x(loc_x[ii]);
              int t_y(loc_y[jj]);
              int t_z(loc_z[kk]);

              int g_x(glob_x[ii]);
              int g_y(glob_y[jj]);
              int g_z(glob_z[kk]);

              for (int n(0); n < ncomp; n++) {
                ccent_sm[idx_vec(t_x,t_y,t_z,n)] = ccent_fab(g_x,g_y,g_z,n);
                state_sm[idx_vec(t_x,t_y,t_z,n)] = state_fab(g_x,g_y,g_z,n);
              }

              flags_sm[idx_sca(t_x,t_y,t_z)] = flag_fab(g_x,g_y,g_z);
            }

            __syncthreads();

            for (int n(0); n < ncomp; n++)
            {
              if (flags_sm[idx_sca(tid_x,tid_y,tid_z)].isCovered())
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
                    if((ii!=0 or jj!=0 or kk!=0) and
                       flags_sm[idx_sca(tid_x,tid_y,tid_z)].isConnected(ii,jj,kk)) {
                      // Not multiplying by dx to be consistent with how the
                      // slope is stored. Also not including the global shift
                      // wrt plo or i,j,k. We only need relative distance.
                      A[lc][0] =
                        ii + ccent_sm[idx_vec(tid_x+ii,tid_y+jj,tid_z+kk,0)] -
                        ccent_sm[idx_vec(tid_x,tid_y,tid_z,0)];
                      A[lc][1] =
                        jj + ccent_sm[idx_vec(tid_x+ii,tid_y+jj,tid_z+kk,1)] -
                        ccent_sm[idx_vec(tid_x,tid_y,tid_z,1)];
                      A[lc][2] =
                        kk + ccent_sm[idx_vec(tid_x+ii,tid_y+jj,tid_z+kk,2)]-
                        ccent_sm[idx_vec(tid_x,tid_y,tid_z,2)];

                      du[lc] = state_sm[idx_vec(tid_x+ii,tid_y+jj,tid_z+kk,n)] -
                        state_sm[idx_vec(tid_x,tid_y,tid_z,n)];
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
                if(flags_sm[idx_sca(tid_x,tid_y,tid_z)].isSingleValued() or
                   (flags_sm[idx_sca(tid_x-1,tid_y,tid_z)].isSingleValued() or
                    (not flags_sm[idx_sca(tid_x,tid_y,tid_z)].isConnected(-1,0,0))) or
                   (flags_sm[idx_sca(tid_x+1,tid_y,tid_z)].isSingleValued() or
                    (not flags_sm[idx_sca(tid_x,tid_y,tid_z)].isConnected( 1,0,0)))) {

                  amrex::Real detAtA_x =
                    Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
                    AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
                    AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );

                  // Slope at centroid of (i,j,k)
                  xs_fab(i,j,k,slopes_comp+n) = detAtA_x / detAtA;
                }
                else {
                  Real du_xl = 2.0*(state_sm[idx_vec(tid_x,tid_y,tid_z,n)] -
                      state_sm[idx_vec(tid_x-1,tid_y,tid_z,n)]);
                  Real du_xr = 2.0*(state_sm[idx_vec(tid_x+1,tid_y,tid_z,n)] -
                      state_sm[idx_vec(tid_x,tid_y,tid_z,n)]);
                  Real du_xc = 0.5*(state_sm[idx_vec(tid_x+1,tid_y,tid_z,n)] -
                      state_sm[idx_vec(tid_x-1,tid_y,tid_z,n)]);

                  Real xslope = amrex::min(std::abs(du_xl), std::abs(du_xc), std::abs(du_xr));
                  xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                  xs_fab(i,j,k,slopes_comp+n) = (du_xc > 0.0) ? xslope : -xslope;
                }

                // Y direction
                if(flags_sm[idx_sca(tid_x,tid_y,tid_z)].isSingleValued() or
                   (flags_sm[idx_sca(tid_x,tid_y-1,tid_z)].isSingleValued() or
                    (not flags_sm[idx_sca(tid_x,tid_y,tid_z)].isConnected(0,-1,0))) or
                   (flags_sm[idx_sca(tid_x,tid_y+1,tid_z)].isSingleValued() or
                    (not flags_sm[idx_sca(tid_x,tid_y,tid_z)].isConnected(0, 1,0)))) {

                  amrex::Real detAtA_y =
                    AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
                    Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
                    AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);

                  // Slope at centroid of (i,j,k)
                  ys_fab(i,j,k,slopes_comp+n) = detAtA_y / detAtA;
                }
                else {
                  Real du_yl = 2.0*(state_sm[idx_vec(tid_x,tid_y,tid_z,n)] -
                      state_sm[idx_vec(tid_x,tid_y-1,tid_z,n)]);
                  Real du_yr = 2.0*(state_sm[idx_vec(tid_x,tid_y+1,tid_z,n)] -
                      state_sm[idx_vec(tid_x,tid_y,tid_z,n)]);
                  Real du_yc = 0.5*(state_sm[idx_vec(tid_x,tid_y+1,tid_z,n)] -
                      state_sm[idx_vec(tid_x,tid_y-1,tid_z,n)]);

                  Real yslope = amrex::min(std::abs(du_yl), std::abs(du_yc), std::abs(du_yr));
                  yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                  ys_fab(i,j,k,slopes_comp+n) = (du_yc > 0.0) ? yslope : -yslope;
                }

                // Z direction
                if(flags_sm[idx_sca(tid_x,tid_y,tid_z)].isSingleValued() or
                   (flags_sm[idx_sca(tid_x,tid_y,tid_z-1)].isSingleValued() or
                    (not flags_sm[idx_sca(tid_x,tid_y,tid_z)].isConnected(0,0,-1))) or
                   (flags_sm[idx_sca(tid_x,tid_y,tid_z+1)].isSingleValued() or
                    (not flags_sm[idx_sca(tid_x,tid_y,tid_z)].isConnected(0,0, 1)))) {

                  amrex::Real detAtA_z =
                    AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
                    AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
                    Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

                  zs_fab(i,j,k,slopes_comp+n) = detAtA_z / detAtA;
                }
                else {
                  Real du_zl = 2.0*(state_sm[idx_vec(tid_x,tid_y,tid_z,n)] -
                      state_sm[idx_vec(tid_x,tid_y,tid_z-1,n)]);
                  Real du_zr = 2.0*(state_sm[idx_vec(tid_x,tid_y,tid_z+1,n)] -
                      state_sm[idx_vec(tid_x,tid_y,tid_z,n)]);
                  Real du_zc = 0.5*(state_sm[idx_vec(tid_x,tid_y,tid_z+1,n)] -
                      state_sm[idx_vec(tid_x,tid_y,tid_z-1,n)]);

                  Real zslope = amrex::min(std::abs(du_zl), std::abs(du_zc), std::abs(du_zr));
                  zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                  zs_fab(i,j,k,slopes_comp+n) = (du_zc > 0.0) ? zslope : -zslope;
                }
              }
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

      if(domain.smallEnd(0) >= bx.smallEnd(0) and domain.smallEnd(0) <= bx.bigEnd(0))
      {
        Box bx_x_lo(IntVect(domain.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)),
                    IntVect(domain.smallEnd(0), bx.bigEnd(1), bx.bigEnd(2)));

        amrex::ParallelFor(bx_x_lo, ncomp,
          [domain,flag_fab,ilo_ifab,state_fab,xs_fab,minf,slopes_comp]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          const Real state = state_fab(i,j,k,n);
          const Real state_mns = state_fab(i-1,j,k,n);
          const Real state_pls = state_fab(i+1,j,k,n);

          if (!flag_fab(i,j,k).isCovered() and ilo_ifab(i-1,j,k,0) == minf)
          {
            Real du_xl = 2.0*(state - state_mns);
            Real du_xr = 2.0*(state_pls - state);
            Real du_xc = (state_pls+3.0*state-4.0*state_mns)/3.0;

            Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
            xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
            xs_fab(i,j,k,slopes_comp+n) = (du_xc > 0.0) ? xslope : -xslope;
          }
        });
      }

      if(domain.bigEnd(0) >= bx.smallEnd(0) and domain.bigEnd(0) <= bx.bigEnd(0)) {
        Box bx_x_hi(IntVect(domain.bigEnd(0), bx.smallEnd(1), bx.smallEnd(2)),
                    IntVect(domain.bigEnd(0), bx.bigEnd(1), bx.bigEnd(2)));

        amrex::ParallelFor(bx_x_hi, ncomp,
          [domain,flag_fab,ihi_ifab,state_fab,xs_fab,minf,slopes_comp]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          const Real state = state_fab(i,j,k,n);
          const Real state_mns = state_fab(i-1,j,k,n);
          const Real state_pls = state_fab(i+1,j,k,n);

          if (!flag_fab(i,j,k).isCovered() and ihi_ifab(i+1,j,k,0) == minf)
          {
            Real du_xl = 2.0*(state - state_mns);
            Real du_xr = 2.0*(state_pls - state);
            Real du_xc = -(state_mns+3.0*state-4.0*state_pls)/3.0;

            Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
            xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
            xs_fab(i,j,k,slopes_comp+n) = (du_xc > 0.0) ? xslope : -xslope;
          }
        });
      }

      if(domain.smallEnd(1) >= bx.smallEnd(1) and domain.smallEnd(1) <= bx.bigEnd(1)) 
      {
        Box bx_y_lo(IntVect(bx.smallEnd(0), domain.smallEnd(1), bx.smallEnd(2)),
                    IntVect(bx.bigEnd(0), domain.smallEnd(1), bx.bigEnd(2)));

        amrex::ParallelFor(bx_y_lo, ncomp,
          [domain,flag_fab,jlo_ifab,state_fab,ys_fab,minf,slopes_comp]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          const Real state = state_fab(i,j,k,n);
          const Real state_mns = state_fab(i,j-1,k,n);
          const Real state_pls = state_fab(i,j+1,k,n);

          if (!flag_fab(i,j,k).isCovered() and jlo_ifab(i,j-1,k,0) == minf)
          {
            Real du_yl = 2.0*(state - state_mns);
            Real du_yr = 2.0*(state_pls - state);
            Real du_yc = (state_pls+3.0*state-4.0*state_mns)/3.0;

            Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
            yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
            ys_fab(i,j,k,slopes_comp+n) = (du_yc > 0.0) ? yslope : -yslope;
          }
        });
      }

      if(domain.bigEnd(1) >= bx.smallEnd(1) and domain.bigEnd(1) <= bx.bigEnd(1))
      {
        Box bx_y_hi(IntVect(bx.smallEnd(0), domain.bigEnd(1), bx.smallEnd(2)),
                    IntVect(bx.bigEnd(0), domain.bigEnd(1), bx.bigEnd(2)));

        amrex::ParallelFor(bx_y_hi, ncomp,
          [domain,flag_fab,jhi_ifab,state_fab,ys_fab,minf,slopes_comp]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          const Real state = state_fab(i,j,k,n);
          const Real state_mns = state_fab(i,j-1,k,n);
          const Real state_pls = state_fab(i,j+1,k,n);

          if (!flag_fab(i,j,k).isCovered() and jhi_ifab(i,j+1,k,0) == minf)
          {
            Real du_yl = 2.0*(state - state_mns);
            Real du_yr = 2.0*(state_pls - state);
            Real du_yc = -(state_mns+3.0*state-4.0*state_pls)/3.0;

            Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
            yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
            ys_fab(i,j,k,slopes_comp+n) = (du_yc > 0.0) ? yslope : -yslope;
          }
        });
      }

      if(domain.smallEnd(2) >= bx.smallEnd(2) and domain.smallEnd(2) <= bx.bigEnd(2))
      {
        Box bx_z_lo(IntVect(bx.smallEnd(0), bx.smallEnd(1), domain.smallEnd(2)),
                    IntVect(bx.bigEnd(0), bx.bigEnd(1), domain.smallEnd(2)));

        amrex::ParallelFor(bx_z_lo, ncomp,
          [domain,flag_fab,klo_ifab,state_fab,zs_fab,minf,slopes_comp]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          const Real state = state_fab(i,j,k,n);
          const Real state_mns = state_fab(i,j,k-1,n);
          const Real state_pls = state_fab(i,j,k+1,n);

          if (!flag_fab(i,j,k).isCovered() and klo_ifab(i,j,k-1,0) == minf)
          {
            Real du_zl = 2.0*(state - state_mns);
            Real du_zr = 2.0*(state_pls - state);
            Real du_zc = (state_pls+3.0*state-4.0*state_mns)/3.0;

            Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
            zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
            zs_fab(i,j,k,slopes_comp+n) = (du_zc > 0.0) ? zslope : -zslope;
          }
        });
      }

      if(domain.bigEnd(2) >= bx.smallEnd(2) and domain.bigEnd(2) <= bx.bigEnd(2))
      {
        Box bx_z_hi(IntVect(bx.smallEnd(0), bx.smallEnd(1), domain.bigEnd(2)),
                    IntVect(bx.bigEnd(0), bx.bigEnd(1), domain.bigEnd(2)));

        amrex::ParallelFor(bx_z_hi, ncomp,
          [domain,flag_fab,khi_ifab,state_fab,zs_fab,minf,slopes_comp]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          const Real state = state_fab(i,j,k,n);
          const Real state_mns = state_fab(i,j,k-1,n);
          const Real state_pls = state_fab(i,j,k+1,n);

          if (!flag_fab(i,j,k).isCovered() and khi_ifab(i,j,k+1,0) == minf)
          {
            Real du_zl = 2.0*(state - state_mns);
            Real du_zr = 2.0*(state_pls - state);
            Real du_zc = -(state_mns+3.0*state-4.0*state_pls)/3.0;

            Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
            zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
            zs_fab(i,j,k,slopes_comp+n) = (du_zc > 0.0) ? zslope : -zslope;
          }
        });
      }
    } // not covered
  } // MFIter

  xslopes_in[lev]->FillBoundary(geom[lev].periodicity());
  yslopes_in[lev]->FillBoundary(geom[lev].periodicity());
  zslopes_in[lev]->FillBoundary(geom[lev].periodicity());
}

#endif
