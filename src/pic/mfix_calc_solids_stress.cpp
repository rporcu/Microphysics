#include <AMReX_BC_TYPES.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>

#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_pic_K.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_K.H>
#include <mfix_mf_helpers.H>

void mfix::MFIX_CalcSolidsStress (amrex::Vector< amrex::MultiFab* >& ep_s_in,
                                  amrex::Vector< amrex::MultiFab* >& avg_prop_in)
{
  BL_PROFILE("mfix::MFIX_CalcSolidsStress()");

#define DO_VISMF 0

  // We copy the value inside the domain to the outside to avoid
  // unphysical volume fractions.
  const int dir_bc_in = 2;
  mfix_set_epg_bcs(ep_s_in, dir_bc_in);

  const amrex::Real covered_val = 9.8765e300;

  for( int lev(0); lev < nlev; lev ++)
  {
    const int interp_comp = 4;  // Four components (3 grad_tau + 1 ep_s)

    const BoxArray&            pba = pc->ParticleBoxArray(lev);
    const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

    // Solids stress gradient. Note that it has an extra component to
    // store a copy of the solids volume fraction for interpolation.
    MultiFab grad_tau;
    grad_tau.define(pba, pdm, interp_comp, 1, MFInfo(), *particle_ebfactory[lev]);
    grad_tau.setVal(0.0, 0, interp_comp, grad_tau.nGrow());

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(grad_tau.Factory());
    const auto& flags = factory.getMultiEBCellFlagFab();
    const auto& cellcent = factory.getCentroid();

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom[lev].CellSizeArray();
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom[lev].InvCellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(grad_tau, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      // Tilebox
      const Box& bx = mfi.tilebox();

      Array4<const Real> const& ep_s_arr = ep_s_in[lev]->const_array(mfi);
      Array4<Real> const& grad_tau_arr = grad_tau.array(mfi);

      const Real Ps = PIC::Ps;
      const Real beta = PIC::beta;
      const Real ep_cp = PIC::ep_cp;
      const Real small_number = PIC::small_number;

      // const auto& area  = particle_ebfactory.getAreaFrac();


      const auto& flagfab = flags[mfi];

      // Fully covered FAB
      if (flagfab.getType(amrex::grow(bx,0)) == FabType::covered )
      {
        amrex::ParallelFor(bx,[grad_tau_arr, covered_val]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
           grad_tau_arr(i,j,k,0) = covered_val;
           grad_tau_arr(i,j,k,1) = covered_val;
           grad_tau_arr(i,j,k,2) = covered_val;
           grad_tau_arr(i,j,k,3) = covered_val;
        });
      }

      // No cut cells in this FAB
      else if (flagfab.getType(amrex::grow(bx,1)) == FabType::regular )
      {

        amrex::ParallelFor(bx,
          [ep_s_arr,grad_tau_arr,Ps,beta,ep_cp,small_number,dxinv]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

          grad_tau_arr(i,j,k,3) = ep_s_arr(i,j,k);

          amrex::Real tau_lo;
          amrex::Real tau_hi;

          // x-direction
          tau_lo = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i-1,j  ,k  ));
          tau_hi = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i+1,j  ,k  ));

          grad_tau_arr(i,j,k,0) = 0.5 * dxinv[0] * (tau_hi - tau_lo);

          // Y-direction
          tau_lo = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i  ,j-1,k  ));
          tau_hi = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i  ,j+1,k  ));

          grad_tau_arr(i,j,k,1) = 0.5 * dxinv[1] * (tau_hi - tau_lo);

          // Z-direction
          tau_lo = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i  ,j  ,k-1));
          tau_hi = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i  ,j  ,k+1));

          grad_tau_arr(i,j,k,2) = 0.5 * dxinv[2] * (tau_hi - tau_lo);

        });

      }

      // A mix of regular and cut cells
      else
      {
        const auto& flagsarr = flagfab.array();
        Array4<Real const> const& ccent = cellcent.const_array(mfi);

        amrex::ParallelFor(bx,
          [ep_s_arr,grad_tau_arr,Ps,beta,ep_cp,small_number,flagsarr,
          covered_val,ccent,dx, dxinv]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (flagsarr(i,j,k).isCovered())
          {
            grad_tau_arr(i,j,k,0) = covered_val;
            grad_tau_arr(i,j,k,1) = covered_val;
            grad_tau_arr(i,j,k,2) = covered_val;
            grad_tau_arr(i,j,k,3) = covered_val;
          }
          else
          {

            grad_tau_arr(i,j,k,3) = ep_s_arr(i,j,k);

            amrex::IntVect do_least_squares = {0,0,0};

            if (flagsarr(i,j,k).isSingleValued()) {
              do_least_squares = {1,1,1};

            } else {
              // Checking x-dir
              if (not flagsarr(i,j,k).isConnected(-1,0,0) or flagsarr(i-1,j  ,k  ).isSingleValued() or
                  not flagsarr(i,j,k).isConnected( 1,0,0) or flagsarr(i+1,j  ,k  ).isSingleValued()) {
                do_least_squares[0] = 1;

              } else {

                const amrex::Real tau_lo = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i-1,j  ,k  ));
                const amrex::Real tau_hi = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i+1,j  ,k  ));

                grad_tau_arr(i,j,k,0) = 0.5 * dxinv[0] * (tau_hi - tau_lo);
              }

              // Checking y-dir
              if (not flagsarr(i,j,k).isConnected(0,-1,0) or flagsarr(i  ,j-1,k  ).isSingleValued() or
                  not flagsarr(i,j,k).isConnected(0, 1,0) or flagsarr(i  ,j+1,k  ).isSingleValued()) {
                do_least_squares[1] = 1;

              } else {

                const amrex::Real tau_lo = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i  ,j-1,k  ));
                const amrex::Real tau_hi = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i  ,j+1,k  ));

                grad_tau_arr(i,j,k,1) = 0.5 * dxinv[1] * (tau_hi - tau_lo);
              }

              // Checking z-dir
              if (not flagsarr(i,j,k).isConnected(0,0,-1) or flagsarr(i  ,j  ,k-1).isSingleValued() or
                  not flagsarr(i,j,k).isConnected(0,0, 1) or flagsarr(i  ,j  ,k+1).isSingleValued()) {
                do_least_squares[2] = 1;

              } else {

                const amrex::Real tau_lo = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i  ,j  ,k-1));
                const amrex::Real tau_hi = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i  ,j  ,k+1));

                grad_tau_arr(i,j,k,2) = 0.5 * dxinv[2] * (tau_hi - tau_lo);
              }
            }

            if( do_least_squares[0] == 1 or
                do_least_squares[1] == 1 or
                do_least_squares[2] == 1 ) {

              GpuArray< GpuArray <Real, 3>, 27> A;
              GpuArray< Real, 27> du;

              {
                const amrex::Real tau_lo = solids_pressure(Ps, beta, ep_cp, small_number, ep_s_arr(i,j,k));

                int lc=0;
                for(int kk(-1); kk<=1; kk++){
                  for(int jj(-1); jj<=1; jj++){
                    for(int ii(-1); ii<=1; ii++){
                      if( flagsarr(i,j,k).isConnected(ii,jj,kk) ) {

                        A[lc][0] = (ii + ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0))*dx[0];
                        A[lc][1] = (jj + ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1))*dx[1];
                        A[lc][2] = (kk + ccent(i+ii,j+jj,k+kk,2) - ccent(i,j,k,2))*dx[2];

                        const amrex::Real tau_hi = solids_pressure(Ps, beta, ep_cp,
                            small_number, ep_s_arr(i+ii,j+jj,k+kk));

                        du[lc] = tau_hi - tau_lo;

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

              GpuArray< GpuArray< Real, 3>, 3> AtA;
              GpuArray< Real, 3> Atb;

              for(int jj(0); jj<3; ++jj){
                for(int ii(0); ii<3; ++ii){
                  AtA[ii][jj] = 0.0;
                }
                Atb[jj] = 0.0;
              }

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

              // Fill in symmetric
              AtA[1][0] = AtA[0][1];
              AtA[2][0] = AtA[0][2];
              AtA[2][1] = AtA[1][2];

              const amrex::Real detAtA =
                AtA[0][0]*(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
                AtA[0][1]*(AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
                AtA[0][2]*(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

              AMREX_ASSERT_WITH_MESSAGE(detAtA > 0.,
                "Negative determinate in solids stress term");

              if( do_least_squares[0] == 1 ) {

                const amrex::Real detAtA_x =
                  Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
                  AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
                  AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );

                grad_tau_arr(i,j,k,0) = detAtA_x / detAtA;
              }

              if( do_least_squares[1] == 1 ) {

                amrex::Real detAtA_y =
                  AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
                  Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
                  AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);

                grad_tau_arr(i,j,k,1) = detAtA_y / detAtA;
              }

              if( do_least_squares[2] == 1 ) {

                amrex::Real detAtA_z =
                  AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
                  AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
                  Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

                grad_tau_arr(i,j,k,2) = detAtA_z / detAtA;
              }

            } // do least squares
          }
        });
      }
    }

    grad_tau.FillBoundary(geom[lev].periodicity());


#if DO_VISMF > 0
    {
      EB_set_covered(grad_tau, 0, 4, 1, 0.0);
      const std::string& tfile = "grad_tau";
      VisMF::Write(grad_tau, tfile);
    }
#endif

    // We now have the solids stress on the field. We now need to interpolate
    // this to the parcel positions. We are going to store this for the
    // parcel in "omega" in the particle container.

    using MFIXParIter = MFIXParticleContainer::MFIXParIter;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {

      const auto dx_array  = geom[lev].CellSizeArray();
      const auto dxi_array = geom[lev].InvCellSizeArray();
      const auto plo_array = geom[lev].ProbLoArray();

      const amrex::RealVect dxv(dx_array[0], dx_array[1], dx_array[2]);
      const amrex::RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
      const amrex::RealVect plo(plo_array[0], plo_array[1], plo_array[2]);

      const auto cellcent =  &(particle_ebfactory[lev]->getCentroid());
      const auto bndrycent = &(particle_ebfactory[lev]->getBndryCent());
      const auto areafrac =  particle_ebfactory[lev]->getAreaFrac();

      for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        const int np = particles.size();

        Box bx = pti.tilebox ();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  ep_s_fab = static_cast<EBFArrayBox const&>((*ep_s_in[lev])[pti]);
        const EBCellFlagFab&  flags = ep_s_fab.getEBCellFlagFab();

        const auto& grad_tau_array = grad_tau.const_array(pti);
        const auto& flags_array = flags.array();

        if (flags.getType(amrex::grow(bx,1)) == FabType::covered)
        {
          // We shouldn't have this case but if we do -- zero the stress.
          amrex::ParallelFor(np, [pstruct]
            AMREX_GPU_DEVICE (int ip) noexcept
            {
              MFIXParticleContainer::ParticleType& particle = pstruct[ip];
              particle.rdata(realData::oneOverI) = 0.0;
              particle.rdata(realData::omegax) = 0.0;
              particle.rdata(realData::omegay) = 0.0;
              particle.rdata(realData::omegaz) = 0.0;
            });
        }
        else if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
        {
          amrex::ParallelFor(np,
            [pstruct,grad_tau_array,plo,dxi]
            AMREX_GPU_DEVICE (int ip) noexcept
          {

            MFIXParticleContainer::ParticleType& particle = pstruct[ip];

            // Local array storing interpolated values
            GpuArray<Real, interp_comp> interp_loc;

            // Remember we cheated an stored ep_s as a 4th component
            trilinear_interp(particle.pos(), interp_loc.data(),
              grad_tau_array, plo, dxi, interp_comp);

            const Real ep_s_loc = interp_loc[3];

            // Store local solids volume fraction in unused array location.
            particle.rdata(realData::oneOverI) = ep_s_loc;

            if(ep_s_loc > 1.e-4) {

              const amrex::Real inv_rops = 1.0 / (ep_s_loc*particle.rdata(realData::density));

              particle.rdata(realData::omegax) = interp_loc[0]*inv_rops;
              particle.rdata(realData::omegay) = interp_loc[1]*inv_rops;
              particle.rdata(realData::omegaz) = interp_loc[2]*inv_rops;

            } else {

              particle.rdata(realData::omegax) = 0.0;
              particle.rdata(realData::omegay) = 0.0;
              particle.rdata(realData::omegaz) = 0.0;

            }

          });
        }
        else // FAB not all regular
        {

          // Cell centroids
          const auto& ccent_fab = cellcent->array(pti);
          // Centroid of EB
          const auto& bcent_fab = bndrycent->array(pti);
          // Area fractions
          const auto& apx_fab = areafrac[0]->array(pti);
          const auto& apy_fab = areafrac[1]->array(pti);
          const auto& apz_fab = areafrac[2]->array(pti);


          amrex::ParallelFor(np,
            [pstruct,grad_tau_array,flags_array,plo,dxi,dxv,ccent_fab,
             bcent_fab, apx_fab, apy_fab, apz_fab]
            AMREX_GPU_DEVICE (int pid) noexcept
          {
            // Local array storing interpolated values
            GpuArray<Real, interp_comp> interp_loc;

            MFIXParticleContainer::ParticleType& particle = pstruct[pid];

            // Cell containing particle centroid
            const int ip = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
            const int jp = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
            const int kp = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

            if(flags_array(ip,jp,kp).isCovered())
            {
              particle.rdata(realData::omegax) = 0.0;
              particle.rdata(realData::omegay) = 0.0;
              particle.rdata(realData::omegaz) = 0.0;

            // Cut or regular cell and none of the cells in the stencil is covered
            // (Note we can't assume regular cell has no covered cells in the stencil
            //      because of the diagonal case)
            } else {

              // Upper cell in trilinear stencil
              const int i = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0] + 0.5));
              const int j = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1] + 0.5));
              const int k = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2] + 0.5));

              // All cells in the stencil are regular. Use
              // traditional trilinear interpolation
              if (flags_array(i-1,j-1,k-1).isRegular() and
                  flags_array(i  ,j-1,k-1).isRegular() and
                  flags_array(i-1,j  ,k-1).isRegular() and
                  flags_array(i  ,j  ,k-1).isRegular() and
                  flags_array(i-1,j-1,k  ).isRegular() and
                  flags_array(i  ,j-1,k  ).isRegular() and
                  flags_array(i-1,j  ,k  ).isRegular() and
                  flags_array(i  ,j  ,k  ).isRegular()) {

                  trilinear_interp(particle.pos(), interp_loc.data(),
                                    grad_tau_array, plo, dxi, interp_comp);

              // At least one of the cells in the stencil is cut or covered
              } else {

                const int scomp = 0;
                fe_interp(particle.pos(), ip, jp, kp, dxv, dxi, plo,
                          flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                          grad_tau_array, interp_loc.data(), interp_comp, scomp);

              } // Cut cell

              // Store local solids volume fraction in unused array location.
              const Real ep_s_loc = interp_loc[3];
              particle.rdata(realData::oneOverI) = ep_s_loc;

              if(ep_s_loc > 1.e-4) {

                const amrex::Real inv_rops = 1.0 / (ep_s_loc*particle.rdata(realData::density));
                particle.rdata(realData::omegax) = interp_loc[0]*inv_rops;
                particle.rdata(realData::omegay) = interp_loc[1]*inv_rops;
                particle.rdata(realData::omegaz) = interp_loc[2]*inv_rops;

              } else {
                particle.rdata(realData::omegax) = 0.0;
                particle.rdata(realData::omegay) = 0.0;
                particle.rdata(realData::omegaz) = 0.0;
              }
            }
          }); // part loop
        } // if box not all regular
      } // FAB not covered
    } // pti
  } // omp region

}
