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
                                  amrex::Vector< amrex::MultiFab* >& avg_prop_in,
                                  amrex::Vector< amrex::MultiFab* >& cost,
                                  std::string& knapsack_weight_type)
{
  BL_PROFILE("mfix::MFIX_CalcSolidsStress()");

  //const amrex::Real covered_val = 9.8765e300; UNUSED

  for( int lev(0); lev < nlev; lev ++)
  {

    const BoxArray&            pba = pc->ParticleBoxArray(lev);
    const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

    // Solids stress gradient. Note that it has an extra component to
    // store a copy of the solids volume fraction for interpolation.
    MultiFab Ps;
    Ps.define(pba, pdm, 1, 1, MFInfo(), *particle_ebfactory[lev]);
    Ps.setVal(-1.0e100, 0, 1, Ps.nGrow());

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(Ps.Factory());
    const auto& flags = factory.getMultiEBCellFlagFab();

    const Real Ps0 = PIC::Ps;
    const Real beta = PIC::beta;
    const Real ep_cp = PIC::ep_cp;
    const Real small_number = PIC::small_number;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Ps, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      // Tilebox
      const Box& bx = mfi.tilebox();

      Array4<const Real> const& ep_s_arr = ep_s_in[lev]->const_array(mfi);
      Array4<      Real> const& Ps_arr   = Ps.array(mfi);


      const auto& flagfab = flags[mfi];

      // Fully covered FAB
      if (flagfab.getType(amrex::grow(bx,0)) == FabType::covered ) {
        amrex::ParallelFor(bx,[Ps_arr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
           Ps_arr(i,j,k) = 0.0;
        });

      } else if (flagfab.getType(amrex::grow(bx,1)) == FabType::regular ) {

        amrex::ParallelFor(bx,
          [ep_s_arr,Ps_arr,Ps0,beta,ep_cp,small_number]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

          Ps_arr(i,j,k) = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));

        });

      } else {
        const auto& flagsarr = flagfab.array();

        amrex::ParallelFor(bx,
          [ep_s_arr,Ps_arr,Ps0,beta,ep_cp,small_number,flagsarr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (flagsarr(i,j,k).isCovered()) {
            Ps_arr(i,j,k) = 0.0;

          } else {

            Ps_arr(i,j,k) = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));

          }
        });
      }
    }

    // Set all values outside the domain to a crazy value to make
    // catching errors easier to find.
    Ps.setDomainBndry(-2.0e100,geom[lev]);

    // Calculate the pressure in flow boundary cells. This is needed
    // for inflow faces that may need to support the weight of the bed.

    Box domain(geom[lev].Domain());

    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());

    for (MFIter mfi(Ps, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      Array4<Real> const& ep_s_arr = ep_s_in[lev]->array(mfi);
      Array4<Real> const& Ps_arr = Ps.array(mfi);

      IntVect ps_lo(Ps[mfi].loVect());
      IntVect ps_hi(Ps[mfi].hiVect());

      const Box ps_bx(ps_lo, ps_hi);

      Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
      Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
      Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
      Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
      Array4<const int> const& bct_klo = bc_klo[lev]->array();
      Array4<const int> const& bct_khi = bc_khi[lev]->array();

      const int nlft = amrex::max(0, dom_lo[0]-ps_lo[0]);
      const int nbot = amrex::max(0, dom_lo[1]-ps_lo[1]);
      const int ndwn = amrex::max(0, dom_lo[2]-ps_lo[2]);

      const int nrgt = amrex::max(0, ps_hi[0]-dom_hi[0]);
      const int ntop = amrex::max(0, ps_hi[1]-dom_hi[1]);
      const int nup  = amrex::max(0, ps_hi[2]-dom_hi[2]);

      // Create InVects for following 2D Boxes
      IntVect bx_yz_lo_lo_2D(ps_lo), bx_yz_lo_hi_2D(ps_hi);
      IntVect bx_yz_hi_lo_2D(ps_lo), bx_yz_hi_hi_2D(ps_hi);
      IntVect bx_xz_lo_lo_2D(ps_lo), bx_xz_lo_hi_2D(ps_hi);
      IntVect bx_xz_hi_lo_2D(ps_lo), bx_xz_hi_hi_2D(ps_hi);
      IntVect bx_xy_lo_lo_2D(ps_lo), bx_xy_lo_hi_2D(ps_hi);
      IntVect bx_xy_hi_lo_2D(ps_lo), bx_xy_hi_hi_2D(ps_hi);

      // Fix lo and hi limits
      bx_yz_lo_lo_2D[0] = dom_lo[0]-1;
      bx_yz_lo_hi_2D[0] = dom_lo[0]-1;
      bx_yz_hi_lo_2D[0] = dom_hi[0]+1;
      bx_yz_hi_hi_2D[0] = dom_hi[0]+1;

      bx_xz_lo_lo_2D[1] = dom_lo[1]-1;
      bx_xz_lo_hi_2D[1] = dom_lo[1]-1;
      bx_xz_hi_lo_2D[1] = dom_hi[1]+1;
      bx_xz_hi_hi_2D[1] = dom_hi[1]+1;

      bx_xy_lo_lo_2D[2] = dom_lo[2]-1;
      bx_xy_lo_hi_2D[2] = dom_lo[2]-1;
      bx_xy_hi_lo_2D[2] = dom_hi[2]+1;
      bx_xy_hi_hi_2D[2] = dom_hi[2]+1;

      // Create 2D boxes for GPU loops
      const Box bx_yz_lo_2D(bx_yz_lo_lo_2D, bx_yz_lo_hi_2D);
      const Box bx_yz_hi_2D(bx_yz_hi_lo_2D, bx_yz_hi_hi_2D);

      const Box bx_xz_lo_2D(bx_xz_lo_lo_2D, bx_xz_lo_hi_2D);
      const Box bx_xz_hi_2D(bx_xz_hi_lo_2D, bx_xz_hi_hi_2D);

      const Box bx_xy_lo_2D(bx_xy_lo_lo_2D, bx_xy_lo_hi_2D);
      const Box bx_xy_hi_2D(bx_xy_hi_lo_2D, bx_xy_hi_hi_2D);

      const int minf = bc_list.get_minf();
      const int pinf = bc_list.get_pinf();
      const int pout = bc_list.get_pout();


      if (nlft > 0) {

        amrex::ParallelFor(bx_yz_lo_2D,
        [bct_ilo,dom_lo,minf,pinf,pout,Ps_arr,ep_s_arr,
         Ps0,beta,ep_cp, small_number]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

          if((bct == pinf) or (bct == pout) or (bct == minf)){

            const amrex::Real ep_s_extrap = 2.0*ep_s_arr(dom_lo[0],j,k) - ep_s_arr(dom_lo[0]+1,j,k);

            ep_s_arr(i,j,k) = amrex::max(0.0,
                              amrex::min(ep_s_extrap,
                              amrex::max(ep_cp,ep_s_arr(dom_lo[0],j,k))));

            Ps_arr(i,j,k) = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));

          }
        });
      } // nlft


      if (nrgt > 0) {

        amrex::ParallelFor(bx_yz_hi_2D,
        [bct_ihi,dom_hi,minf,pinf,pout,Ps_arr,ep_s_arr,
         Ps0,beta,ep_cp,small_number]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

          if((bct == pinf) or (bct == pout) or (bct == minf)){
            const amrex::Real ep_s_extrap = 2.0*ep_s_arr(dom_hi[0],j,k) - ep_s_arr(dom_hi[0]-1,j,k);

            ep_s_arr(i,j,k) = amrex::max(0.0,
                              amrex::min(ep_s_extrap,
                              amrex::max(ep_cp,ep_s_arr(dom_hi[0],j,k))));

            Ps_arr(i,j,k) = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));
          }
        });
      } // nrgt


      if (nbot > 0) {

        amrex::ParallelFor(bx_xz_lo_2D,
        [bct_jlo,dom_lo,minf,pinf,pout,Ps_arr,ep_s_arr,
         Ps0,beta,ep_cp, small_number]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

          if((bct == pinf) or (bct == pout) or (bct == minf)){

            const amrex::Real ep_s_extrap = 2.0*ep_s_arr(i,dom_lo[1],k) - ep_s_arr(i,dom_lo[1]+1,k);

            ep_s_arr(i,j,k) = amrex::max(0.0,
                              amrex::min(ep_s_extrap,
                              amrex::max(ep_cp,ep_s_arr(i,dom_lo[1],k))));

            Ps_arr(i,j,k) = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));

          }
        });
      } // nbot


      if (ntop > 0) {

        amrex::ParallelFor(bx_xz_hi_2D,
        [bct_jhi,dom_hi,minf,pinf,pout,Ps_arr,ep_s_arr,
         Ps0,beta,ep_cp,small_number]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

          if((bct == pinf) or (bct == pout) or (bct == minf)){
            const amrex::Real ep_s_extrap = 2.0*ep_s_arr(i,dom_hi[1],k) - ep_s_arr(i,dom_hi[1]-1,k);

            ep_s_arr(i,j,k) = amrex::max(0.0,
                              amrex::min(ep_s_extrap,
                              amrex::max(ep_cp,ep_s_arr(i,dom_hi[1],k))));

            Ps_arr(i,j,k) = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));
          }
        });
      } // ntop


      if (ndwn > 0) {

        amrex::ParallelFor(bx_xy_lo_2D,
        [bct_klo,dom_lo,minf,pinf,pout,Ps_arr,ep_s_arr,
         Ps0,beta,ep_cp, small_number]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_klo(i,j,dom_lo[2]-1,0);

          if((bct == pinf) or (bct == pout) or (bct == minf)){

            const amrex::Real ep_s_extrap = 2.0*ep_s_arr(i,j,dom_lo[2]) - ep_s_arr(i,j,dom_lo[2]+1);

            ep_s_arr(i,j,k) = amrex::max(0.0,
                              amrex::min(ep_s_extrap,
                              amrex::max(ep_cp,ep_s_arr(i,j,dom_lo[2]))));

            Ps_arr(i,j,k) = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));

          }
        });
      } // ndwn


      if (nup  > 0) {

        amrex::ParallelFor(bx_xy_hi_2D,
        [bct_khi,dom_hi,minf,pinf,pout,Ps_arr,ep_s_arr,
         Ps0,beta,ep_cp,small_number]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_khi(i,j,dom_hi[2]+1,0);

          if((bct == pinf) or (bct == pout) or (bct == minf)){
            const amrex::Real ep_s_extrap = 2.0*ep_s_arr(i,j,dom_hi[2]) - ep_s_arr(i,j,dom_hi[2]-1);

            ep_s_arr(i,j,k) = amrex::max(0.0,
                              amrex::min(ep_s_extrap,
                              amrex::max(ep_cp,ep_s_arr(i,j,dom_hi[2]))));

            Ps_arr(i,j,k) = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_arr(i,j,k));
          }
        });
      } // nup

    }

    Ps.FillBoundary(geom[lev].periodicity());



    // We now have the solids stress on the field. We now need to interpolate
    // this to the parcel positions. We are going to store this for the
    // parcel in "omega" in the particle container.

    using MFIXParIter = MFIXParticleContainer::MFIXParIter;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {

      const auto dxi_array = geom[lev].InvCellSizeArray();
      const auto plo_array = geom[lev].ProbLoArray();

      const amrex::RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
      const amrex::RealVect plo(plo_array[0], plo_array[1], plo_array[2]);


      for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
      {
        // Timer used for load-balancing
        amrex::Real wt = ParallelDescriptor::second();

        auto& particles = pti.GetArrayOfStructs();
        MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        auto& soa = pti.GetStructOfArrays();
        auto p_realarray = soa.realarray();

        const int np = particles.size();

        Box bx = pti.tilebox ();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  ep_s_fab = static_cast<EBFArrayBox const&>((*ep_s_in[lev])[pti]);
        const EBCellFlagFab&  flags = ep_s_fab.getEBCellFlagFab();

        Array4<const Real> const& ep_s_arr = ep_s_in[lev]->const_array(pti);
        Array4<const Real> const& Ps_arr   = Ps.array(pti);
        const auto& flags_array = flags.array();

        if (flags.getType(amrex::grow(bx,1)) == FabType::covered)
        {
          // We shouldn't have this case but if we do -- zero the stress.
          amrex::ParallelFor(np, [p_realarray]
            AMREX_GPU_DEVICE (int ip) noexcept
            {
              p_realarray[SoArealData::oneOverI][ip] = 0.0;

              p_realarray[SoArealData::omegax][ip] = 0.0;
              p_realarray[SoArealData::omegay][ip] = 0.0;
              p_realarray[SoArealData::omegaz][ip] = 0.0;
            });
        }
        else if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
        {
          amrex::ParallelFor(np, [pstruct,p_realarray,ep_s_arr,Ps_arr,plo,dxi]
            AMREX_GPU_DEVICE (int ip) noexcept
          {

            MFIXParticleContainer::ParticleType& particle = pstruct[ip];

            const amrex::Real lx = (particle.pos(0) - plo[0]) * dxi[0] + 0.5;
            const amrex::Real ly = (particle.pos(1) - plo[1]) * dxi[1] + 0.5;
            const amrex::Real lz = (particle.pos(2) - plo[2]) * dxi[2] + 0.5;

            const int i = static_cast<int>(amrex::Math::floor(lx));
            const int j = static_cast<int>(amrex::Math::floor(ly));
            const int k = static_cast<int>(amrex::Math::floor(lz));

            const amrex::Real wx_hi(lx - static_cast<amrex::Real>(i));
            const amrex::Real wy_hi(ly - static_cast<amrex::Real>(j));
            const amrex::Real wz_hi(lz - static_cast<amrex::Real>(k));

            const amrex::Real wx_lo(1.0 - wx_hi);
            const amrex::Real wy_lo(1.0 - wy_hi);
            const amrex::Real wz_lo(1.0 - wz_hi);

            AMREX_ASSERT(ep_s_arr(i-1,j-1, k-1) >= 0.);
            AMREX_ASSERT(ep_s_arr(i-1,j-1, k  ) >= 0.);
            AMREX_ASSERT(ep_s_arr(i-1,j  , k-1) >= 0.);
            AMREX_ASSERT(ep_s_arr(i-1,j  , k  ) >= 0.);
            AMREX_ASSERT(ep_s_arr(i  ,j-1, k-1) >= 0.);
            AMREX_ASSERT(ep_s_arr(i  ,j-1, k  ) >= 0.);
            AMREX_ASSERT(ep_s_arr(i  ,j  , k-1) >= 0.);
            AMREX_ASSERT(ep_s_arr(i  ,j  , k  ) >= 0.);

            AMREX_ASSERT(Ps_arr(i-1,j-1, k-1) >= 0.);
            AMREX_ASSERT(Ps_arr(i-1,j-1, k  ) >= 0.);
            AMREX_ASSERT(Ps_arr(i-1,j  , k-1) >= 0.);
            AMREX_ASSERT(Ps_arr(i-1,j  , k  ) >= 0.);
            AMREX_ASSERT(Ps_arr(i  ,j-1, k-1) >= 0.);
            AMREX_ASSERT(Ps_arr(i  ,j-1, k  ) >= 0.);
            AMREX_ASSERT(Ps_arr(i  ,j  , k-1) >= 0.);
            AMREX_ASSERT(Ps_arr(i  ,j  , k  ) >= 0.);

            const Real ep_s_loc =
              + wx_lo*wy_lo*wz_lo*ep_s_arr(i-1,j-1, k-1)
              + wx_lo*wy_lo*wz_hi*ep_s_arr(i-1,j-1, k  )
              + wx_lo*wy_hi*wz_lo*ep_s_arr(i-1,j  , k-1)
              + wx_lo*wy_hi*wz_hi*ep_s_arr(i-1,j  , k  )
              + wx_hi*wy_lo*wz_lo*ep_s_arr(i  ,j-1, k-1)
              + wx_hi*wy_lo*wz_hi*ep_s_arr(i  ,j-1, k  )
              + wx_hi*wy_hi*wz_lo*ep_s_arr(i  ,j  , k-1)
              + wx_hi*wy_hi*wz_hi*ep_s_arr(i  ,j  , k  );

            const Real dPsdx = dxi[0]*(
              - wy_lo*wz_lo*Ps_arr(i-1,j-1, k-1)
              - wy_lo*wz_hi*Ps_arr(i-1,j-1, k  )
              - wy_hi*wz_lo*Ps_arr(i-1,j  , k-1)
              - wy_hi*wz_hi*Ps_arr(i-1,j  , k  )
              + wy_lo*wz_lo*Ps_arr(i  ,j-1, k-1)
              + wy_lo*wz_hi*Ps_arr(i  ,j-1, k  )
              + wy_hi*wz_lo*Ps_arr(i  ,j  , k-1)
              + wy_hi*wz_hi*Ps_arr(i  ,j  , k  ));

            const Real dPsdy = dxi[1]*(
              - wx_lo*wz_lo*Ps_arr(i-1,j-1, k-1)
              - wx_lo*wz_hi*Ps_arr(i-1,j-1, k  )
              + wx_lo*wz_lo*Ps_arr(i-1,j  , k-1)
              + wx_lo*wz_hi*Ps_arr(i-1,j  , k  )
              - wx_hi*wz_lo*Ps_arr(i  ,j-1, k-1)
              - wx_hi*wz_hi*Ps_arr(i  ,j-1, k  )
              + wx_hi*wz_lo*Ps_arr(i  ,j  , k-1)
              + wx_hi*wz_hi*Ps_arr(i  ,j  , k  ));

            const Real dPsdz = dxi[2]*(
              - wx_lo*wy_lo*Ps_arr(i-1,j-1, k-1)
              + wx_lo*wy_lo*Ps_arr(i-1,j-1, k  )
              - wx_lo*wy_hi*Ps_arr(i-1,j  , k-1)
              + wx_lo*wy_hi*Ps_arr(i-1,j  , k  )
              - wx_hi*wy_lo*Ps_arr(i  ,j-1, k-1)
              + wx_hi*wy_lo*Ps_arr(i  ,j-1, k  )
              - wx_hi*wy_hi*Ps_arr(i  ,j  , k-1)
              + wx_hi*wy_hi*Ps_arr(i  ,j  , k  ));


            // Store local solids volume fraction in unused array location.
            p_realarray[SoArealData::oneOverI][ip] = ep_s_loc;

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_loc >= 0.,"Check valid ep_s_loc");

            p_realarray[SoArealData::omegax][ip] = dPsdx;
            p_realarray[SoArealData::omegay][ip] = dPsdy;
            p_realarray[SoArealData::omegaz][ip] = dPsdz;
          });
        }
        else // FAB not all regular
        {

          amrex::ParallelFor(np,
          [pstruct,p_realarray,flags_array,ep_s_arr,Ps_arr,plo,dxi,
           Ps0,beta,ep_cp,small_number]
          AMREX_GPU_DEVICE (int pid) noexcept
          {

            MFIXParticleContainer::ParticleType& particle = pstruct[pid];

            const int ip = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
            const int jp = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
            const int kp = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

            if(flags_array(ip,jp,kp).isCovered())
            {
              p_realarray[SoArealData::omegax][ip] = 0.0;
              p_realarray[SoArealData::omegay][ip] = 0.0;
              p_realarray[SoArealData::omegaz][ip] = 0.0;

            // Cut or regular cell and none of the cells in the stencil is covered
            // (Note we can't assume regular cell has no covered cells in the stencil
            //      because of the diagonal case)
            } else {

              const amrex::Real lx = (particle.pos(0) - plo[0]) * dxi[0] + 0.5;
              const amrex::Real ly = (particle.pos(1) - plo[1]) * dxi[1] + 0.5;
              const amrex::Real lz = (particle.pos(2) - plo[2]) * dxi[2] + 0.5;

              // Upper cell in trilinear stencil
              const int i = static_cast<int>(amrex::Math::floor(lx));
              const int j = static_cast<int>(amrex::Math::floor(ly));
              const int k = static_cast<int>(amrex::Math::floor(lz));

              const int di = i - ip; // 0 or 1
              const int dj = j - jp; // 0 or 1
              const int dk = k - kp; // 0 or 1

              amrex::Real ep_s_loc;

              amrex::Real dPsdx;
              amrex::Real dPsdy;
              amrex::Real dPsdz;

              // Check that all the cells in the stencil are connected to the
              // that contains the parcels centroid. We don't need to check
              // cell (ip,jp,kp) as we've already tested if it is covered.

              if (flags_array(ip,jp,kp).isConnected(di-1,dj-1,dk-1) and
                  flags_array(ip,jp,kp).isConnected(di  ,dj-1,dk-1) and
                  flags_array(ip,jp,kp).isConnected(di-1,dj  ,dk-1) and
                  flags_array(ip,jp,kp).isConnected(di  ,dj  ,dk-1) and
                  flags_array(ip,jp,kp).isConnected(di-1,dj-1,dk  ) and
                  flags_array(ip,jp,kp).isConnected(di  ,dj-1,dk  ) and
                  flags_array(ip,jp,kp).isConnected(di-1,dj  ,dk  ) and
                  flags_array(ip,jp,kp).isConnected(di  ,dj  ,dk  )) {

                // After sufficient testing, these should be changed to
                // standard AMREX_ASSERT checks
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_arr(i-1,j-1, k-1) >= 0.,"Check valid ep_s A");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_arr(i-1,j-1, k  ) >= 0.,"Check valid ep_s B");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_arr(i-1,j  , k-1) >= 0.,"Check valid ep_s C");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_arr(i-1,j  , k  ) >= 0.,"Check valid ep_s D");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_arr(i  ,j-1, k-1) >= 0.,"Check valid ep_s E");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_arr(i  ,j-1, k  ) >= 0.,"Check valid ep_s F");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_arr(i  ,j  , k-1) >= 0.,"Check valid ep_s G");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_arr(i  ,j  , k  ) >= 0.,"Check valid ep_s H");

                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Ps_arr(i-1,j-1, k-1) >= 0.,"Check valid Ps A");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Ps_arr(i-1,j-1, k  ) >= 0.,"Check valid Ps B");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Ps_arr(i-1,j  , k-1) >= 0.,"Check valid Ps C");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Ps_arr(i-1,j  , k  ) >= 0.,"Check valid Ps D");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Ps_arr(i  ,j-1, k-1) >= 0.,"Check valid Ps E");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Ps_arr(i  ,j-1, k  ) >= 0.,"Check valid Ps F");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Ps_arr(i  ,j  , k-1) >= 0.,"Check valid Ps G");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Ps_arr(i  ,j  , k  ) >= 0.,"Check valid Ps H");

                const amrex::Real wx_hi(lx - static_cast<amrex::Real>(i));
                const amrex::Real wy_hi(ly - static_cast<amrex::Real>(j));
                const amrex::Real wz_hi(lz - static_cast<amrex::Real>(k));

                const amrex::Real wx_lo(1.0 - wx_hi);
                const amrex::Real wy_lo(1.0 - wy_hi);
                const amrex::Real wz_lo(1.0 - wz_hi);

                ep_s_loc =
                  + wx_lo*wy_lo*wz_lo*ep_s_arr(i-1,j-1, k-1)
                  + wx_lo*wy_lo*wz_hi*ep_s_arr(i-1,j-1, k  )
                  + wx_lo*wy_hi*wz_lo*ep_s_arr(i-1,j  , k-1)
                  + wx_lo*wy_hi*wz_hi*ep_s_arr(i-1,j  , k  )
                  + wx_hi*wy_lo*wz_lo*ep_s_arr(i  ,j-1, k-1)
                  + wx_hi*wy_lo*wz_hi*ep_s_arr(i  ,j-1, k  )
                  + wx_hi*wy_hi*wz_lo*ep_s_arr(i  ,j  , k-1)
                  + wx_hi*wy_hi*wz_hi*ep_s_arr(i  ,j  , k  );

                dPsdx = dxi[0]*(
                  - wy_lo*wz_lo*Ps_arr(i-1,j-1, k-1)
                  - wy_lo*wz_hi*Ps_arr(i-1,j-1, k  )
                  - wy_hi*wz_lo*Ps_arr(i-1,j  , k-1)
                  - wy_hi*wz_hi*Ps_arr(i-1,j  , k  )
                  + wy_lo*wz_lo*Ps_arr(i  ,j-1, k-1)
                  + wy_lo*wz_hi*Ps_arr(i  ,j-1, k  )
                  + wy_hi*wz_lo*Ps_arr(i  ,j  , k-1)
                  + wy_hi*wz_hi*Ps_arr(i  ,j  , k  ));

                dPsdy = dxi[1]*(
                  - wx_lo*wz_lo*Ps_arr(i-1,j-1, k-1)
                  - wx_lo*wz_hi*Ps_arr(i-1,j-1, k  )
                  + wx_lo*wz_lo*Ps_arr(i-1,j  , k-1)
                  + wx_lo*wz_hi*Ps_arr(i-1,j  , k  )
                  - wx_hi*wz_lo*Ps_arr(i  ,j-1, k-1)
                  - wx_hi*wz_hi*Ps_arr(i  ,j-1, k  )
                  + wx_hi*wz_lo*Ps_arr(i  ,j  , k-1)
                  + wx_hi*wz_hi*Ps_arr(i  ,j  , k  ));

                dPsdz = dxi[2]*(
                  - wx_lo*wy_lo*Ps_arr(i-1,j-1, k-1)
                  + wx_lo*wy_lo*Ps_arr(i-1,j-1, k  )
                  - wx_lo*wy_hi*Ps_arr(i-1,j  , k-1)
                  + wx_lo*wy_hi*Ps_arr(i-1,j  , k  )
                  - wx_hi*wy_lo*Ps_arr(i  ,j-1, k-1)
                  + wx_hi*wy_lo*Ps_arr(i  ,j-1, k  )
                  - wx_hi*wy_hi*Ps_arr(i  ,j  , k-1)
                  + wx_hi*wy_hi*Ps_arr(i  ,j  , k  ));

              // At least one of the cells in the stencil is not connected
              } else {

                ep_s_loc = ep_s_arr(ip,jp,kp);

                // 1D approx of dPsdx

                if (flags_array(ip,jp,kp).isConnected(di-1,0,0) and
                    flags_array(ip,jp,kp).isConnected(di  ,0,0)){

                  dPsdx = dxi[0]*(Ps_arr(i,jp,kp) - Ps_arr(i-1,jp,kp));

                } else {

                  if(di == 0){

                    // Don't use the low side Ps because i==ip and ip is not covered.

                    const amrex::Real ep_s_extrap = amrex::max(0.0,
                      amrex::min(2.0*ep_s_arr(i,jp,kp) - ep_s_arr(i+1,jp,kp),
                      amrex::max(ep_cp, ep_s_arr(ip,jp,kp))));

                    const amrex::Real Ps_lo = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_extrap);

                    dPsdx = dxi[0]*(Ps_arr(ip,jp,kp) - Ps_lo);

                  } else {

                    // Don't use the high side Ps because i==ip+1 and ip is not covered.

                    const amrex::Real ep_s_extrap = amrex::max(0.0,
                      amrex::min(2.0*ep_s_arr(ip,jp,kp) - ep_s_arr(ip-1,jp,kp),
                      amrex::max(ep_cp, ep_s_arr(ip,jp,kp))));

                    const amrex::Real Ps_hi = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_extrap);

                    dPsdx = dxi[0]*(Ps_hi - Ps_arr(ip,jp,kp));

                  }

                }

                // 1D approx of dPsdy

                if (flags_array(ip,jp,kp).isConnected(0,dj-1,0) and
                    flags_array(ip,jp,kp).isConnected(0,dj  ,0)){

                  dPsdy = dxi[1]*(Ps_arr(ip,j,kp) - Ps_arr(ip,j-1,kp));

                } else {

                  if(dj == 0){ // Low side is NOT connected

                    const amrex::Real ep_s_extrap = amrex::max(0.0,
                      amrex::min(2.0*ep_s_arr(ip,j,kp) - ep_s_arr(ip,j+1,kp),
                      amrex::max(ep_cp, ep_s_arr(ip,jp,kp))));

                    const amrex::Real Ps_lo = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_extrap);

                    dPsdy = dxi[1]*(Ps_arr(ip,jp,kp) - Ps_lo);

                  } else {

                    // Don't use the high side Ps because i==ip+1 and ip is not covered.

                    const amrex::Real ep_s_extrap = amrex::max(0.0,
                      amrex::min(2.0*ep_s_arr(ip,jp,kp) - ep_s_arr(ip,jp-1,kp),
                      amrex::max(ep_cp, ep_s_arr(ip,jp,kp))));

                    const amrex::Real Ps_hi = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_extrap);

                    dPsdy = dxi[1]*(Ps_hi - Ps_arr(ip,jp,kp));

                  }

                }

                // 1D approx of dPsdz

                if (flags_array(ip,jp,kp).isConnected(0,0,dk-1) and
                    flags_array(ip,jp,kp).isConnected(0,0,dk  )){

                  dPsdz = dxi[2]*(Ps_arr(ip,jp,k) - Ps_arr(ip,jp,k-1));

                } else {

                  if(dk == 0){ // Low side is NOT connected

                    const amrex::Real ep_s_extrap = amrex::max(0.0,
                      amrex::min(2.0*ep_s_arr(ip,jp,k) - ep_s_arr(ip,jp,k+1),
                      amrex::max(ep_cp, ep_s_arr(ip,jp,kp))));

                    const amrex::Real Ps_lo = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_extrap);

                    dPsdz = dxi[2]*(Ps_arr(ip,jp,kp) - Ps_lo);

                  } else { // High side is NOT connected

                    const amrex::Real ep_s_extrap = amrex::max(0.0,
                      amrex::min(2.0*ep_s_arr(ip,jp,kp) - ep_s_arr(ip,jp,kp-1),
                      amrex::max(ep_cp, ep_s_arr(ip,jp,kp))));

                    const amrex::Real Ps_hi = solids_pressure(Ps0, beta, ep_cp, small_number, ep_s_extrap);

                    dPsdz = dxi[2]*(Ps_hi - Ps_arr(ip,jp,kp));

                  }

                }

              } // Cut cell

              // Store local solids volume fraction in unused array location.
              p_realarray[SoArealData::oneOverI][ip] = ep_s_loc;

              AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ep_s_loc >= 0.,"Check valid ep_s_loc");

              p_realarray[SoArealData::omegax][ip] = dPsdx;
              p_realarray[SoArealData::omegay][ip] = dPsdy;
              p_realarray[SoArealData::omegaz][ip] = dPsdz;
            }

          }); // part loop
        } // if box not all regular


        /********************************************************************
         * Update runtime cost (used in load-balancing)                     *
         *******************************************************************/
        if (cost[lev])
        {
          // Runtime cost is either (weighted by tile box size):
          //   * time spent
          //   * number of particles
          const Box& tbx = pti.tilebox();
          if (knapsack_weight_type == "RunTimeCosts")
          {
            wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
          }
          else if (knapsack_weight_type == "NumParticles")
          {
            wt = np / tbx.d_numPts();
          }
          (*cost[lev])[pti].plus<RunOn::Device>(wt, tbx);
        }

      } // pti
    }  // omp region
  }// lev

}
