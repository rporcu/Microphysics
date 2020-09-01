#include <AMReX.H>
#include <mfix.H>
#include <mfix_mf_helpers.H>

using namespace amrex;

void mfix::MFIX_CalcAvgSolidsVel (amrex::Vector< amrex::MultiFab* >& avg_prop_in)
{
  BL_PROFILE("MFIX_CalcAvgSolidsProp()");

  const amrex::FabArray<EBCellFlagFab>* flags = nullptr;

  for (int lev = 0; lev < nlev; lev ++ )
  {

    const amrex::BoxArray&            pba = pc->ParticleBoxArray(lev);
    const amrex::DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

    // Use level 0 to define the EB factory. If we are not on level 0
    // then create a copy of the coarse factory to use.
   if (lev == 0) {
      flags   = &(particle_ebfactory[lev]->getMultiEBCellFlagFab());

    } else {

      amrex::Vector<int> ngrow = {1,1,1};
      amrex::EBFArrayBoxFactory* crse_factory;

      crse_factory = (amrex::makeEBFabFactory(geom[0], pba, pdm, ngrow, EBSupport::volume)).release();

      flags   = &(crse_factory->getMultiEBCellFlagFab());

      delete crse_factory;
    }

    // We deposit vel*pmass and pmass. Later on we need to divide out
    // the accumulated mass to get the mass-averaged grid velocities.
    pc->MFIX_PC_SolidsVelocityDeposition(lev, *avg_prop_in[lev], flags);
  }


  {
    // The deposition occurred on level 0, thus the next few operations
    // only need to be carried out on level 0.
    int lev(0);

    // Sum grid boundaries to capture any material that was deposited into
    // your grid from an adjacent grid.
    avg_prop_in[lev]->SumBoundary(geom[lev].periodicity());

  // Compute the mass-averaged velocity.
    const Real tolerance = std::numeric_limits<Real>::epsilon();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*avg_prop_in[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      // Tilebox
      const Box& bx = mfi.tilebox();

      Array4<Real> const& avg_prop_arr = avg_prop_in[lev]->array(mfi);

      if ( (*flags)[mfi].getType(amrex::grow(bx,0)) == FabType::regular)
      {
        amrex::ParallelFor(bx, [tolerance, avg_prop_arr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

          if(avg_prop_arr(i,j,k,3) > tolerance)
            avg_prop_arr(i,j,k,0) *= (1.0/avg_prop_arr(i,j,k,3));
          else
            avg_prop_arr(i,j,k,0) = 0.0;

          if(avg_prop_arr(i,j,k,4) > tolerance)
            avg_prop_arr(i,j,k,1) *= (1.0/avg_prop_arr(i,j,k,4));
          else
            avg_prop_arr(i,j,k,1) = 0.0;

          if(avg_prop_arr(i,j,k,5) > tolerance)
            avg_prop_arr(i,j,k,2) *= (1.0/avg_prop_arr(i,j,k,5));
          else
            avg_prop_arr(i,j,k,2) = 0.0;


        });
        amrex::Gpu::synchronize();
      }

      // A mix of regular and cut-cells.
      else if ( (*flags)[mfi].getType(amrex::grow(bx,0)) != FabType::covered)
      {

        const auto& flagsarr = (*flags)[mfi].array();

        amrex::ParallelFor(bx,
          [tolerance, avg_prop_arr, flagsarr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (not flagsarr(i,j,k).isCovered()) {

            if(avg_prop_arr(i,j,k,3) > tolerance)
              avg_prop_arr(i,j,k,0) *= (1.0/avg_prop_arr(i,j,k,3));
            else
              avg_prop_arr(i,j,k,0) = 0.0;

            if(avg_prop_arr(i,j,k,4) > tolerance)
              avg_prop_arr(i,j,k,1) *= (1.0/avg_prop_arr(i,j,k,4));
            else
              avg_prop_arr(i,j,k,1) = 0.0;

            if(avg_prop_arr(i,j,k,5) > tolerance)
              avg_prop_arr(i,j,k,2) *= (1.0/avg_prop_arr(i,j,k,5));
            else
              avg_prop_arr(i,j,k,2) = 0.0;
          }
        });
        amrex::Gpu::synchronize();
      }

    }

    Box domain(geom[lev].Domain());

    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());

    // Set the face velocities for mass and pressure inflows to zero
    // so that if used in calculating the average solids velocity, the
    // result is zero (a non-moving wall).
    for (MFIter mfi(*avg_prop_in[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      Array4<Real> const& avg_prop_arr = avg_prop_in[lev]->array(mfi);

      IntVect arr_lo((*avg_prop_in[lev])[mfi].loVect());
      IntVect arr_hi((*avg_prop_in[lev])[mfi].hiVect());

      const Box arr_bx(arr_lo, arr_hi);

      Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
      Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
      Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
      Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
      Array4<const int> const& bct_klo = bc_klo[lev]->array();
      Array4<const int> const& bct_khi = bc_khi[lev]->array();

      const int nlft = amrex::max(0, dom_lo[0]-arr_lo[0]);
      const int nbot = amrex::max(0, dom_lo[1]-arr_lo[1]);
      const int ndwn = amrex::max(0, dom_lo[2]-arr_lo[2]);

      const int nrgt = amrex::max(0, arr_hi[0]-dom_hi[0]);
      const int ntop = amrex::max(0, arr_hi[1]-dom_hi[1]);
      const int nup  = amrex::max(0, arr_hi[2]-dom_hi[2]);

      // Create InVects for following 2D Boxes
      IntVect bx_yz_lo_lo_2D(arr_lo), bx_yz_lo_hi_2D(arr_hi);
      IntVect bx_yz_hi_lo_2D(arr_lo), bx_yz_hi_hi_2D(arr_hi);
      IntVect bx_xz_lo_lo_2D(arr_lo), bx_xz_lo_hi_2D(arr_hi);
      IntVect bx_xz_hi_lo_2D(arr_lo), bx_xz_hi_hi_2D(arr_hi);
      IntVect bx_xy_lo_lo_2D(arr_lo), bx_xy_lo_hi_2D(arr_hi);
      IntVect bx_xy_hi_lo_2D(arr_lo), bx_xy_hi_hi_2D(arr_hi);

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
      //const int pout = bc_list.get_pout(); UNUSED

      if (nlft > 0) {
        amrex::ParallelFor(bx_yz_lo_2D,
        [bct_ilo,dom_lo,minf,pinf,avg_prop_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

          if((bct == pinf) or (bct == minf))
            avg_prop_arr(dom_lo[0],j,k,0) = 0.0;

        });
        amrex::Gpu::synchronize();
      } // nlft

      if (nrgt > 0 ) {
        amrex::ParallelFor(bx_yz_hi_2D,
        [bct_ihi,dom_hi,minf,pinf,avg_prop_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

          if((bct == pinf) or (bct == minf))
            avg_prop_arr(dom_hi[0]+1,j,k,0) = 0.0;

        });
        amrex::Gpu::synchronize();
      } // nrgt

      if (nbot > 0) {
        amrex::ParallelFor(bx_xz_lo_2D,
        [bct_jlo,dom_lo,minf,pinf,avg_prop_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

          if((bct == pinf) or (bct == minf))
            avg_prop_arr(i,dom_lo[1],k,1) = 0.0;

        });
        amrex::Gpu::synchronize();
      } // nbot

      if (ntop > 0 ) {
        amrex::ParallelFor(bx_xz_hi_2D,
        [bct_jhi,dom_hi,minf,pinf,avg_prop_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

          if((bct == pinf) or (bct == minf))
            avg_prop_arr(i,dom_hi[1]+1,k,1) = 0.0;

        });
        amrex::Gpu::synchronize();
      } // ntop

      if (ndwn > 0) {
        amrex::ParallelFor(bx_xy_lo_2D,
        [bct_klo,dom_lo,minf,pinf,avg_prop_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_klo(i,j,dom_lo[2]-1,0);

          if((bct == pinf) or (bct == minf))
            avg_prop_arr(i,j,dom_lo[2],2) = 0.0;

        });
        amrex::Gpu::synchronize();
      } // ndwn

      if (nup  > 0 ) {
        amrex::ParallelFor(bx_xy_hi_2D,
        [bct_khi,dom_hi,minf,pinf,avg_prop_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

          const int bct = bct_khi(i,j,dom_hi[2]+1,0);

          if((bct == pinf) or (bct == minf))
            avg_prop_arr(i,j,dom_hi[2]+1,2) = 0.0;

        });
        amrex::Gpu::synchronize();
      } // nup

      amrex::Gpu::synchronize();
    }

    // Set covered cells and any remaining domain boundary cells to zero.
    // By setting these cells to zero, a parcel will not "see" any motion
    // in walls cells.

    const Real avg_covered_val = 0.0;
    EB_set_covered(*avg_prop_in[lev], 0, 3, 1, avg_covered_val);
    avg_prop_in[lev]->setDomainBndry(avg_covered_val,geom[lev]);

    // Impose periodic BC's at domain boundaries and fine-fine copies in the interior
    avg_prop_in[lev]->FillBoundary(geom[lev].periodicity());

  }

  // Interpolate from the coarse grid to define the fine grid arrays
  if (nlev > 1) {

    amrex::Abort(" HACK: For right now the PIC implementation only works"
                 " with a single level. If we want to make sure it works"
                 " for more, then we need to fix this.");
  }

}
