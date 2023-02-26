#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>

#include <mfix_bc.H>

using namespace amrex;

void MFIXBoundaryConditions::
calc_bc_areas (int const lev,
               BoxArray            const& a_grids,
               DistributionMapping const& a_dmap,
               EBFArrayBoxFactory*        a_factory)
{
  // Nothing to do for fully periodic systems!
  if ( m_geom[lev].isAllPeriodic() ) { return; }

  MultiFab tmpMF(a_grids, a_dmap, 1, 0, MFInfo(), *a_factory);

  int const num_bcs(m_bc.size());

  std::vector<Real> bc_areas(num_bcs, 0.);

  Box domain(m_geom[lev].Domain());

  if( !m_geom[lev].isPeriodic(0) ) { // begin x direction

    int const dir(0);

    { // x-lo side of the domain

      Box box_ilo = amrex::adjCellLo(domain,dir,1);

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_xlo().size(); ++lc) {

        const int bcv  = this->bc_xlo(lc);
        const BC_t& bc = this->bc(bcv);

        bc_areas[bcv] = calc_bc_area(lev, dir, box_ilo, bc, tmpMF, a_factory);
      } // loop over bc_xlo

    } // end x-lo side of the domain

    { // x-hi side of the domain
      Box box_ihi = amrex::adjCellHi(domain,dir,1);

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_xhi().size(); ++lc) {

        const int bcv  = this->bc_xhi(lc);
        const BC_t& bc = this->bc(bcv);

        bc_areas[bcv] = calc_bc_area(lev, dir, box_ihi, bc, tmpMF, a_factory);
      }
    } // end x-hi side of the domain
  } // end x-direction


  if( !m_geom[lev].isPeriodic(1) ) { // begin y-direction

    const int dir = 1;

    { // y-lo side of the domain

      Box box_jlo = amrex::adjCellLo(domain,dir,1);

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_ylo().size(); ++lc) {

        const int bcv  = this->bc_ylo(lc);
        const BC_t& bc = this->bc(bcv);

        bc_areas[bcv] = calc_bc_area(lev, dir, box_jlo, bc, tmpMF, a_factory);
      }

    }// end y-lo side of the domain

    { // y-hi side of the domain

      Box box_jhi = amrex::adjCellHi(domain,dir,1);

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_yhi().size(); ++lc) {

        const int bcv  = this->bc_yhi(lc);
        const BC_t& bc = this->bc(bcv);

        bc_areas[bcv] = calc_bc_area(lev, dir, box_jhi, bc, tmpMF, a_factory);
      }

    } // end y-hi side of the domain
  } // end y-direction

  { // begin z-direction

    const int dir(2);

    { // z-lo side of the domain

      Box box_klo = amrex::adjCellLo(domain,dir,1);

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_zlo().size(); ++lc) {

        const int bcv  = this->bc_zlo(lc);
        const BC_t& bc = this->bc(bcv);

        bc_areas[bcv] = calc_bc_area(lev, dir, box_klo, bc, tmpMF, a_factory);
      }


    } // end z-lo side of the domain


    { // z-hi side of the domain

      Box box_khi = amrex::adjCellHi(domain,dir,1);

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_zhi().size(); ++lc) {

        const int bcv  = this->bc_zhi(lc);
        const BC_t& bc = this->bc(bcv);

        bc_areas[bcv] = calc_bc_area(lev, dir, box_khi, bc, tmpMF, a_factory);
      }

    } // end z-hi side of the domain
  }// end z-direction

  // Do a global reduce sum and store in data container

  ParallelDescriptor::ReduceRealSum(bc_areas.data(), bc_areas.size());

  for(int bcv(0); bcv < num_bcs; bcv++) {
    this->m_area[bcv] = bc_areas[bcv];
  }

}


Real
MFIXBoundaryConditions::
calc_bc_area( int      const       a_lev,
              int      const       a_dir,
              Box      const&      a_bx,
              BC_t     const&      a_bc,
              MultiFab const&      a_MF,
              EBFArrayBoxFactory*  a_factory)
{

 Real dx = m_geom[a_lev].CellSize(0);
 Real dy = m_geom[a_lev].CellSize(1);
 Real dz = m_geom[a_lev].CellSize(2);

 // This should be caught elsewhere but just in case...
 AMREX_ASSERT(dx == dy && dy == dz);
 Real const da( dx*dx );

  const GpuArray<Real, 3> plo = m_geom[a_lev].ProbLoArray();

  IntVect bx_lo(a_bx.loVect());
  IntVect bx_hi(a_bx.hiVect());

  if (a_dir != 0){
    bx_lo[0] = static_cast<int>(amrex::Math::floor((a_bc.region->lo(0)-plo[0])/dx + 0.5));
    bx_hi[0] = static_cast<int>(amrex::Math::floor((a_bc.region->hi(0)-plo[0])/dx + 0.5)-1);
  }
  if (a_dir != 1){
    bx_lo[1] = static_cast<int>(amrex::Math::floor((a_bc.region->lo(1)-plo[1])/dy + 0.5));
    bx_hi[1] = static_cast<int>(amrex::Math::floor((a_bc.region->hi(1)-plo[1])/dy + 0.5)-1);
  }
  if (a_dir != 2){
    bx_lo[2] = static_cast<int>(amrex::Math::floor((a_bc.region->lo(2)-plo[2])/dz + 0.5));
    bx_hi[2] = static_cast<int>(amrex::Math::floor((a_bc.region->hi(2)-plo[2])/dz + 0.5)-1);
  }

  const Box bc_bx(bx_lo, bx_hi);

  Real regCellArea(0.);

  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  // No tiling
  for (MFIter mfi(a_MF, false); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.growntilebox(IntVect::TheDimensionVector(a_dir));

    // Ensure lo_box intersects bx
    if (bc_bx.intersects(bx)) {

      const Box bx_int = bx&(bc_bx);

      EBCellFlagFab const& flagfab = a_factory->getMultiEBCellFlagFab()[mfi];
      FabType t = flagfab.getType(bx);

      if (t == FabType::regular ) {

        regCellArea += static_cast<Real>(bx_int.numPts());

      } else if (t == FabType::singlevalued ) {

        Array4<Real const> const& afrac = (a_factory->getAreaFrac()[a_dir])->const_array(mfi);

        reduce_op.eval(bx_int, reduce_data, [afrac]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        { return afrac(i,j,k); });
      }
    } // a_bx and the MFIter box intersect
  } // MFIter loop

  ReduceTuple host_tuple = reduce_data.value(reduce_op);
  Real const cutCellArea = amrex::get<0>(host_tuple);

  Real totalArea = da*(regCellArea + cutCellArea);

  return totalArea;

}
