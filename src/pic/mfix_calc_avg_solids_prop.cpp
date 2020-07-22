#include <AMReX.H>
#include <mfix.H>
#include <mfix_mf_helpers.H>

using namespace amrex;

void mfix::MFIX_CalcAvgSolidsProp (amrex::Vector< amrex::MultiFab* >& avg_prop_in)
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
          amrex::Real mass = avg_prop_arr(i,j,k,3);
          const amrex::Real inv_mass = (mass > tolerance) ? 1.0 / mass : 0.0;
          avg_prop_arr(i,j,k,0) *= inv_mass;
          avg_prop_arr(i,j,k,1) *= inv_mass;
          avg_prop_arr(i,j,k,2) *= inv_mass;
        });
      }

      // A mix of regular and cut-cells.
      else if ( (*flags)[mfi].getType(amrex::grow(bx,0)) != FabType::covered)
      {

        const auto& flagsarr = (*flags)[mfi].array();

        amrex::ParallelFor(bx,
          [tolerance, avg_prop_arr, flagsarr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          amrex::Real mass = avg_prop_arr(i,j,k,3);
          if (not flagsarr(i,j,k).isCovered()) {
            const amrex::Real inv_mass = (mass > tolerance) ? 1.0 / mass : 0.0;
            avg_prop_arr(i,j,k,0) *= inv_mass;
            avg_prop_arr(i,j,k,1) *= inv_mass;
            avg_prop_arr(i,j,k,2) *= inv_mass;
          }
        });
      }
    }

    const Real avg_covered_val = 1.0e20;
    EB_set_covered(*avg_prop_in[lev], 0, 3, 1, avg_covered_val);

    // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
    avg_prop_in[lev]->FillBoundary(geom[lev].periodicity());

  }


  // Interpolate from the coarse grid to define the fine grid arrays

  if (nlev > 1) {

    amrex::Abort(" HACK: For right now the PIC implementation only works"
                 " with a single level. If we want to make sure it works"
                 " for more, then we need to fix this.");
  }

}
