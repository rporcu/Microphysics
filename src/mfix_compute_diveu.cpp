#include <mfix_F.H>
#include <mfix.H>
#include <MFIX_BC_Parms.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLNodeLaplacian.H>

//
// Compute div(ep_g * u)
//
void
mfix::mfix_compute_diveu (Real time)
{
  // Note that the solver imposes the boundary conditions with the right scalings so we don't
  //      fill any ghost cells here.
  Vector<std::unique_ptr<MultiFab> > epu;
  epu.resize(nlev);

  for (int lev = 0; lev < nlev; lev++)
    {
      // We only need one ghost cell here -- so no need to make it bigger
      epu[lev].reset(new MultiFab(vel_g[lev]->boxArray(), vel_g[lev]->DistributionMap(),
                                  vel_g[lev]->nComp(), 1 , MFInfo(), *ebfactory[lev]));

      epu[lev]->setVal(1.e200);

      Box domain(geom[lev].Domain());

      MultiFab::Copy(*epu[lev], *vel_g[lev], 0, 0, 3, epu[lev]->nGrow());

      for (int n = 0; n < 3; n++)
        MultiFab::Multiply(*epu[lev], *ep_g[lev], 0, n, 1, epu[lev]->nGrow());


      epu[lev]->FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      // Extrapolate Dirichlet values to ghost cells -- but do it differently in that
      //  no-slip walls are treated exactly like slip walls --
      // Note that this routine is essential to impose the correct inflow bc's on
      //  the product ep_g * vel_g
      for (MFIter mfi(*epu[lev], false); mfi.isValid(); ++mfi)
      {
        set_vec_bcs(lev, (*epu[lev])[mfi], domain);
      }

      epu[lev]->FillBoundary(geom[lev].periodicity());

      // We set these to zero because if the values in the covered cells are undefined,
      //   even though they are multiplied by zero in the divu computation, we can still get NaNs
      EB_set_covered(*epu[lev], 0, epu[lev]->nComp(), 1, 0.0);
    }

  // Define the operator in order to compute the multi-level divergence
  //
  //        (del dot b sigma grad)) phi
  //
  LPInfo          info;
  MLNodeLaplacian matrix(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));

  // Set domain BCs for Poisson's solver
  // The domain BCs refer to level 0 only
  int bc_lo[3], bc_hi[3];
  Box domain(geom[0].Domain());

  matrix.setDomainBC(BC::ppe_lobc, BC::ppe_hibc);

  matrix.compDivergence(GetVecOfPtrs(diveu), GetVecOfPtrs(epu));

  // Restore velocities to carry Dirichlet values on faces
  int extrap_dir_bcs = 0;
  mfix_set_velocity_bcs(time, vel_g, extrap_dir_bcs);
}
