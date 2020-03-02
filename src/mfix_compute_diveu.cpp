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
  Vector< MultiFab* > epu(nlev, nullptr);

  for (int lev = 0; lev < nlev; lev++)
    {
      MultiFab& vel_g = *(m_leveldata[lev]->vel_g);
      // We only need one ghost cell here -- so no need to make it bigger
      epu[lev] = new MultiFab(vel_g.boxArray(), vel_g.DistributionMap(),
                              vel_g.nComp(), 1 , MFInfo(), *ebfactory[lev]);

      epu[lev]->setVal(1.e200);

      Box domain(geom[lev].Domain());

      MultiFab::Copy(*epu[lev], vel_g, 0, 0, 3, epu[lev]->nGrow());

      for (int n = 0; n < 3; n++)
        MultiFab::Multiply(*epu[lev], *(m_leveldata[lev]->ep_g), 0, n, 1, epu[lev]->nGrow());

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
  LPInfo info;
  MLNodeLaplacian matrix(geom, grids, dmap, info, ebfactory);

  // Set domain BCs for Poisson's solver
  // The domain BCs refer to level 0 only
  Box domain(geom[0].Domain());

  matrix.setDomainBC(BC::ppe_lobc, BC::ppe_hibc);

  Vector< MultiFab* > diveu(m_leveldata.size(), nullptr);
  for (int lev(0); lev < m_leveldata.size(); ++lev)
    diveu[lev] = m_leveldata[lev]->diveu;

  matrix.compDivergence(diveu, epu);

  for(int lev(0); lev < nlev; lev++)
    delete epu[lev];

  Vector< MultiFab* > vel_g(m_leveldata.size(), nullptr);
  for(int lev(0); lev < m_leveldata.size(); ++lev)
    vel_g[lev] = m_leveldata[lev]->vel_g;

  // Restore velocities to carry Dirichlet values on faces
  int extrap_dir_bcs = 0;
  mfix_set_velocity_bcs(time, vel_g, extrap_dir_bcs);
}
