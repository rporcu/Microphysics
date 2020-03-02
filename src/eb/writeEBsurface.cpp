#include <mfix.H>
#include <AMReX_WriteEBSurface.H>

void mfix::WriteMyEBSurface ()
{
  if (Geom(0).isAllPeriodic()) return;

  // Only write at the finest level!
  int lev = nlev-1;

  BoxArray & ba            = grids[lev];
  DistributionMapping & dm = dmap[lev];

  const EBFArrayBoxFactory * ebf;

  if (particle_ebfactory[lev] != nullptr) {
      ebf = particle_ebfactory[lev];
      ba  = pc->ParticleBoxArray(lev);
      dm  = pc->ParticleDistributionMap(lev);
  } else {
      ebf = ebfactory[lev];
  }

  WriteEBSurface(ba,dm,Geom(lev),ebf);
}
