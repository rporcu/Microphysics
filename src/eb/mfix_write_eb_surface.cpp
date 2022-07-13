#include <mfix_rw.H>
#include <mfix_dem.H>

#include <AMReX_WriteEBSurface.H>

using namespace amrex;


namespace MfixIO {

void
MfixRW::WriteMyEBSurface () const
{
  if (geom[0].isAllPeriodic()) return;

  // Only write at the finest level!
  int lev = nlev-1;

  BoxArray & ba            = grids[lev];
  DistributionMapping & dm = dmap[lev];

  const EBFArrayBoxFactory * ebf;

  if (particle_ebfactory[lev] != nullptr && m_dem.solve()) {
      ebf = particle_ebfactory[lev].get();
      ba  = pc->ParticleBoxArray(lev);
      dm  = pc->ParticleDistributionMap(lev);
  } else {
      ebf = ebfactory[lev].get();
  }

  WriteEBSurface(ba,dm,geom[lev],ebf);
}

}
