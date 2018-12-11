#include <algorithm>
#include <mfix.H>
#include <mfix_eb_F.H>

void mfix::WriteEBSurface()
{
  if (Geom(0).isAllPeriodic()) return;

  // Only write at the finest level!
  int lev = nlev-1;

  const Real* dx = Geom(lev).CellSize();

  BoxArray & ba            = grids[lev];
  DistributionMapping & dm = dmap[lev];

  // This creates the associated Distribution Mapping
  // DistributionMapping dm(ba, ParallelDescriptor::NProcs());

  const EBFArrayBoxFactory * ebf;

  if (particle_ebfactory[lev] != nullptr) {
      ebf = particle_ebfactory[lev].get();
      ba  = pc->ParticleBoxArray(lev);
      dm  = pc->ParticleDistributionMap(lev);
  } else {
      ebf = ebfactory[lev].get();
  }

  MultiFab mf_ba(ba, dm, 1, 0, MFInfo(), * ebf);


  // // // Deliberately didn't time this loop.
  for (MFIter mfi(mf_ba); mfi.isValid(); ++mfi) {

    const auto& sfab = static_cast<EBFArrayBox const&>((mf_ba)[mfi]);
    const auto& my_flag = sfab.getEBCellFlagFab();

    const Box& bx = mfi.validbox();

    if (my_flag.getType(bx) == FabType::covered or my_flag.getType(bx) == FabType::regular) continue;

    std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
    const MultiCutFab * bndrycent;

    areafrac  =  ebf->getAreaFrac();
    bndrycent = &(ebf->getBndryCent());

    mfix_eb_to_polygon(dx, BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_3D(my_flag),
                       BL_TO_FORTRAN_3D((* bndrycent)[mfi]),
                       BL_TO_FORTRAN_3D((*areafrac[0])[mfi]),
                       BL_TO_FORTRAN_3D((*areafrac[1])[mfi]),
                       BL_TO_FORTRAN_3D((*areafrac[2])[mfi])   );
  }

  int cpu = ParallelDescriptor::MyProc();
  int nProcs = ParallelDescriptor::NProcs();

  mfix_write_eb_vtp(&cpu);

  if(ParallelDescriptor::IOProcessor())
    mfix_write_pvtp(&nProcs);


  // // // Deliberately didn't time this loop.
  for (MFIter mfi(mf_ba); mfi.isValid(); ++mfi) {

    const auto& sfab = static_cast<EBFArrayBox const&>((mf_ba)[mfi]);
    const auto& my_flag = sfab.getEBCellFlagFab();

    const Box& bx = mfi.validbox();

    if (my_flag.getType(bx) == FabType::covered or my_flag.getType(bx) == FabType::regular) continue;

    mfix_eb_grid_coverage(&cpu, dx, bx.loVect(), bx.hiVect(),
         my_flag.dataPtr(), my_flag.loVect(), my_flag.hiVect());
  }
}
