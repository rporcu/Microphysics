/********************************************************************************
 *                                                                              *
 * NOTE: after no longer using GeometryShop, this file will be deprecated       *
 *       -> delete after it is no longer needed as reference                    *
 *                                                                              *
 ********************************************************************************/



#include <algorithm>
#include <AMReX_GeometryShop.H>
#include <mfix_level.H>
#include <mfix_eb_F.H>

void
mfix_level::fix_up_geometry(EBFArrayBoxFactory* my_eb_factory, int lev)
{
    BL_PROFILE("mfix_level::fix_up_geometry()");

    const auto& domain = geom[lev].Domain();

    const amrex::MultiFab* volfrac   = &(my_eb_factory->getVolFrac());
    auto& cell_flag =  const_cast<FabArray<EBCellFlagFab>&> (my_eb_factory->getMultiEBCellFlagFab());

    int ng = volfrac->nGrow();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*volfrac); mfi.isValid(); ++mfi) 
    {
      const Box& bx = mfi.growntilebox(ng);

      auto& flag_fab = cell_flag[mfi];

      if (flag_fab.getType(bx) == FabType::singlevalued)
      {
          mfix_fixup_eb_geom(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_ANYD(flag_fab),
                             BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                             BL_TO_FORTRAN_BOX(domain));
      }
    }
}
