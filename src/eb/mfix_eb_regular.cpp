#include <AMReX_EB2.H>

#include <AMReX_GeometryShop.H>
#include <AMReX_AllRegularService.H>
#include <AMReX_EBISLayout.H>
#include <AMReX_EBGraph.H>
#include <AMReX_EBDebugOut.H>
#include <AMReX_EBCellFAB.H>
#include <AMReX_EBCellFactory.H>
#include <AMReX_EBIndexSpace.H>

#include <algorithm>
#include <mfix_level.H>
#include <mfix_eb_F.H>


/****************************************************************************
 * Function to create a simple rectangular box with EB walls.               *
 *                                                                          *
 ****************************************************************************/
void
mfix_level::make_eb_regular(int lev)
{
    int max_level_here = 0;

    /***************************************************************************
     *                                                                         *
     * Build standard EB Factories                                             *
     *                                                                         *
     ***************************************************************************/

    // set up ebfactory
    int m_eb_basic_grow_cells = nghost;
    int m_eb_volume_grow_cells = nghost;
    int m_eb_full_grow_cells = nghost;
    EBSupport m_eb_support_level = EBSupport::full;

    EB2::useEB2(true);

    int max_coarsening_level = 100;
    Real offset    = 1.0e-8;

    amrex::Print() << " " << std::endl;
    amrex::Print() << "Now making the ebfactory's ..." << std::endl;

    EB2::AllRegularIF my_regular;
    auto gshop = EB2::makeShop(my_regular);
    EB2::Build(gshop, geom.back(), max_level_here,
               max_level_here+max_coarsening_level);

    const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
    const EB2::Level& eb_level = eb_is.getLevel(geom[lev]);

    if (solve_fluid)
       ebfactory[lev].reset(new EBFArrayBoxFactory(eb_level,
                geom[lev], grids[lev], dmap[lev],
                {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                 m_eb_full_grow_cells}, m_eb_support_level));

    if (solve_dem)
       particle_ebfactory.reset(new EBFArrayBoxFactory(eb_level,
                geom[lev], grids[lev], dmap[lev],
                {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                 m_eb_full_grow_cells}, m_eb_support_level));
}
