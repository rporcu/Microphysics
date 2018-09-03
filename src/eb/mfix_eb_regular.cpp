#include <AMReX_EB2.H>

#include <mfix_level.H>
#include <mfix_eb_F.H>


/********************************************************************************
 *                                                                              *
 * Placeholder: create a simulation box _without_ EB walls.                     *
 *                                                                              *
 ********************************************************************************/
void
mfix_level::make_eb_regular(int lev)
{
    int max_level_here = 0;

    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ****************************************************************************/

    // set up ebfactory
    int m_eb_basic_grow_cells = nghost;
    int m_eb_volume_grow_cells = nghost;
    int m_eb_full_grow_cells = nghost;
    EBSupport m_eb_support_level = EBSupport::full;

    int max_coarsening_level = 100;

    amrex::Print() << " " << std::endl;
    amrex::Print() << "Now making the ebfactory's ..." << std::endl;

    EB2::AllRegularIF my_regular;
    auto gshop = EB2::makeShop(my_regular);
    EB2::Build(gshop, geom.back(), max_level_here,
               max_level_here + max_coarsening_level);

    const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
    eb_level_fluid     = & eb_is.getLevel(geom[lev]);
    eb_level_particles =   eb_level_fluid;

    if (solve_fluid)
       ebfactory[lev].reset(new EBFArrayBoxFactory(
                                    * eb_level_fluid,
                                    geom[lev], grids[lev], dmap[lev],
                                    {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                     m_eb_full_grow_cells}, m_eb_support_level
                                )
                            );

    if (solve_dem) {
       particle_ebfactory[lev].reset(new EBFArrayBoxFactory(
                                            * eb_level_particles,
                                            geom[lev], grids[lev], dmap[lev],
                                            {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                             m_eb_full_grow_cells}, m_eb_support_level
                                        )
                                     );

       eb_normals = pc->EBNormals(lev, particle_ebfactory[lev].get(), dummy.get());

       /*************************************************************************
        *                                                                       *
        * Fill level-set:                                                       *
        * NOTE: this is necessary so that the ls_data MultiFab (as well as the  *
        *       level_set LSFactory is not full of junk. This will break if     *
        *       particle radius > 1                                             *
        *                                                                       *
        *************************************************************************/

       GShopLSFactory<EB2::AllRegularIF> reg_lsfactory(gshop, * level_set);
       std::unique_ptr<MultiFab> mf_impfunc = reg_lsfactory.fill_impfunc();

       level_set->intersection_impfunc( * mf_impfunc);

       // store copy of level set (for plotting).
       std::unique_ptr<MultiFab> ls_data = level_set->coarsen_data();
       ls[lev]->copy(* ls_data, 0, 0, 1, 0, 0);
       ls[lev]->FillBoundary(geom[lev].periodicity());
    }
}
