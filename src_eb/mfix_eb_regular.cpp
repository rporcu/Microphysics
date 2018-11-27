#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <mfix_eb_if.H>

#include <AMReX_EB_levelset.H>


#include <mfix.H>
#include <mfix_eb_F.H>


/********************************************************************************
 *                                                                              *
 * Placeholder: create a simulation box _without_ EB walls.                     *
 *                                                                              *
 ********************************************************************************/
void
mfix::make_eb_regular()
{
    int max_level_here = 0;

    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ****************************************************************************/

    // set up ebfactory
    EBSupport m_eb_support_level = EBSupport::full;

    int max_coarsening_level = 100;

    amrex::Print() << " " << std::endl;
    amrex::Print() << "Now making the ebfactories ..." << std::endl;

    // If filling level-set: this is used to store the implicit function (due to
    // any walls defined in mfix.dat). It is filled while after EB2::Build.
    // NOTE: this pointer is undefined if and of:
    //     * ! solve_dem
    //     * levelset__restart
    //     * ! has_walls
    // are true
    std::unique_ptr<MultiFab> mf_impfunc;

    if (solve_fluid) {
        bool has_walls = false;
        std::unique_ptr<UnionListIF<EB2::PlaneIF>> impfunc_walls = get_real_walls(has_walls);

        if (has_walls){
            auto gshop = EB2::makeShop(* impfunc_walls);
            EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
        } else {
            EB2::AllRegularIF my_regular;
            auto gshop = EB2::makeShop(my_regular);
            EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
        }

        const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();

        for (int lev = 0; lev < nlev; lev++)
        {
           eb_level_fluid = & eb_is.getLevel(geom[lev]);
           ebfactory[lev].reset(new EBFArrayBoxFactory(
                                        * eb_level_fluid,
                                        geom[lev], grids[lev], dmap[lev],
                                        {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                      m_eb_full_grow_cells}, m_eb_support_level));
        }
       }

    // Do _not_ fill level-set with AllRegularIF => if there are no walls, then
    // the level-set function is just huge(amrex_real) => this flag is set to
    // true iff there are walls.
    bool has_walls = false;

    if (solve_dem)
    {
        std::unique_ptr<UnionListIF<EB2::PlaneIF>> impfunc_walls = get_walls(has_walls);

        if (has_walls)
        {
            auto gshop = EB2::makeShop(* impfunc_walls);
            EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);

            if (! levelset__restart) {
                GShopLSFactory<UnionListIF<EB2::PlaneIF>> reg_lsfactory(gshop, * level_set);
                mf_impfunc = reg_lsfactory.fill_impfunc();
            }

        } else {
            EB2::AllRegularIF my_regular;
            auto gshop = EB2::makeShop(my_regular);
            EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
        }

        for (int lev = 0; lev < nlev; lev++)
        {
           const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
           eb_level_particles = & eb_is.getLevel(geom[lev]);

           particle_ebfactory[lev].reset(
               new EBFArrayBoxFactory(* eb_level_particles, geom[lev], grids[lev], dmap[lev],
                                      {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                                      m_eb_support_level)
               );

           eb_normals = pc->EBNormals(lev, particle_ebfactory[lev].get(), dummy.get());
        }

        /*************************************************************************
         *                                                                       *
         * Fill level-set:                                                       *
         * NOTE: this is necessary so that the ls_data MultiFab (as well as the  *
         *       level_set LSFactory is not full of junk. This will break if     *
         *       particle radius > 1                                             *
         *                                                                       *
         *************************************************************************/

        if (has_walls) {
            if (!levelset__restart) level_set->intersection_impfunc( * mf_impfunc);
            else
                amrex::Print() << "Loaded level-set is fine => skipping levelset calculation."
                               << std::endl;
        }

        // store copy of level set (for plotting).
        std::unique_ptr<MultiFab> ls_data = level_set->coarsen_data();
        int ng = ls_data->nGrow();

        for (int lev = 0; lev < nlev; lev++)
            ls[lev]->copy(* ls_data, 0, 0, 1, ng, ng);
    }
}
