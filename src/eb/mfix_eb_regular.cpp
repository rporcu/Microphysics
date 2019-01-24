#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <mfix_eb_if.H>

#include <AMReX_EB_utils.H>
#include <AMReX_EB_levelset.H>


#include <mfix.H>
#include <mfix_eb_F.H>


/********************************************************************************
 *                                                                              *
 * Placeholder: create a simulation box _without_ EB walls.                     *
 *                                                                              *
 *******************************************************************************/
void
mfix::make_eb_regular()
{
    int max_level_here = 0;

    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ***************************************************************************/

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

    // if (solve_fluid)
    {
        bool has_walls = false;
        std::unique_ptr<UnionListIF<EB2::PlaneIF>> impfunc_walls = get_real_walls(has_walls);

        if (has_walls)
        {
            auto gshop = EB2::makeShop(* impfunc_walls);

            build_eb_levels(gshop);

            //___________________________________________________________________
            // Particles need the correct volfrac at the inflow
            bool has_more_walls = false;
            std::unique_ptr<UnionListIF<EB2::PlaneIF>> walls = get_walls(has_more_walls);
            if (has_more_walls)
            {
                auto if_part = EB2::makeUnion(* impfunc_walls, * walls);
                auto gshop_part = EB2::makeShop(if_part);

                build_particle_eb_levels(gshop_part);
            }

        }
        else
        {
            EB2::AllRegularIF my_regular;
            auto gshop = EB2::makeShop(my_regular);

            build_eb_levels(gshop);

            //___________________________________________________________________
            // Particles need the correct volfrac at the inflow
            bool has_more_walls = false;
            std::unique_ptr<UnionListIF<EB2::PlaneIF>> walls = get_walls(has_more_walls);
            if (has_more_walls)
            {
                // auto if_part = EB2::makeUnion(my_regular, * walls);
                auto gshop_part = EB2::makeShop(* walls);

                build_particle_eb_levels(gshop_part);
            }

        }
    }
}
