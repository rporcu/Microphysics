#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <mfix_eb_if.H>

#include <AMReX_EB_utils.H>

#include <mfix.H>


/********************************************************************************
 *                                                                              *
 * Placeholder: create a simulation box _without_ EB walls.                     *
 *                                                                              *
 *******************************************************************************/
void
mfix::make_eb_regular ()
{
    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ***************************************************************************/

    // set up ebfactory

    amrex::Print() << " " << std::endl;
    amrex::Print() << "Now making the ebfactories ..." << std::endl;

    // If filling level-set: this is used to store the implicit function (due to
    // any walls defined in mfix.dat). It is filled while after EB2::Build.
    // NOTE: this pointer is undefined if and of:
    //     * ! m_dem.solve()
    //     * levelset_restart
    //     * ! has_walls
    // are true
    
    //MultiFab* mf_impfunc;

    // if (fluid.solve())
    {
        bool has_walls = false;
        std::shared_ptr<UnionListIF<EB2::PlaneIF>> impfunc_walls = get_real_walls(has_walls);

        if (has_walls)
        {
            auto gshop = EB2::makeShop(*impfunc_walls);

            build_eb_levels(gshop);

            // Since walls are specified in the mfix.dat, setting the
            // `contains_ebs` flag to true
            contains_ebs = true;
        }
        else
        {
            EB2::AllRegularIF my_regular;
            auto gshop = EB2::makeShop(my_regular);

            build_eb_levels(gshop);
        }
    }
}
