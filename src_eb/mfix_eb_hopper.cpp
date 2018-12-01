#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Rotation.H>

#include <algorithm>
#include <AMReX_EB_utils.H>
#include <AMReX_EB_LSCore.H>
#include <AMReX_EB_levelset.H>
#include <mfix.H>
#include <mfix_eb_F.H>


/********************************************************************************
 *                                                                              *
 * Function to create a simple hopper EB.                                       *
 *                                                                              *
 * Comments: The hopper is constructed by "lathing" a union of two planes       *
 * (representing the hopper wall and outlet funnel) around a central axis.      *
 *                                                                              *
 *******************************************************************************/
void
mfix::make_eb_hopper()
{
    ParmParse pp("hopper");

    int max_level_here = 0;

    /****************************************************************************
     * Get hopper information from inputs file.                                 *
     ***************************************************************************/

    Real orifice_radius  = 0.00025;
    Real funnel_radius   = 0.00050;
    Real funnel_height   = 0.00105;

    int direction = 0;
    Vector<Real> centervec(3);

    pp.query("orifice_radius", orifice_radius);
    pp.query("funnel_radius",  funnel_radius);
    pp.query("funnel_height",  funnel_height);

    pp.query("direction", direction);
    pp.getarr("center",   centervec,  0, 3);
    Array<Real,3> center={centervec[0],centervec[1],centervec[2]};

    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ***************************************************************************/

    // set up ebfactory
    EBSupport m_eb_support_level = EBSupport::full;

    amrex::Print() << " " << std::endl;
    amrex::Print() << " Orifice Radius: " << orifice_radius << std::endl;
    amrex::Print() << " Funnel Radius:  " << funnel_radius  << std::endl;
    amrex::Print() << " Funnel Height:  " << funnel_height  << std::endl;
    amrex::Print() << " Direction: "      << direction      << std::endl;
    amrex::Print() << " Center:    "      << center[0]      << ", "
                   << center[1]      << ", "
                   << center[2]      << std::endl;


    // Construct the plane representing the funnel edge
    Array<Real,3> point{0.0, 0.0, 0.0};
    point[0] = orifice_radius;

    Array<Real,3> normal{0.0, 0.0, 0.0};
    normal[0] =  funnel_height;
    normal[1] = -orifice_radius;

    amrex::Print() << " " << std::endl;
    amrex::Print() << " Funnel plane: " << std::endl;
    amrex::Print() << "   Point:  " << point[0]  << ", "
                   << point[1]  << ", "
                   << point[2]  << std::endl;

    amrex::Print() << "   Normal: " << normal[0] << ", "
                   << normal[1] << ", "
                   << normal[2] << std::endl;

    EB2::PlaneIF funnel(point, normal);


    // Construct the plane representing the hopper-body's edge
    point = {funnel_radius, 0., 0.};
    normal = {1.0, 0.0, 0.0};

    amrex::Print() << " " << std::endl;
    amrex::Print() << "Hopper body plane: " << std::endl;
    amrex::Print() << "   Point:  " << point[0]  << ", "
                   << point[1]  << ", "
                   << point[2]  << std::endl;

    amrex::Print() << "   Normal: " << normal[0] << ", "
                   << normal[1] << ", "
                   << normal[2] << std::endl;

    EB2::PlaneIF hbody(point, normal);

    // Basic hopper shape (by lathing the union of the body and hopper edges)
    auto hopper0 = EB2::lathe(EB2::makeUnion(hbody, funnel));

    // Align to correct axis
    auto hopper1 = EB2::rotate(hopper0, 2.*std::atan(1.0), 1);

    // Translate to correct location
    auto my_hopper = EB2::translate(hopper1, center);

    // Construct EB2 Index Space
    amrex::Print() << "Building the hopper geometry ..." << std::endl;

    auto gshop = EB2::makeShop(my_hopper);
    int max_coarsening_level = 100;
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);

    const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();

    for (int lev = 0; lev < nlev; lev++)
    {

        eb_level_fluid     = & eb_is.getLevel(geom[lev]);
        eb_level_particles =   eb_level_fluid;


        amrex::Print() << "Done building the hopper geometry" << std::endl;

        /************************************************************************
         *                                                                      *
         * THIS FILLS PARTICLE EBFACTORY                                        *
         *                                                                      *
         ***********************************************************************/
        if (solve_dem)
        {
            amrex::Print() << " " << std::endl;
            amrex::Print() << "Now  making the particle ebfactory ..." << std::endl;

            particle_ebfactory[lev].reset(
                new EBFArrayBoxFactory(* eb_level_particles, geom[lev], grids[lev], dmap[lev],
                                       {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                        m_eb_full_grow_cells},m_eb_support_level)
                );

            eb_normals[lev]->define(grids[lev], dmap[lev], 3, 2, MFInfo(), *particle_ebfactory[lev]);
            amrex::FillEBNormals( * eb_normals[lev], * particle_ebfactory[lev], geom[lev]);


            /********************************************************************
             *                                                                  *
             * Fill level-set:                                                  *
             *                                                                  *
             *******************************************************************/

            if (!levelset__restart) {
                amrex::Print() << "Creating the levelset ..." << std::endl;

                // neat trick: `decltype(my_hopper)` returns the IF data type of
                // the hopper (which otherwise would be a long template string
                // of compounded IFs ...). Together with `auto`, this lets use
                // be more suspect.
                GShopLSFactory<decltype(my_hopper)> gshop_lsfactory(gshop, * level_set);
                std::unique_ptr<MultiFab> mf_impfunc = gshop_lsfactory.fill_impfunc();

                // Construct EB2 Index space based on the refined geometry
                // (level_set->get_eb_geom()). The IndexSpace's geometry needs to
                // match the one used by the eb_factory later.

                auto gshop = EB2::makeShop(my_hopper);
                int max_coarsening_level = 100;
                EB2::Build(gshop, level_set->get_eb_geom(), max_level_here,
                           max_level_here + max_coarsening_level);

                const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
                const EB2::Level & eb_is_ref  = eb_is.getLevel(level_set->get_eb_geom());

                // Construct EBFABFactory based on the refined EB geometry (built above).
                int eb_grow = level_set->get_eb_pad();
                EBFArrayBoxFactory eb_factory(eb_is_ref,
                                              level_set->get_eb_geom(),
                                              level_set->get_eb_ba(),
                                              level_set->get_dm(),
                                              {eb_grow, eb_grow, eb_grow}, EBSupport::full);

                level_set->intersection_ebf(eb_factory, * mf_impfunc);

                amrex::Print() << "Done making the levelset ..." << std::endl;
            } else {
                amrex::Print() << "Loaded level-set is fine => skipping levelset calculation."
                               << std::endl;
            }

            amrex::Print() << "Done making the particle ebfactory ..." << std::endl;
            amrex::Print() << " " << std::endl;
        }

        /************************************************************************
         *                                                                      *
         * THIS FILLS FLUID EBFACTORY                                           *
         *                                                                      *
         ***********************************************************************/

        if (solve_fluid)
        {
            amrex::Print() << "Now  making the fluid ebfactory ..." << std::endl;

            ebfactory[lev].reset(
                new EBFArrayBoxFactory(* eb_level_fluid, geom[lev], grids[lev], dmap[lev],
                                       {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                        m_eb_full_grow_cells}, m_eb_support_level       )
                );

            amrex::Print() << "Done making the fluid ebfactory ..." << std::endl;
        }
    }
}

void mfix::make_amr_hopper()
{
    ParmParse pp("hopper");

    /****************************************************************************
     * Get hopper information from inputs file.                                 *
     ***************************************************************************/

    Real orifice_radius  = 0.00025;
    Real funnel_radius   = 0.00050;
    Real funnel_height   = 0.00105;

    int direction = 0;
    Vector<Real> centervec(3);

    pp.query("orifice_radius", orifice_radius);
    pp.query("funnel_radius",  funnel_radius);
    pp.query("funnel_height",  funnel_height);

    pp.query("direction", direction);
    pp.getarr("center",   centervec,  0, 3);
    Array<Real,3> center={centervec[0],centervec[1],centervec[2]};

    amrex::Print() << " " << std::endl;
    amrex::Print() << " Orifice Radius: " << orifice_radius << std::endl;
    amrex::Print() << " Funnel Radius:  " << funnel_radius  << std::endl;
    amrex::Print() << " Funnel Height:  " << funnel_height  << std::endl;
    amrex::Print() << " Direction: "      << direction      << std::endl;
    amrex::Print() << " Center:    "      << center[0]      << ", "
                   << center[1]      << ", "
                   << center[2]      << std::endl;


    // Construct the plane representing the funnel edge
    Array<Real,3> point{0.0, 0.0, 0.0};
    point[0] = orifice_radius;

    Array<Real,3> normal{0.0, 0.0, 0.0};
    normal[0] =  funnel_height;
    normal[1] = -orifice_radius;

    amrex::Print() << " " << std::endl;
    amrex::Print() << " Funnel plane: " << std::endl;
    amrex::Print() << "   Point:  " << point[0]  << ", "
                   << point[1]  << ", "
                   << point[2]  << std::endl;

    amrex::Print() << "   Normal: " << normal[0] << ", "
                   << normal[1] << ", "
                   << normal[2] << std::endl;

    EB2::PlaneIF funnel(point, normal);


    // Construct the plane representing the hopper-body's edge
    point = {funnel_radius, 0., 0.};
    normal = {1.0, 0.0, 0.0};

    amrex::Print() << " " << std::endl;
    amrex::Print() << "Hopper body plane: " << std::endl;
    amrex::Print() << "   Point:  " << point[0]  << ", "
                   << point[1]  << ", "
                   << point[2]  << std::endl;

    amrex::Print() << "   Normal: " << normal[0] << ", "
                   << normal[1] << ", "
                   << normal[2] << std::endl;

    EB2::PlaneIF hbody(point, normal);

    // Basic hopper shape (by lathing the union of the body and hopper edges)
    auto hopper0 = EB2::lathe(EB2::makeUnion(hbody, funnel));

    // Align to correct axis
    auto hopper1 = EB2::rotate(hopper0, 2.*std::atan(1.0), 1);

    // Translate to correct location
    auto my_hopper = EB2::translate(hopper1, center);

    // Construct EB2 Index Space
    amrex::Print() << "Building the hopper geometry ..." << std::endl;

    auto gshop = EB2::makeShop(my_hopper);


    /****************************************************************************
     * Construct AMR level-set.                                                 *
     ***************************************************************************/

    if (use_amr_ls)
    {
        int lev_lowest = 0;

        const RealBox & rb = geom[lev_lowest].ProbDomain();
        Box domain = geom[lev_lowest].Domain();
        domain.coarsen(amr_ls_crse);

        const IntVect & dom_lo = domain.smallEnd();
        const IntVect & dom_hi = domain.bigEnd();
        // Picket-fence principle
        IntVect n_cells = dom_hi - dom_lo + IntVect{1, 1, 1};
        Vector<int> v_cells = {
            AMREX_D_DECL(n_cells[0], n_cells[1], n_cells[2])
        };


        amrex::Print() << "Declaring AMR levelset:" << std::endl
                       << "coarsest level: " << domain << " n_cells: " << n_cells << std::endl;

        amr_level_set.reset(new LSCore<decltype(my_hopper)>(gshop, & rb, amr_ls_max_level, v_cells));

        amrex::Print() << "... done declaring AMR levelset" << std::endl;
    }
}
