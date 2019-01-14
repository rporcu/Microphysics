#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Rotation.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Complement.H>


#include <algorithm>
#include <AMReX_EB_utils.H>
#include <AMReX_EB_LSCore.H>
#include <AMReX_EB_levelset.H>
#include <mfix.H>
#include <mfix_eb_F.H>


/********************************************************************************
 * Function to create a simple hopper EB.                                       *
 *                                                                              *
 * Comments: The hopper is constructed by "lathing" a union of two planes       *
 * (representing the hopper wall and outlet funnel) around a central axis.      *
 *                                                                              *
 *******************************************************************************/
void
mfix::make_eb_cyclone()
{
    ParmParse pp("hopper");

    int max_level_here = 0;

    // set up ebfactory

    EBSupport m_eb_support_level = EBSupport::full;

    /****************************************************************************
     * Build the cyclone's inlet tube.                                          *
     ***************************************************************************/

    Real xlo = 0.0017, ylo = 0.0020, zlo = 0.0005;
    Real xhi = 0.0035, yhi = 0.0060, zhi = 0.0012;

    Array<Real,3>  point_lox{ xlo, 0.0, 0.0};
    Array<Real,3> normal_lox{-1.0, 0.0, 0.0};
    Array<Real,3>  point_hix{ xhi, 0.0, 0.0};
    Array<Real,3> normal_hix{ 1.0, 0.0, 0.0};

    Array<Real,3>  point_loy{0.0, ylo, 0.0};
    Array<Real,3> normal_loy{0.0,-1.0, 0.0};
    Array<Real,3>  point_hiy{0.0, yhi, 0.0};
    Array<Real,3> normal_hiy{0.0, 1.0, 0.0};

    Array<Real,3>  point_loz{0.0, 0.0, zlo};
    Array<Real,3> normal_loz{0.0, 0.0,-1.0};
    Array<Real,3>  point_hiz{0.0, 0.0, zhi};
    Array<Real,3> normal_hiz{0.0, 0.0, 1.0};

    EB2::PlaneIF plane_lox(point_lox,normal_lox);
    EB2::PlaneIF plane_hix(point_hix,normal_hix);

    EB2::PlaneIF plane_loy(point_loy,normal_loy);
    EB2::PlaneIF plane_hiy(point_hiy,normal_hiy);

    EB2::PlaneIF plane_loz(point_loz,normal_loz);
    EB2::PlaneIF plane_hiz(point_hiz,normal_hiz);

    auto inlet = EB2::makeUnion(plane_lox, plane_hix,
                                plane_loy, plane_hiy,
                                plane_loz, plane_hiz );


    Real radius    = 0.0015;
    Real height    = 0.0040;
    Array<Real,3> center={0.0015, 0.0020, 0.0020};

    int direction  = 0;
    bool inside    = true;


    EB2::CylinderIF main_body(radius, height, direction, center, inside);


    radius    = 0.0008;
    height    = 0.0040;
    center={0.0035, 0.0020, 0.0020};

    auto od = EB2::makeComplement(
        EB2::CylinderIF(radius, height, direction, center, inside));


    radius    = 0.00055;
    height    = 0.00400;
    center={0.0035, 0.0020, 0.0020};

    EB2::CylinderIF id(radius, height, direction, center, inside);


    auto my_cyclone = EB2::makeIntersection(id,
                                            EB2::makeUnion(od,
                                                           EB2::makeIntersection(inlet,main_body)));


    // Create a plane to block the bottom for fluid.

    Array<Real,3> point{0.0, 0.0, 0.0};
    Array<Real,3> normal{0.0, 0.0, 0.0};

    point[direction] = geom[0].ProbLo(direction) + 1.0e-8;
    normal[direction] = -1.0;

    EB2::PlaneIF false_bottom(point, normal);

    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ***************************************************************************/

    // Construct EB2 Index Space
    amrex::Print() << "Building the cyclone geometry ..." << std::endl;

    auto cyc_if = EB2::makeUnion(my_cyclone, false_bottom);

    // auto gshop_cyc = EB2::makeShop(EB2::makeUnion(my_cyclone, false_bottom));
    auto gshop_cyc = EB2::makeShop(cyc_if);

    build_eb_levels(gshop_cyc);

    // int max_coarsening_level = 100;
    // EB2::Build(gshop_cyc, geom.back(), max_level_here, max_level_here + max_coarsening_level);

    // const EB2::IndexSpace & ebis_cyc = EB2::IndexSpace::top();

    // for (int lev = 0; lev < nlev; lev++)
    // {

    //     const EB2::Level & ebis_lev_cyc  = ebis_cyc.getLevel(geom[lev]);

    //     amrex::Print() << "Done building the cyclone geometry" << std::endl;

    //     /************************************************************************
    //      *                                                                      *
    //      * THIS FILLS PARTICLE EBFACTORY                                        *
    //      * NOTE: the "bottom" plane is only applied to particles => fluid       *
    //      *       EBFactory only needs the cyclone walls.                        *
    //      *                                                                      *
    //      ***********************************************************************/

    //     if (solve_dem)
    //     {
    //         amrex::Print() << " " << std::endl;
    //         amrex::Print() << "Now making the particle ebfactory ..." << std::endl;

    //         auto gshop = EB2::makeShop(my_cyclone);
    //         EB2::Build(gshop, geom.back(), max_level_here,
    //                    max_level_here + max_coarsening_level);

    //         // intercept the cyclone + plane level (if it is built)
    //         const EB2::IndexSpace & eb_is = EB2::IndexSpace::top();
    //         eb_level_particles = & eb_is.getLevel(geom[lev]);

    //         particle_ebfactory[lev].reset(
    //             new EBFArrayBoxFactory(* eb_level_particles, geom[lev], grids[lev], dmap[lev],
    //                                    {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
    //                                     m_eb_full_grow_cells}, m_eb_support_level)
    //             );

    //         /********************************************************************
    //          *                                                                  *
    //          * Fill level-set:                                                  *
    //          *                                                                  *
    //          *******************************************************************/

    //         if (!levelset__restart)
    //         {
    //             amrex::Print() << "Creating the levelset ..." << std::endl;

    //             GShopLSFactory<decltype(my_cyclone)> cyc_lsfactory(gshop, * level_set);
    //             std::unique_ptr<MultiFab> mf_impfunc_cyc = cyc_lsfactory.fill_impfunc();

    //             // Construct EB2 Index space based on the refined geometry
    //             // (level_set->get_eb_geom()). The IndexSpace's geometry needs to
    //             // match the one used by the eb_factory later.

    //             EB2::Build(gshop_cyc, level_set->get_eb_geom(), max_level_here,
    //                        max_level_here + max_coarsening_level);

    //             const EB2::IndexSpace & ebis_cyc    = EB2::IndexSpace::top();
    //             const EB2::Level & ebis_lev_cyc_ref = ebis_cyc.getLevel(level_set->get_eb_geom());

    //             // Construct EBFABFactory based on the refined EB geometry (built above).
    //             int eb_grow = level_set->get_eb_pad();
    //             EBFArrayBoxFactory eb_factory_cyclone(ebis_lev_cyc_ref,
    //                                                   level_set->get_eb_geom(),
    //                                                   level_set->get_eb_ba(),
    //                                                   level_set->get_dm(),
    //                                                   {eb_grow, eb_grow, eb_grow}, EBSupport::full);

    //             level_set->Intersect(eb_factory_cyclone, * mf_impfunc_cyc );

    //             amrex::Print() << "Done making the levelset ..." << std::endl;

    //         } else {
    //             amrex::Print() << "Loaded level-set is fine => skipping levelset calculation."
    //                            << std::endl;
    //         }

    //         amrex::Print() << "Done making the particle ebfactories ..." << std::endl;
    //         amrex::Print() << " " << std::endl;
    //     }

    //     /************************************************************************
    //      *                                                                      *
    //      * THIS FILLS FLUID EBFACTORY                                           *
    //      *                                                                      *
    //      ***********************************************************************/

    //     if (solve_fluid)
    //     {
    //         amrex::Print() << "Now making the fluid ebfactory ..." << std::endl;

    //         eb_level_fluid = & ebis_lev_cyc;

    //         ebfactory[lev].reset(
    //             new EBFArrayBoxFactory(* eb_level_fluid, geom[lev], grids[lev], dmap[lev],
    //                                  {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
    //                                   m_eb_full_grow_cells}, m_eb_support_level)
    //             );

    //         amrex::Print() << "Done making the fluid ebfactory ..." << std::endl;
    //     }
    // } // lev
}


void mfix::make_amr_cyclone()
{
    ParmParse pp("hopper");

    /****************************************************************************
     * Build the cyclone's inlet tube.                                          *
     ***************************************************************************/

    Real xlo = 0.0017, ylo = 0.0020, zlo = 0.0005;
    Real xhi = 0.0035, yhi = 0.0060, zhi = 0.0012;

    Array<Real,3>  point_lox{ xlo, 0.0, 0.0};
    Array<Real,3> normal_lox{-1.0, 0.0, 0.0};
    Array<Real,3>  point_hix{ xhi, 0.0, 0.0};
    Array<Real,3> normal_hix{ 1.0, 0.0, 0.0};

    Array<Real,3>  point_loy{0.0, ylo, 0.0};
    Array<Real,3> normal_loy{0.0,-1.0, 0.0};
    Array<Real,3>  point_hiy{0.0, yhi, 0.0};
    Array<Real,3> normal_hiy{0.0, 1.0, 0.0};

    Array<Real,3>  point_loz{0.0, 0.0, zlo};
    Array<Real,3> normal_loz{0.0, 0.0,-1.0};
    Array<Real,3>  point_hiz{0.0, 0.0, zhi};
    Array<Real,3> normal_hiz{0.0, 0.0, 1.0};

    EB2::PlaneIF plane_lox(point_lox,normal_lox);
    EB2::PlaneIF plane_hix(point_hix,normal_hix);

    EB2::PlaneIF plane_loy(point_loy,normal_loy);
    EB2::PlaneIF plane_hiy(point_hiy,normal_hiy);

    EB2::PlaneIF plane_loz(point_loz,normal_loz);
    EB2::PlaneIF plane_hiz(point_hiz,normal_hiz);

    auto inlet = EB2::makeUnion(plane_lox, plane_hix,
                                plane_loy, plane_hiy,
                                plane_loz, plane_hiz );


    Real radius    = 0.0015;
    Real height    = 0.0040;
    Array<Real,3> center={0.0015, 0.0020, 0.0020};

    int direction  = 0;
    bool inside    = true;


    EB2::CylinderIF main_body(radius, height, direction, center, inside);


    radius    = 0.0008;
    height    = 0.0040;
    center={0.0035, 0.0020, 0.0020};

    auto od = EB2::makeComplement(
        EB2::CylinderIF(radius, height, direction, center, inside));


    radius    = 0.00055;
    height    = 0.00400;
    center={0.0035, 0.0020, 0.0020};

    EB2::CylinderIF id(radius, height, direction, center, inside);


    auto my_cyclone = EB2::makeIntersection(id,
                                            EB2::makeUnion(od,
                                                           EB2::makeIntersection(inlet,main_body)));


    // Create a plane to block the bottom for fluid.

    Array<Real,3> point{0.0, 0.0, 0.0};
    Array<Real,3> normal{0.0, 0.0, 0.0};

    point[direction] = geom[0].ProbLo(direction) + 1.0e-8;
    normal[direction] = -1.0;

    EB2::PlaneIF false_bottom(point, normal);

    auto impfunc = EB2::makeUnion(my_cyclone, false_bottom);
    auto gshop = EB2::makeShop(impfunc);

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

        amr_level_set.reset(new LSCore<decltype(impfunc)>(gshop, & rb, amr_ls_max_level, v_cells));

        amrex::Print() << "... done declaring AMR levelset" << std::endl;
    }
}

