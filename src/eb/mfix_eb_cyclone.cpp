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
    Print() << "Building the cyclone geometry ..." << std::endl;

    auto cyc_if = EB2::makeUnion(my_cyclone, false_bottom);

    // auto gshop_cyc = EB2::makeShop(EB2::makeUnion(my_cyclone, false_bottom));
    auto gshop_cyc = EB2::makeShop(cyc_if);

    build_eb_levels(gshop_cyc);

    //_______________________________________________________________________
    // Particles need the correct volfrac at the inflow
    bool has_walls = false;
    std::unique_ptr<UnionListIF<EB2::PlaneIF>> walls = get_walls(has_walls);
    if (has_walls)
    {
        auto if_part = EB2::makeUnion(cyc_if, * walls);
        auto gshop_part = EB2::makeShop(if_part);

        build_particle_eb_levels(gshop_part);
    }

    Print() << "Done making the fluid eb levels ..." << std::endl;
}
