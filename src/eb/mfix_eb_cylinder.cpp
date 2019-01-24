#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <algorithm>
#include <AMReX_EB_levelset.H>
#include <mfix.H>
#include <mfix_eb_F.H>


/********************************************************************************
 *                                                                              *
 * Function to create a simple cylinder EB.                                     *
 *                                                                              *
 * Comments: The EB for the cylinder is open at the top and bottom; if the      *
 * boundary condition is mass inflow then we modify the level set so that       *
 * the particles see it as a solid boundary.                                    *
 *                                                                              *
 *******************************************************************************/
void mfix::make_eb_cylinder()
{
    ParmParse pp("cylinder");

    int max_level_here = 0;

    /****************************************************************************
     * Get cylinder information from inputs file.                               *
     ***************************************************************************/
    bool inside       = true;

    Real radius    = 0.0002;
    Real height    = -1.;

    int direction  = 0;
    Vector<Real> centervec(3);

    pp.query("internal_flow", inside);

    pp.query("radius",     radius);
    pp.query("height",     height);
    pp.query("direction",  direction);
    pp.getarr("center",    centervec,  0, 3);
    Array<Real,3> center={centervec[0], centervec[1], centervec[2]};

    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ***************************************************************************/

    // set up ebfactory
    EBSupport m_eb_support_level = EBSupport::full;

    amrex::Print() << " " << std::endl;
    amrex::Print() << " Internal Flow: " << inside << std::endl;
    amrex::Print() << " Radius:    " << radius    << std::endl;
    amrex::Print() << " Height:    " << height    << std::endl;
    amrex::Print() << " Direction: " << direction << std::endl;
    amrex::Print() << " Center:    " << center[0] << ", "
                   << center[1] << ", "
                   << center[2] << std::endl;


    // Create the cylinder -- used for both fluid and particles
    amrex::Print() << "Building the cylinder (side wall) geometry ..." << std::endl;

    // Build the Cylinder geometry first representing the curved walls (this is
    // always present regardless of user input).
    EB2::CylinderIF my_cyl(radius, height, direction, center, inside);

    auto gshop_cyl = EB2::makeShop(my_cyl);

    build_eb_levels(gshop_cyl);

    //_______________________________________________________________________
    // Particles need the correct volfrac at the inflow
    bool has_walls = false;
    std::unique_ptr<UnionListIF<EB2::PlaneIF>> walls = get_walls(has_walls);
    auto if_part = EB2::makeUnion(my_cyl, * walls);
    auto gshop_part = EB2::makeShop(if_part);

    build_particle_eb_levels(gshop_part);


    Print() << "Done making the fluid eb levels ..." << std::endl;
}
