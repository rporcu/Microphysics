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
mfix::make_eb_air_reactor()
{
    ParmParse pp("airreactor");

    Real reactor_radius  = 0.0762;
    Real reactor_height  = 0.7720;
    Real riser_radius    = 0.03175;
    Real riser_height    = 3.70840;
    Real outlet_height   = 3.6322;
    Real c2c_height      = 0.0762;

    pp.query("radius0", reactor_radius);
    pp.query("height0", reactor_height);

    pp.query("c2c",     c2c_height    );

    pp.query("radius1", riser_radius  );
    pp.query("height1", riser_height  );

    pp.query("outlet",  outlet_height );

    // set up ebfactory

    Array<Real,3> point, normal, center;


    /****************************************************************************
     * Build the air reactor's reactor section                                  *
     ***************************************************************************/

    // Define point+normal to define a cylinder
    point  = {reactor_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF reactor_body(point, normal);

    // Define point+normal to define the cone cap
    point  = {reactor_radius, reactor_height, 0.0};
    normal = {c2c_height, reactor_radius - riser_radius, 0.0};

    EB2::PlaneIF reactor_cap(point, normal);

    // Create reactor with cap
    auto reactor_wcap = EB2::makeUnion(reactor_body, reactor_cap);



    /****************************************************************************
     * Build the air reactor's riser section                                    *
     ***************************************************************************/

    // Define point+normal to define a cylinder
    point  = {riser_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF riser_body(point, normal);

    // Define point+normal to define a reactor top
    point  = {0.0, riser_height, 0.0};
    normal = {0.0, 1.0, 0.0};

    EB2::PlaneIF riser_top(point, normal);

    // Create reactor by clipping off top
    auto riser = EB2::makeUnion(riser_body, riser_top);


    /****************************************************************************
     * Combine the reactor and riser sections and make 3D                       *
     ***************************************************************************/

    auto airreactor0 = EB2::lathe(EB2::makeIntersection(reactor_wcap, riser));

    // Align to correct axis
    auto airreactor1 = EB2::rotate(airreactor0, 2.0*std::atan(1.0), 1);

    // Translate to correct location
    center = {0.0, 0.1392, 0.1392};
    auto airreactor2 = EB2::translate(airreactor1, center);


    /****************************************************************************
     * Build the outflow pipe and connect to the air reactor                    *
     ***************************************************************************/

    // Define point+normal to define a cylinder
    point  = {riser_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF outlet_body(point, normal);

    // Define point+normal to define a reactor top
    point  = {0.0, 0.1392, 0.0};
    normal = {0.0, -1.0, 0.0};

    EB2::PlaneIF outlet_side(point, normal);

    auto outlet1 = EB2::lathe(EB2::makeUnion(outlet_body, outlet_side));

    // Translate to correct location
    center = {outlet_height, 0.1392, 0.0};
    auto outlet2 = EB2::translate(outlet1, center);


    /****************************************************************************
     * Attached the outflow pipe to the air-reactor's riser                     *
     ***************************************************************************/

    auto airreactor = EB2::makeIntersection(airreactor2, outlet2);


    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ***************************************************************************/

    // Construct EB2 Index Space
    auto gshop_cyc = EB2::makeShop(airreactor);

    build_eb_levels(gshop_cyc);



    Print() << "Done making the fluid eb levels ..." << std::endl;
}
