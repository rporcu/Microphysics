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
mfix::make_eb_proto_clr()
{
    ParmParse pp("airreactor");

    // This is used a padding on for the Y and Z directions.
    Real offset          = 0.16;

    // Air-reactor parameters
    Real reactor_radius  = 0.10;
    Real reactor_height  = 0.75;

    // Riser parameters
    Real riser_radius    = 0.05;
    Real riser_height    = 3.75;

    // Height of the cone connecting the air-reactor and riser
    Real c2c_height      = 0.10;

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
    center = {0.0, offset, offset};
    auto airreactorIF = EB2::translate(airreactor1, center);


    /****************************************************************************
     * Build a horizontal pipe. This is used to connect the riser to the        *
     * cyclone as well as a portion of the l-valve.                             *
     ***************************************************************************/

    // Riser parameters
    Real pipe_radius    = 0.05;

    // Define point+normal to define a cylinder
    point  = {pipe_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF pipeIF_body(point, normal);

    // Define point+normal to define a reactor top
    point  = {0.0, offset, 0.0};
    normal = {0.0, -1.0, 0.0};

    EB2::PlaneIF pipeIF_lo(point, normal);

    // Define point+normal to define a reactor top
    point  = {0.0, offset+0.636, 0.0};
    normal = {0.0, 1.0, 0.0};

    EB2::PlaneIF pipeIF_hi(point, normal);


    auto pipeIF = EB2::lathe(EB2::makeUnion(pipeIF_body,
                                            pipeIF_lo,
                                            pipeIF_hi));


    /****************************************************************************
     * Build crossover to connect air-reactor riser to the cyclone.             *
     ***************************************************************************/

    Real crossover_height   = 3.60;

    // Translate to correct location
    center = {crossover_height, offset, 0.0};

    auto crossoverIF = EB2::translate(pipeIF, center);


    /****************************************************************************
     * Build the cyclone                                                        *
     ***************************************************************************/

    Real cyclone_top              = 3.80;
    Real cyclone_radius           = 0.10;
    Real cyclone_height           = 0.35;
    Real cyclone_dipleg_radius    = 0.04;
    Real cyclone_dipleg_bottom    = 0.15;
    Real cyclone_to_dipleg_height = 0.20;

    Real cyclone_bottom = cyclone_top - cyclone_height;

    amrex::Print() << "cyclone bottom " << cyclone_bottom << std::endl;

    // Define point+normal to define the top of the cyclone
    point  = {0.0, cyclone_top, 0.0};
    normal = {0.0, 1.0, 0.0};

    EB2::PlaneIF cyclone_cap(point, normal);

    // Define point+normal to define a cylinder
    point  = {cyclone_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF cyclone_body(point, normal);

    // Define point+normal to define the cone cap
    point  = {cyclone_radius, cyclone_bottom, 0.0};
    normal = {cyclone_to_dipleg_height, cyclone_dipleg_radius - cyclone_radius, 0.0};

    EB2::PlaneIF cyclone_c2c(point, normal);

    // Create reactor with cap
    auto cyclone_wc2c = EB2::makeUnion(cyclone_body, cyclone_cap, cyclone_c2c);

    // Define point+normal to define the top of the cyclone dipleg
    point  = {0.0, cyclone_bottom, 0.0};
    normal = {0.0, 1.0, 0.0};

    EB2::PlaneIF cyclone_dipleg_top(point, normal);

    // Define point+normal to define the bottom of the cyclone dipleg
    point  = {0.0, cyclone_dipleg_bottom, 0.0};
    normal = {0.0, -1.0, 0.0};

    EB2::PlaneIF cyclone_dipleg_bottomIF(point, normal);

    // Define point+normal to define the radius of the cyclone dipleg
    point  = {cyclone_dipleg_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF cyclone_dipleg_body(point, normal);

    // Create reactor with cap
    auto cyclone_dipleg = EB2::makeUnion(cyclone_dipleg_bottomIF,
                                         cyclone_dipleg_body);

    auto cyclone_wdipleg = EB2::makeIntersection(cyclone_wc2c,
                                                 cyclone_dipleg);


    /****************************************************************************
     * Combine the reactor and riser sections and make 3D                       *
     ***************************************************************************/

    auto cyclone0 = EB2::lathe(cyclone_wdipleg);

    // Align to correct axis
    auto cyclone1 = EB2::rotate(cyclone0, 2.0*std::atan(1.0), 1);

    // Translate to correct location
    center = {0.0, offset, offset+0.636};
    auto cyclone2 = EB2::translate(cyclone1, center);


    /****************************************************************************
     * Create the loop-seal                                                     *
     ***************************************************************************/

    Real loopseal_radius = 0.12;
    Real loopseal_bottom = 2.24;
    Real loopseal_top    = loopseal_bottom + 0.50;

    // Define point+normal to define the top of the cyclone dipleg
    point  = {0.0, loopseal_top, 0.0};
    normal = {0.0, 1.0, 0.0};

    EB2::PlaneIF loopsealIF_top(point, normal);

    // Define point+normal to define the bottom of the cyclone dipleg
    point  = {0.0, loopseal_bottom, 0.0};
    normal = {0.0, -1.0, 0.0};

    EB2::PlaneIF loopsealIF_bottom(point, normal);

    // Define point+normal to define the radius of the cyclone dipleg
    point  = {loopseal_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF loopsealIF_body(point, normal);

    // Create reactor with cap
    auto loopsealIF = EB2::makeUnion(loopsealIF_bottom,
                                     loopsealIF_top,
                                     loopsealIF_body);

    auto loopseal0 = EB2::lathe(loopsealIF);

    // Align to correct axis
    auto loopseal1 = EB2::rotate(loopseal0, 2.0*std::atan(1.0), 1);

    // Translate to correct location
    center = {0.0, offset, 0.96-offset};
    auto loopseal2 = EB2::translate(loopseal1, center);

    /****************************************************************************
     * Create the loop-seal                                                     *
     ***************************************************************************/

    Real fuelreactor_radius = 0.12;
    Real fuelreactor_bottom = 0.52;
    Real fuelreactor_height = 1.54;
    Real fuelreactor_top    = fuelreactor_bottom + fuelreactor_height;

    // Define point+normal to define the top of the cyclone dipleg
    point  = {0.0, fuelreactor_top, 0.0};
    normal = {0.0, 1.0, 0.0};

    EB2::PlaneIF fuelreactorIF_top(point, normal);

    // Define point+normal to define the bottom of the cyclone dipleg
    point  = {0.0, fuelreactor_bottom, 0.0};
    normal = {0.0, -1.0, 0.0};

    EB2::PlaneIF fuelreactorIF_bottom(point, normal);

    // Define point+normal to define the radius of the cyclone dipleg
    point  = {fuelreactor_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF fuelreactorIF_body(point, normal);

    // Create reactor with cap
    auto fuelreactorIF = EB2::makeUnion(fuelreactorIF_bottom,
                                        fuelreactorIF_top,
                                        fuelreactorIF_body);

    auto fuelreactor0 = EB2::lathe(fuelreactorIF);

    // Align to correct axis
    auto fuelreactor1 = EB2::rotate(fuelreactor0, 2.0*std::atan(1.0), 1);

    // Translate to correct location
    center = {0.0, offset, offset+0.636};
    auto fuelreactor2 = EB2::translate(fuelreactor1, center);

    /****************************************************************************
     * Create the loop-seal                                                     *
     ***************************************************************************/

    center = {0.2, offset, 0.0};
    auto lvalve = EB2::translate(pipeIF, center);


    // auto full_clr = cyclone2;
    auto full_clr = EB2::makeIntersection(airreactorIF,
                                          crossoverIF,
                                          cyclone2,
                                          loopseal2,
                                          fuelreactor2,
                                          lvalve);


    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ***************************************************************************/

    // Construct EB2 Index Space
    // auto gshop_cyc = EB2::makeShop(airreactor);
    auto gshop_cyc = EB2::makeShop(full_clr);




    build_eb_levels(gshop_cyc);



    Print() << "Done making the fluid eb levels ..." << std::endl;
}
