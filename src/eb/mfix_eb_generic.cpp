#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>

#include <mfix.H>


/********************************************************************************
 *                                                                              *
 * Function to create a generic EB.                                       *
 *                                                                              *
 *******************************************************************************/

void
mfix::make_eb_generic ()
{

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false, "If you called this, you need to code this");
    // Change to user defined EB geometry definition. This is here so that the
    // code compiles.
    EB2::CylinderIF user_eb(0., 0., 0, {0.,0.,0.}, 0);

    // Uncomment and pass final EB geometry
    //auto gshop_cyc = EB2::makeShop(user_eb);
    //build_eb_levels(gshop_cyc);

    //Print() << "Done making the fluid eb levels ..." << std::endl;
}
