#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

//#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)
//#include <sstream>

#include <algorithm>
#include <AMReX_EB_levelset.H>
#include <mfix_level.H>
#include <mfix_eb_F.H>

void
mfix_level::make_eb_geometry (int lev)
{

  /******************************************************************************
   * mfix.geometry=<string> specifies the EB geometry. <string> can be on of    *
   * box, cylinder, hopper, clr, clr_riser (or nothing)                         *
   ******************************************************************************/

  ParmParse pp("mfix");

  std::string geom_type;
  pp.query("geometry", geom_type);

  /******************************************************************************
   * Legacy inputs:                                                             *
   *   -- mfix.hourglass = true <=> mfix.geometry=box                           *
   *   -- mfix.clr       = true <=> mfix.geometry=clr                           *
   *   -- mfix.clr_riser = true <=> mfix.geometry=clr_riser                     *
   ******************************************************************************/

  bool hourglass    = false;
  bool clr          = false;
  bool clr_riser    = false;

  pp.query("hourglass", hourglass);
  pp.query("clr", clr);
  pp.query("clr_riser", clr_riser);

  // Avoid multiple (ambiguous) inputs
  if (hourglass || clr || clr_riser) {
      if (! geom_type.empty()) {
          amrex::Abort("The input file cannot specify both:\n"
                       "mfix.<geom_type>=true and mfix.geometry=<geom_type>\n"
                       "at the same time."                                     );
      }
  }

  if (hourglass) geom_type = "hourglass";
  if (clr)       geom_type = "clr";
  if (clr_riser) geom_type = "clr_riser";


  /******************************************************************************
   *                                                                            *
   *  CONSTRUCT EB                                                              *
   *                                                                            *
   ******************************************************************************/


  if (geom_type == "box") {
    amrex::Print() << "\n Building box geometry." << std::endl;
    make_eb_box(lev);
  }
  else if (geom_type == "cylinder") {
    amrex::Print() << "\n Building cylinder geometry." << std::endl;
    make_eb_cylinder(lev);

  }
  else if (geom_type == "hopper") {
    amrex::Print() << "\n Building hopper geometry." << std::endl;
    make_eb_hopper(lev);
  }
  else {
    amrex::Print() << "\n No EB geometry in this problem." << std::endl;
    make_eb_regular(lev);
  }
}
