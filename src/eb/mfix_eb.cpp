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

  ParmParse pp("mfix");

  std::string geom_type;
  pp.query("geometry", geom_type);

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
