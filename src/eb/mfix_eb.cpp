#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <AMReX_GeometryShop.H>
#include <AMReX_SphereIF.H>
#include <AMReX_PlaneIF.H>
#include <AMReX_AllRegularService.H>
#include <AMReX_FlatPlateGeom.H>
#include <AMReX_EBISLayout.H>
#include <AMReX_EBGraph.H>
#include <AMReX_EBDebugOut.H>
#include <AMReX_EBCellFAB.H>
#include <AMReX_EBCellFactory.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_UnionIF.H>
#include <AMReX_TransformIF.H>
#include <AMReX_ComplementIF.H>
#include <AMReX_IntersectionIF.H>
#include <AMReX_LatheIF.H>
#include <AMReX_PolynomialIF.H>
#include <AMReX_AnisotropicDxPlaneIF.H>
#include <AMReX_AnisotropicIF.H>

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
