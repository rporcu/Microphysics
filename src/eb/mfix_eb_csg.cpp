#include <mfix.H>

#include <csg.hpp>

void mfix::make_eb_csg(std::string geom_file) {
  bool is_internal_flow = true;
  ParmParse pp("csg");
  pp.query("internal_flow", is_internal_flow);
  amrex::Print() << "\n Building geometry with is_internal_flow:  " << is_internal_flow << std::endl;

  auto csg_if = csg::get_csgif(geom_file, is_internal_flow);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(csg_if, "Unable to create CsgIF from geometry file");
  auto gshop = EB2::makeShop(*csg_if);
  build_eb_levels(gshop);
}
