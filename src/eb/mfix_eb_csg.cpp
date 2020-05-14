#include <mfix.H>

#include <csg.hpp>

void mfix::make_eb_csg(std::string geom_file) {
  auto csg_if = csg::get_csgif_from_filename(geom_file);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(csg_if, "Unable to create CsgIF from geometry file");
  auto gshop = EB2::makeShop(*csg_if);
  build_eb_levels(gshop);
}
