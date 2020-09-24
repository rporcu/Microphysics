#include <mfix.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB2_IF_Scale.H>
#include <AMReX_EB2_IF_Translation.H>

#include <csg.hpp>

void mfix::make_eb_csg(std::string geom_file) {
  bool is_internal_flow = true;
  Vector<Real> scaling_factor_vec(3, 1.0);
  Vector<Real> translation_vec(3, 0.0);

  ParmParse pp("csg");
  pp.query("internal_flow", is_internal_flow);
  if(pp.queryarr("scaling_factor", scaling_factor_vec, 0, 3)) {
    amrex::Print() << "WARNING: The implicit function magnitudes will not be scaled" << std::endl;
  }
  pp.queryarr("translation", translation_vec, 0, 3);
  Array<Real,3> scaling_factor = {scaling_factor_vec[0], scaling_factor_vec[1], scaling_factor_vec[2]};
  Array<Real,3> translation = {translation_vec[0], translation_vec[1], translation_vec[2]};

  amrex::Print() << "\n Building geometry with is_internal_flow:  " << is_internal_flow << std::endl;
  auto csg_if = csg::get_csgif(geom_file, is_internal_flow);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(csg_if, "Unable to create CsgIF from geometry file");

  auto final_csg_if = EB2::translate(EB2::scale(*csg_if, scaling_factor), translation);
  
  auto gshop = EB2::makeShop(final_csg_if);
  build_eb_levels(gshop);
}
