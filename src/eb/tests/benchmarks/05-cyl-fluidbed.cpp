#include "catch2/catch.hpp"

#include <AMReX_EB2_IF_Cylinder.H>

#include "main.hpp"

namespace {

using namespace amrex;

const Real radius = 0.00045;
const Real Cx = 0.002;
const Real Cy = 0.0005;
const Real Cz = 0.0005;

void check(const csg::CsgIF &my_if, const EB2::CylinderIF &their_cyl,
           std::tuple<double, double, double> point) {
  auto xx = std::get<0>(point);
  auto yy = std::get<1>(point);
  auto zz = std::get<2>(point);
  CHECK(Approx(my_if(xx, yy, zz)) == their_cyl(xx, yy, zz));
}

TEST_CASE("05-cyl-fluidbed", "[benchmarks]") {
  Real height = -1.; // Infinite length cylinder
  int direction = 0;

  Array<Real, 3> center = {Cx, Cy, Cz};
  EB2::CylinderIF their_cyl(radius, height, direction, center, true);

  auto csg_str = std::string(R"(
difference() {
	cube(size = [3e+06, 3e+06, 3e+06], center = true);
	multmatrix([[0, 0, 1, 0.002], [0, 1, 0, 0.0005], [-1, 0, 0, 0.0005], [0, 0, 0, 1]]) {
		cylinder(h = 1e+06, r1 = 0.00045, r2 = 0.00045, center = true);
	}
}
  )");
  auto maybe_tree = csg::parse_csg(csg_str);
  REQUIRE(maybe_tree);
  csg::CsgIF my_if(maybe_tree);

  SECTION("Inside the hole") {
    check(my_if, their_cyl, {Cx + 1000, Cy, Cz});
    check(my_if, their_cyl, {Cx, Cy, Cz});
    check(my_if, their_cyl, {Cx - 1000, Cy, Cz});

    check(my_if, their_cyl, {Cx, Cy + radius * 0.7, Cz - radius * 0.7});
    check(my_if, their_cyl, {Cx, Cy - radius * 0.7, Cz + radius * 0.7});
    check(my_if, their_cyl, {Cx, Cy, Cz + radius * 0.9});
    check(my_if, their_cyl, {Cx, Cy - radius * 0.9, Cz});
  }

  SECTION("Outside the hole") {
    check(my_if, their_cyl, {Cx, Cy + radius * 0.9, Cz - radius * 0.9});
    check(my_if, their_cyl, {Cx, Cy - radius * 0.9, Cz + radius * 0.9});
    check(my_if, their_cyl, {Cx, Cy, Cz + radius * 1.1});
    check(my_if, their_cyl, {Cx, Cy - radius * 1.1, Cz});
  }
}

TEST_CASE("05-cyl-fluidbed-finite", "[benchmarks]") {
  Real height = 1.0;
  int direction = 0;

  Array<Real, 3> center = {Cx, Cy, Cz};
  EB2::CylinderIF their_cyl(radius, height, direction, center, true);

  SECTION("Using difference") {
    auto csg_str = std::string(R"(
difference() {
	cube(size = [3e+06, 3e+06, 3e+06], center = true);
	multmatrix([[0, 0, 1, 0.002], [0, 1, 0, 0.0005], [-1, 0, 0, 0.0005], [0, 0, 0, 1]]) {
		cylinder(h = 1.0, r1 = 0.00045, r2 = 0.00045, center = true);
	}
}
  )");
    auto maybe_tree = csg::parse_csg(csg_str);
    REQUIRE(maybe_tree);
    csg::CsgIF my_if(maybe_tree);

    SECTION("Inside the hole") {
      check(my_if, their_cyl, {Cx + 0.49, Cy, Cz});
      check(my_if, their_cyl, {Cx - 0.49, Cy, Cz});
      check(my_if, their_cyl, {Cx, Cy, Cz});

      check(my_if, their_cyl, {Cx, Cy + radius * 0.7, Cz - radius * 0.7});
      check(my_if, their_cyl, {Cx, Cy - radius * 0.7, Cz + radius * 0.7});
      check(my_if, their_cyl, {Cx, Cy, Cz + radius * 0.9});
      check(my_if, their_cyl, {Cx, Cy - radius * 0.9, Cz});
    }

    SECTION("Outside the hole") {
      check(my_if, their_cyl, {Cx + 0.51, Cy, Cz});
      check(my_if, their_cyl, {Cx - 0.51, Cy, Cz});
      check(my_if, their_cyl, {Cx, Cy + radius * 0.9, Cz - radius * 0.9});
      check(my_if, their_cyl, {Cx, Cy - radius * 0.9, Cz + radius * 0.9});
      check(my_if, their_cyl, {Cx, Cy, Cz + radius * 1.1});
      check(my_if, their_cyl, {Cx, Cy - radius * 1.1, Cz});
    }
  }

  SECTION("Using complement") {
    auto csg_str = std::string(R"(
	multmatrix([[0, 0, 1, 0.002], [0, 1, 0, 0.0005], [-1, 0, 0, 0.0005], [0, 0, 0, 1]]) {
		cylinder(h = 1.0, r1 = 0.00045, r2 = 0.00045, center = true);
	}
  )");
    auto maybe_tree = csg::parse_csg(csg_str);
    REQUIRE(maybe_tree);
    csg::CsgIF my_if(maybe_tree, true);

    SECTION("Inside the hole") {
      check(my_if, their_cyl, {Cx + 0.49, Cy, Cz});
      check(my_if, their_cyl, {Cx - 0.49, Cy, Cz});
      check(my_if, their_cyl, {Cx, Cy, Cz});

      check(my_if, their_cyl, {Cx, Cy + radius * 0.7, Cz - radius * 0.7});
      check(my_if, their_cyl, {Cx, Cy - radius * 0.7, Cz + radius * 0.7});
      check(my_if, their_cyl, {Cx, Cy, Cz + radius * 0.9});
      check(my_if, their_cyl, {Cx, Cy - radius * 0.9, Cz});
    }

    SECTION("Outside the hole") {
      check(my_if, their_cyl, {Cx + 0.51, Cy, Cz});
      check(my_if, their_cyl, {Cx - 0.51, Cy, Cz});
      check(my_if, their_cyl, {Cx, Cy + radius * 0.9, Cz - radius * 0.9});
      check(my_if, their_cyl, {Cx, Cy - radius * 0.9, Cz + radius * 0.9});
      check(my_if, their_cyl, {Cx, Cy, Cz + radius * 1.1});
      check(my_if, their_cyl, {Cx, Cy - radius * 1.1, Cz});
    }
  }
}

} // namespace
