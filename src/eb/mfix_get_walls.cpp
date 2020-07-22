#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <mfix_eb_if.H>

#include <algorithm>
#include <mfix.H>

#include <mfix_bc_parms.H>

using namespace std;

shared_ptr<UnionListIF<EB2::PlaneIF>>
mfix::get_walls (bool & has_walls)
{
  // Extracts all walls from the mfix.dat

  has_walls = (BC::flow_planes.size() > 0);  // will be set to true if there are any walls

  shared_ptr<UnionListIF<EB2::PlaneIF>> ret =
    shared_ptr<UnionListIF<EB2::PlaneIF>>(new UnionListIF<EB2::PlaneIF>(BC::flow_planes));
  return ret;
}


shared_ptr<UnionListIF<EB2::PlaneIF>>
mfix::get_real_walls (bool & has_real_walls)
{
  // Extracts all walls from the mfix.dat

  has_real_walls = (BC::wall_planes.size() > 0);  // will be set to true if there are any walls

  shared_ptr<UnionListIF<EB2::PlaneIF>> ret =
    shared_ptr<UnionListIF<EB2::PlaneIF>>(new UnionListIF<EB2::PlaneIF>(BC::wall_planes));
  return ret;
}
