#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <mfix_eb_if.H>

#include <algorithm>
#include <mfix.H>

#include <mfix_bc.H>


std::shared_ptr<UnionListIF<EB2::PlaneIF>>
mfix::get_walls (bool & has_walls)
{
  // Extracts all walls from the mfix.dat

  has_walls = (m_boundary_conditions.flow_planes().size() > 0);  // will be set to true if there are any walls

  std::shared_ptr<UnionListIF<EB2::PlaneIF>> ret =
    std::shared_ptr<UnionListIF<EB2::PlaneIF>>(new UnionListIF<EB2::PlaneIF>(m_boundary_conditions.flow_planes()));
  return ret;
}


std::shared_ptr<UnionListIF<EB2::PlaneIF>>
mfix::get_real_walls (bool & has_real_walls)
{
  // Extracts all walls from the mfix.dat

  has_real_walls = (m_boundary_conditions.wall_planes().size() > 0);  // will be set to true if there are any walls

  std::shared_ptr<UnionListIF<EB2::PlaneIF>> ret =
    std::shared_ptr<UnionListIF<EB2::PlaneIF>>(new UnionListIF<EB2::PlaneIF>(m_boundary_conditions.wall_planes()));
  return ret;
}
