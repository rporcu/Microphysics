#include <AMReX_LO_BCTYPES.H>

#include <mfix.H>
#include <MFIX_BCHandler.H>

#include <algorithm>

using namespace amrex;

//
// Set the boundary condition for diffusion solve
//
// MLMG expects the BC type to be the uniform on each domain wall.
// Since here we allow for BC patches on each wall, we first check that
// the user-provided BCs are uniform, and then return a single BC type for
// each domain wall.
//
void
BCHandler::set_scal_diff_bc (amrex::IntVect& bc_lo,
                             amrex::IntVect& bc_hi,
                             const amrex::IntVect&,
                             const amrex::IntVect&,
                             const int)
{
  //
  // By default, all the BCs are Neumann
  //
  std::fill(bc_lo.begin(), bc_lo.end(), AMREX_LO_NEUMANN);
  std::fill(bc_hi.begin(), bc_hi.end(), AMREX_LO_NEUMANN);

  const int cyclic_x = m_geom[0].isPeriodic(0);
  const int cyclic_y = m_geom[0].isPeriodic(1);
  const int cyclic_z = m_geom[0].isPeriodic(2);

  //
  // BC -- X direction
  //
  if (cyclic_x) {
    bc_lo[0] = AMREX_LO_PERIODIC;
    bc_hi[0] = AMREX_LO_PERIODIC;
  }

  //
  // BC -- Y direction
  //
  if (cyclic_y) {
    bc_lo[1] = AMREX_LO_PERIODIC;
    bc_hi[1] = AMREX_LO_PERIODIC;
  }

  //
  // BC -- Z direction
  //
  if (cyclic_z) {
    bc_lo[2] = AMREX_LO_PERIODIC;
    bc_hi[2] = AMREX_LO_PERIODIC;
  }
}
