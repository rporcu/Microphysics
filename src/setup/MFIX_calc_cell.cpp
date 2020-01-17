//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: calc_cell_ic                                            !
//  Purpose: calculate the i, j or k cell index for IC regions.         !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

#include <MFIX_calc_cell.H>

#include <cmath>

using namespace amrex;

void calc_cell_ic(const Real dx,
                  const Real dy,
                  const Real dz,
                  const amrex::Real* lo,
                  const amrex::Real* hi,
                  int& i_w,
                  int& i_e,
                  int& j_s,
                  int& j_n,
                  int& k_b,
                  int& k_t)
{
  i_w = std::floor(lo[0]/dx + .5);
  i_e = std::floor(hi[0]/dx + .5) - 1;

  j_s = std::floor(lo[1]/dy + .5);
  j_n = std::floor(hi[1]/dy + .5) - 1;

  k_b = std::floor(lo[2]/dz + .5);
  k_t = std::floor(hi[2]/dz + .5) - 1;
}
