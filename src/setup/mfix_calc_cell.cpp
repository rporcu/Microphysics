//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: calc_cell_ic                                            !
//  Purpose: calculate the i, j or k cell index for IC regions.         !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

#include <AMReX_Math.H>
#include <mfix_calc_cell.H>


using namespace amrex;

void calc_cell_ic(const Real dx,
                  const Real dy,
                  const Real dz,
                  const Real* lo,
                  const Real* hi,
                  const Real* plo,
                  int& i_w,
                  int& i_e,
                  int& j_s,
                  int& j_n,
                  int& k_b,
                  int& k_t)
{
  i_w = static_cast<int>(amrex::Math::floor((lo[0]-plo[0])/dx + .5));
  i_e = static_cast<int>(amrex::Math::floor((hi[0]-plo[0])/dx + .5)) - 1;

  j_s = static_cast<int>(amrex::Math::floor((lo[1]-plo[1])/dy + .5));
  j_n = static_cast<int>(amrex::Math::floor((hi[1]-plo[1])/dy + .5)) - 1;

  k_b = static_cast<int>(amrex::Math::floor((lo[2]-plo[2])/dz + .5));
  k_t = static_cast<int>(amrex::Math::floor((hi[2]-plo[2])/dz + .5)) - 1;
}
