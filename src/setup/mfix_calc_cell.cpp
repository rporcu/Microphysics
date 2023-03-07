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


const amrex::Box calc_ic_box(const Geometry& a_geom, const RealBox* a_region)
{

  const GpuArray<Real,3> dxi = a_geom.InvCellSizeArray();
  const GpuArray<Real,3> plo = a_geom.ProbLoArray();

  const Real* lo = a_region->lo();
  const Real* hi = a_region->hi();

  const amrex::IntVect ic_lo(AMREX_D_DECL(
      static_cast<int>(amrex::Math::floor((lo[0]-plo[0])*dxi[0] + .5)),
      static_cast<int>(amrex::Math::floor((lo[1]-plo[1])*dxi[1] + .5)),
      static_cast<int>(amrex::Math::floor((lo[2]-plo[2])*dxi[2] + .5))));

  const amrex::IntVect ic_hi(AMREX_D_DECL(
      static_cast<int>(amrex::Math::floor((hi[0]-plo[0])*dxi[0] + .5)) - 1,
      static_cast<int>(amrex::Math::floor((hi[1]-plo[1])*dxi[1] + .5)) - 1,
      static_cast<int>(amrex::Math::floor((hi[2]-plo[2])*dxi[2] + .5)) - 1));

  Box box(ic_lo, ic_hi);

  // The IC box starts as the full domain, then we take the min
  // with the box we computed to make sure we stay inside the domain.
  Box ic_box(a_geom.Domain());
  ic_box.minBox(box);

  return ic_box;

}

const amrex::Box calc_bc_box(Geometry const& a_geom,
                             RealBox  const* a_region,
                             int const       a_face)
{
  Box ic_box = calc_ic_box(a_geom, a_region);
  if (a_face == -1) {
    return ic_box;
  }

  Box domain(a_geom.Domain());

  // face: 0:x-lo, 1:x-hi, 2:y-lo, 3:y-hi, 4:z-lo, 5:z-hi
  int const dir((a_face - (a_face%2))/((int)2)); // dir=0,1,2 (x,y,z)

  IntVect lo(domain.smallEnd());
  IntVect hi(domain.bigEnd());

  if ( a_face%2==0 ) {
    hi.setVal(dir,lo[dir]);
  } else {
    lo.setVal(dir,hi[dir]);
  }

  // The bc box starts as the whole side of the domain, then
  // we take the min with the box we computed earlier.
  Box bc_box(lo, hi);
  bc_box.minBox(ic_box);

  return bc_box;
}
