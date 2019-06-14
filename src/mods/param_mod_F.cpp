#include <param_mod_F.H>

#include <limits>

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

AMREX_GPU_HOST_DEVICE
bool is_equal(const amrex::Real& x, const amrex::Real& y)
{
  const amrex::Real epsilon = std::numeric_limits<amrex::Real>::epsilon();

  return (std::abs(x-y) < epsilon);
}
