#include <scales_mod_F.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE
Real scale_pressure_cpp (const Real XXX, const Real P_ref, const Real P_scale)
{
  return (XXX - P_ref) / P_scale;
}

AMREX_GPU_HOST_DEVICE
Real unscale_pressure_cpp (const Real XXX, const Real P_ref, const Real P_scale)
{
  return (XXX * P_scale + P_ref);
}
