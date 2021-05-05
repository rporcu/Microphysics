#include "mfix_pc.H"

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real solids_stress (const amrex::Real Ps,
                           const amrex::Real beta,
                           const amrex::Real ep_cp,
                           const amrex::Real small_number,
                           const amrex::Real eps)
{

  return Ps*std::pow(eps,beta) /
    amrex::max(ep_cp-eps, small_number*(1.-eps));
}
