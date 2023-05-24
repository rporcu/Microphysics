#include <AMReX_Utility.H>

#include <mfix_des_drag_K.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE
Real
ComputeDragUser::operator() (Real EPg, Real Mug, Real ROPg, Real vrel,
                             Real DPM, Real /*DPA*/, Real /*PHIS*/,
                             Real /*fvelx*/, Real /*fvely*/, Real /*fvelz*/,
                             int /*i*/, int /*j*/, int /*k*/, int /*pid*/) const
{
    amrex::Real ROg = ROPg / EPg;
    amrex::Real RE = (Mug > 0.0) ? DPM*vrel*ROg/Mug : large_number;

    amrex::Real Cd = 0.0;
    if (RE > eps) Cd = (24.0/RE)*(1.0 + 0.15*std::pow(RE, 0.687));

    return 0.75*(ROg*vrel/DPM)*Cd;
}
