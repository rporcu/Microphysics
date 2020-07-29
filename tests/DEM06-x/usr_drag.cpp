#include <mfix_des_drag_K.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE
amrex::Real
ComputeDragUser::operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                             amrex::Real DPM, amrex::Real DPA, amrex::Real PHIS,
                             amrex::Real fvelx, amrex::Real fvely, amrex::Real fvelz,
                             int i, int j, int k, int pid) const
{
    amrex::Real ROg = ROPg / EPg;
    amrex::Real RE = (Mug > 0.0) ? DPM*vrel*ROg/Mug : large_number;

    amrex::Real Cd = 0.0;
    if (RE > eps) Cd = (24.0/RE)*(1.0 + 0.15*std::pow(RE, 0.687));

    return 0.75*(ROg*vrel/DPM)*Cd;
}
