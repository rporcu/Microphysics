#include <mfix_des_drag_K.H>


AMREX_GPU_HOST_DEVICE
amrex::Real
ComputeDragUser::operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                             amrex::Real DPM, amrex::Real, amrex::Real,
                             amrex::Real, amrex::Real, amrex::Real,
                             int, int, int, int) const
{
    amrex::Real ROg = ROPg / EPg;
    amrex::Real RE = (Mug > 0.0) ? DPM*vrel*ROg/Mug : large_number;

    amrex::Real Cd = 0.0;
    if (RE > eps) Cd = (24.0/RE)*(1.0 + 0.15*std::pow(RE, 0.687));

    return 0.75*(ROg*vrel/DPM)*Cd;
}
