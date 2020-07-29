#include <AMReX_Utility.H>

#include <mfix_des_drag_K.H>


AMREX_GPU_HOST_DEVICE
amrex::Real
ComputeDragUser::operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                             amrex::Real DPM, amrex::Real DPA, amrex::Real PHIS,
                             amrex::Real fvelx, amrex::Real fvely, amrex::Real fvelz,
                             int i, int j, int k, int pid) const
{
    amrex::Abort("The user-defined drag routine was invoked, but no user-provided drag routine was found.");
    return 0.0;
}
