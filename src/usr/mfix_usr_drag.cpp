#include <AMReX_Utility.H>

#include <mfix_des_drag_K.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE
Real
ComputeDragUser::operator() (Real EPg, Real Mug, Real ROPg, Real vrel,
                             Real DPM, Real DPA, Real PHIS,
                             Real fvelx, Real fvely, Real fvelz,
                             int i, int j, int k, int pid) const
{
    amrex::Abort("The user-defined drag routine was invoked, but no user-provided drag routine was found.");
    return 0.0;
}
