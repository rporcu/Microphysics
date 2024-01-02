#ifndef actual_rhs_H
#define actual_rhs_H

#include <extern_parameters.H>
#include <actual_network.H>
#include <burn_type.H>

using namespace amrex;
using namespace network_rp;

AMREX_INLINE
void actual_rhs_init () {
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    using namespace Species;

    // Note: we are assuming that A = 1 for all species so molar and mass fractions are the same

    ydot(Position) = state.xn[Velocity-1];
    ydot(Velocity) = -1.0_rt*state.xn[Position-1];

    // this network has no energy generation
    ydot(net_ienuc) = 0.0_rt;

}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_jac (burn_t& state, MatrixType& jac)
{
    using namespace Species;

    // initialize
    jac.zero();

    // carbon jacobian elements
    jac(Position, Position) = 0.0_rt;
    jac(Position, Velocity) = 1.0_rt;

    jac(Velocity, Position) = -1.0_rt;
    jac(Velocity, Velocity) = 0.0_rt;

    // all the derivatives involving energy are 0
}

#endif
