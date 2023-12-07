#ifndef actual_rhs_H
#define actual_rhs_H

#include <AMReX_RealVect.H>

#include <extern_parameters.H>
#include <actual_network.H>
#include <eos.H>
#include <burn_type.H>

using namespace amrex;
using namespace network_rp;


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real Psat_H2O(const amrex::Real sat_T);

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real Diff_H2O_Air(const amrex::Real ref_T,
                         const amrex::Real ref_P);

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real cal_NRe(const amrex::Real MUg,
                    const amrex::Real ROg,
                    const amrex::Real v_rel,
                    const amrex::Real PD,
                    const amrex::Real large_number);

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real cal_NSc(const amrex::Real MUg,
                    const amrex::Real ROg,
                    const amrex::Real Diff_coeff);

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real cal_NSh(const amrex::Real Re,
                    const amrex::Real Sc);


AMREX_INLINE
void actual_rhs_init () {}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (burn_t& state, Array1D<Real, 1, INT_NEQS>& ydot)
{
  using namespace Species;
  using namespace AuxZero;

  for (int n(1); n <= INT_NEQS; ++n)
    ydot(n) = 0.;

  if (state.aux[ifluid_homogeneous_flag] == 0)
  {
    namespace NP = NetworkProperties;

    const amrex::Real p_g = 101325.0;

    constexpr amrex::Real Pi(3.14159265358979323846);
    const amrex::Real T_p(state.aux[itemperature_p]);
    constexpr amrex::Real large_number(1.0e32);

    // Here I'm assuming to store fluid velocity into the aux variables
    amrex::RealVect vel_g(state.aux[ivelocity_x],
                          state.aux[ivelocity_y],
                          state.aux[ivelocity_z]);

    amrex::Real v_rel = (vel_g).vectorLength();
    v_rel *= 1.e2;
//    printf("v_rel = %e\n", v_rel);

    constexpr amrex::Real DP = 0.001;

    amrex::Real Sa = Pi * DP * DP;
    Sa *= 1.e4;
//    printf("Sa = %e\n", Sa);

    //amrex::Real Diff_H2OinAir = Diff_H2O_Air(state.T, 10*state.aux[ipressure]);
    amrex::Real Diff_H2OinAir = Diff_H2O_Air(state.T, 10*p_g);
//    printf("Diff_H2OinAir = %e\n", Diff_H2OinAir);

    amrex::Real MUg = 1.86e-5;
    MUg *= 1.e1;
//    printf("MUg = %e\n", MUg);

    amrex::Real ROg = state.rho;
    ROg *= 1.e-3;
//    printf("ROg = %e\n", ROg);

    amrex::Real N_Re = cal_NRe(MUg, ROg, v_rel, DP, large_number);
    N_Re *= 1.e2;
//    printf("N_Re = %e\n", N_Re);

    amrex::Real N_Sc = cal_NSc(MUg, ROg, Diff_H2OinAir);
//    printf("N_Sc = %e\n", N_Sc);

    const amrex::Real N_Sh = cal_NSh(N_Re, N_Sc);
//    printf("N_Sh = %e\n", N_Sh);

    amrex::Real H2O_xfr = (N_Sh * Diff_H2OinAir) / DP;
    H2O_xfr *= 1.e-2;
//    printf("H2O_xfr = %e\n", H2O_xfr);

    amrex::Real MW_g_MIX = 1./(state.xn[Air-1]/(2.e-3*NP::aion(Air)) +
                               state.xn[Vapor-1]/(1.e-3*NP::aion(Vapor)));
//    printf("MW_MIX_g = %e\n", MW_g_MIX);

    const amrex::Real Y_H2O = state.xn[Vapor-1] * (MW_g_MIX / (1.e-3*NP::aion(Vapor)));
//    printf("Y_H2O = %e\n", Y_H2O);

    amrex::Real Cm_H2O = Psat_H2O(T_p) / (8314.56e4*state.T);
//    printf("Cm_H2O = %e\n", Cm_H2O);

    amrex::Real Cg_H2O = (Y_H2O * (10*p_g)) / (8314.56e4*state.T);
//    printf("Cg_H2O = %e\n", Cg_H2O);

    const amrex::Real Cmg_H2O = Cm_H2O - Cg_H2O;
//    printf("Cmg_H2O = %e\n", Cmg_H2O);

    amrex::Real rate = Sa * H2O_xfr * Cmg_H2O;
//    printf("Rrate = %e\n", rate);

    ydot(Air) = 0;
    ydot(Vapor) = (Cmg_H2O > 0) ? rate : 0;
    ydot(Liquid) = (Cmg_H2O > 0) ? -1*rate : 0;

  } else if (state.aux[ifluid_homogeneous_flag] == 1) {

    ydot(Air) = 0;
    ydot(Vapor) = 0;
    ydot(Liquid) = 0;

  } else {

    amrex::Abort("Error");
  }
}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_jac (burn_t& /*state*/, MatrixType& jac)
{
    using namespace Species;

    // initialize
    jac.zero();

    // all the derivatives involving energy are 0
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real Psat_H2O(const amrex::Real sat_T)
{
  const amrex::GpuArray<amrex::Real,8> coeffs({-7.419242e0, 2.97221e-1, -1.155286e-1, 8.68563e-3,
                                               1.094098e-3, -4.39993e-3, 2.520658e-3, -5.218684e-4});

  constexpr amrex::Real cp_T = 647.286;
  constexpr amrex::Real cp_p = 22.089e6;
  constexpr amrex::Real ref_T = 338.15;
  constexpr amrex::Real ref_P = 0.01;

  amrex::Real sum1(0);
  const amrex::Real aux = ref_P*(sat_T - ref_T);
  amrex::Real m(1);

  for (int i(0); i < 8; ++i) {
    sum1 += coeffs[i]*m;
    m*= aux;
  }

  const amrex::Real val = sum1*(cp_T/sat_T - 1);

  return 10*cp_p*std::exp(val);
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real Diff_H2O_Air(const amrex::Real ref_T,
                         const amrex::Real ref_P)
{
  const amrex::Real ref_Patm = ref_P / 1.01325e6;

  return (ref_T < 450) ?
    (1.87e-6/ref_Patm) * std::pow(ref_T, 2.072) :
    (27.5e-6/ref_Patm) * std::pow(ref_T, 1.632);
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real cal_NRe(const amrex::Real MUg,
                    const amrex::Real ROg,
                    const amrex::Real v_rel,
                    const amrex::Real PD,
                    const amrex::Real large_number)
{
  return (MUg > 0) ? PD * v_rel * (ROg / MUg) : large_number;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real cal_NSc(const amrex::Real MUg,
                    const amrex::Real ROg,
                    const amrex::Real Diff_coeff)
{
  return MUg / (ROg * Diff_coeff);
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real cal_NSh(const amrex::Real Re,
                    const amrex::Real Sc)
{
  return 2.0 + 0.6*(std::sqrt(Re) * std::cbrt(Sc));
}


#endif
