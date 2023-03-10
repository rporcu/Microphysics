#ifndef MFIX_SOLVERS_H_
#define MFIX_SOLVERS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuPrint.H>

#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_run_on.H>

#include <cmath>
#include <functional>
#include <exception>
#include <string>

namespace Solvers { namespace Newton {

namespace FluidEnthalpy {

struct Residue
{
  AMREX_GPU_HOST_DEVICE
  Residue (const int& i,
           const int& j,
           const int& k,
           const int& fluid_is_a_mixture,
           const int& cell_is_covered,
           const int& nspecies_g,
           const MFIXFluidParms& fluid_parms,
           const amrex::Array4<const amrex::Real>& X_gk,
           const amrex::Real& h_g)
  : m_i(i)
  , m_j(j)
  , m_k(k)
  , m_fluid_is_a_mixture(fluid_is_a_mixture)
  , m_cell_is_covered(cell_is_covered)
  , m_nspecies_g(nspecies_g)
  , m_fluid_parms(fluid_parms)
  , m_X_gk(X_gk)
  , m_h_g(h_g)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() (const amrex::Real& Tg_arg)
  {
    amrex::Real hg_loc(0);

    if (!m_fluid_is_a_mixture) {
      hg_loc = m_fluid_parms.calc_h_g<run_on>(Tg_arg, m_cell_is_covered);
    } else {
      for (int n_g(0); n_g < m_nspecies_g; ++n_g) {
        hg_loc += m_X_gk(m_i,m_j,m_k,n_g)*m_fluid_parms.calc_h_gk<run_on>(Tg_arg, n_g, m_cell_is_covered);
      }
    }

    return hg_loc - m_h_g;
  }

  const int& m_i; const int& m_j; const int& m_k;
  const int& m_fluid_is_a_mixture;
  const int& m_cell_is_covered;
  const int& m_nspecies_g;
  const MFIXFluidParms& m_fluid_parms;
  const amrex::Array4<const amrex::Real>& m_X_gk;
  const amrex::Real& m_h_g;
};

struct Gradient
{
  AMREX_GPU_HOST_DEVICE
  Gradient (const int& i,
            const int& j,
            const int& k,
            const int& fluid_is_a_mixture,
            const int& nspecies_g,
            const MFIXFluidParms& fluid_parms,
            const amrex::Array4<const amrex::Real>& X_gk)
  : m_i(i)
  , m_j(j)
  , m_k(k)
  , m_fluid_is_a_mixture(fluid_is_a_mixture)
  , m_nspecies_g(nspecies_g)
  , m_fluid_parms(fluid_parms)
  , m_X_gk(X_gk)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() (const amrex::Real& Tg_arg)
  {
    amrex::Real gradient(0);

    if (!m_fluid_is_a_mixture) {
      gradient = m_fluid_parms.calc_partial_h_g<run_on>(Tg_arg);
    } else {
      for (int n_g(0); n_g < m_nspecies_g; ++n_g) {
        gradient += m_X_gk(m_i,m_j,m_k,n_g)*m_fluid_parms.calc_partial_h_gk<run_on>(Tg_arg,n_g);
      }
    }

    return gradient;
  }

  const int& m_i; const int& m_j; const int& m_k;
  const int& m_fluid_is_a_mixture;
  const int& m_nspecies_g;
  const MFIXFluidParms& m_fluid_parms;
  const amrex::Array4<const amrex::Real>& m_X_gk;
};

} // end namespace FluidEnthalpy

namespace SolidsEnthalpy {

struct Residue
{
  AMREX_GPU_HOST_DEVICE
  Residue (const int& solid_is_a_mixture,
           const int& nspecies_s,
           const MFIXSolidsParms& solids_parms,
           const amrex::Real* X_sn,
           const amrex::Real& hp)
    : m_solid_is_a_mixture(solid_is_a_mixture)
    , m_nspecies_s(nspecies_s)
    , m_solids_parms(solids_parms)
    , m_X_sn(X_sn)
    , m_hp(hp)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() (const amrex::Real& Tp_arg)
  {
    amrex::Real hp_loc(0);

    if (!m_solid_is_a_mixture) {
      hp_loc = m_solids_parms.calc_h_s<run_on>(Tp_arg);
    } else {
      for (int n_s(0); n_s < m_nspecies_s; ++n_s)
        hp_loc += m_X_sn[n_s]*m_solids_parms.calc_h_sn<run_on>(Tp_arg,n_s);
    }

    return hp_loc - m_hp;
  }

  const int& m_solid_is_a_mixture;
  const int& m_nspecies_s;
  const MFIXSolidsParms& m_solids_parms;
  const amrex::Real* m_X_sn;
  const amrex::Real& m_hp;
};

struct Gradient
{
  AMREX_GPU_HOST_DEVICE
  Gradient (const int& solid_is_a_mixture,
            const int& nspecies_s,
            const MFIXSolidsParms& solids_parms,
            const amrex::Real* X_sn)
    : m_solid_is_a_mixture(solid_is_a_mixture)
    , m_nspecies_s(nspecies_s)
    , m_solids_parms(solids_parms)
    , m_X_sn(X_sn)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() (const amrex::Real& Tp_arg)
  {
    amrex::Real gradient(0);

    if (!m_solid_is_a_mixture) {
      gradient = m_solids_parms.calc_partial_h_s<run_on>(Tp_arg);
    } else {
      for (int n_s(0); n_s < m_nspecies_s; ++n_s) {
        gradient += m_X_sn[n_s]*m_solids_parms.calc_partial_h_sn<run_on>(Tp_arg,n_s);
      }
    }

    return gradient;
  }

  const int& m_solid_is_a_mixture;
  const int& m_nspecies_s;
  const MFIXSolidsParms& m_solids_parms;
  const amrex::Real* m_X_sn;
};

} // end namespace SolidsEnthalpy

// Solve for scalars
template<typename ResidueFunc, typename GradientFunc>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int solve (amrex::Real& solution,
           ResidueFunc& R,
           GradientFunc& partial_R,
           const amrex::Real /*abs_tol*/,
           const amrex::Real rel_tol,
           const int max_iterations,
           const int is_IOProc = 0)
{
  amrex::Real gradient = partial_R(solution);
  amrex::Real residue  = R(solution);
  amrex::Real update(0.);

  int iter(0);
  const amrex::Real update_rel_tol = rel_tol*std::abs(solution);
  //const amrex::Real residue_rel_tol = rel_tol*std::abs(residue);

  // Do not do anything if we're already on the solution
  if (std::abs(residue) < 1.e-15) {
    return iter;
  }

  if (std::abs(gradient) < 1.e-15) {
    amrex::Abort("Gradient in Newton solver is almost zero");
  }

  do {
    update = residue / gradient;
    solution -= update;
    residue = R(solution);
    gradient = partial_R(solution);

    ++iter;
    if(iter > max_iterations) {

      if (is_IOProc) {
#ifdef AMREX_USE_GPU
        AMREX_DEVICE_PRINTF("Newton total iterations = %d\n", iter);
        AMREX_DEVICE_PRINTF("Newton update norm = %e\n", std::abs(update));
        AMREX_DEVICE_PRINTF("Newton residue norm = %e\n", std::abs(residue));
#else
        std::printf("Newton total iterations = %d\n", iter);
        std::printf("Newton update norm = %e\n", std::abs(update));
        std::printf("Newton residue norm = %e\n", std::abs(residue));
#endif
      }

      amrex::Abort("Newton solver did not converge");

    }
  } while(//(std::abs(residue) > residue_rel_tol) ||
          (std::abs(update) > update_rel_tol));

  return iter;
}

} } // end namespace Solvers::Newton


#endif
