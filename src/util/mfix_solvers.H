#ifndef MFIX_SOLVERS_H_
#define MFIX_SOLVERS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuPrint.H>

#include <cmath>
#include <functional>
#include <exception>
#include <string>

namespace MFIXSolvers {

class Newton
{
  public:
    // Constructor
    AMREX_GPU_HOST_DEVICE
    Newton (const amrex::Real abs_tol,
            const amrex::Real rel_tol,
            const int max_iterations,
            const int is_IOProc)
      : m_abs_tol(abs_tol)
      , m_rel_tol(rel_tol)
      , m_max_iterations(max_iterations)
      , m_is_IOProc(is_IOProc)
    {}

    // Solve for scalars
    template<typename ResidueFunc, typename GradientFunc>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
#if defined(AMREX_USE_DPCPP)
    SYCL_EXTERNAL
#endif
    int solve (amrex::Real& solution,
               ResidueFunc& R,
               GradientFunc& partial_R) const
    {
      amrex::Real gradient = partial_R(solution);
      amrex::Real residue  = R(solution);
      amrex::Real update(0.);

      int iter(0);
      const amrex::Real update_tol = m_rel_tol*std::abs(solution) + m_abs_tol;
      //const amrex::Real residue_tol = m_rel_tol*std::abs(residue) + m_abs_tol;

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
        if(iter > m_max_iterations) {

          if (m_is_IOProc) {
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
      } while(//(std::abs(residue) > residue_tol) ||
              (std::abs(update) > update_tol));

      return iter;
    }

  private:
    const amrex::Real m_abs_tol;
    const amrex::Real m_rel_tol;
    const int m_max_iterations;
    const int m_is_IOProc;
};

} // end namespace MFIXSolvers


#endif
