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


namespace Solvers {


class DampedNewton {

  public:

    // Constructor
    AMREX_GPU_HOST_DEVICE
    DampedNewton() = default;


    // TODO description
    class NonConvergingExc : public std::exception {
      public:
        NonConvergingExc(const int is_IOProc,
                         const int iter,
                         const amrex::Real update_norm,
                         const amrex::Real residue_norm)
          : m_is_IOProc(is_IOProc)
          , m_iter(iter)
          , m_update_norm(update_norm)
          , m_residue_norm(residue_norm)
        {}

        virtual const char* what() const throw()
        {
          if (m_is_IOProc) {
            std::printf("Damped-Newton total iterations = %d\n", m_iter);
            std::printf("Damped-Newton update norm = %e\n", m_update_norm);
            std::printf("Damped-Newton residue norm = %e\n", m_residue_norm);
          }

          return "Damped-Newton solver did not converge";
        }

      private:
        const int m_is_IOProc;
        const int m_iter;
        const amrex::Real m_update_norm;
        const amrex::Real m_residue_norm;
    };


    // TODO description
    class DampingFactor {
      public:
        // Default Constructor
        AMREX_GPU_HOST_DEVICE
        DampingFactor() : m_epg_exp(0), m_vfrac_exp(0) {};

        // Constructor
        AMREX_GPU_HOST_DEVICE
        DampingFactor(const amrex::Real epg_exp,
                      const amrex::Real vfrac_exp)
          : m_epg_exp(epg_exp)
          , m_vfrac_exp(vfrac_exp)
        {};

        // Functor
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::Real operator() (const amrex::Real epg,
                                const amrex::Real vfrac) const
        {
          return std::pow(epg, m_epg_exp)*std::pow(vfrac, m_vfrac_exp);
        }

      private:
        amrex::Real m_epg_exp;
        amrex::Real m_vfrac_exp;
    };


    // Solve for scalars
    template<typename ResidueFunc, typename GradientFunc>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static
    int solve(amrex::Real& solution_new,
              ResidueFunc& R,
              GradientFunc& partial_R,
              const int is_IOProc,
              amrex::Real damping_factor,
              const amrex::Real /*abs_tol*/,
              const amrex::Real rel_tol,
              const int max_iterations)
    {
      amrex::Real solution_old = solution_new;
      const amrex::Real original_solution(solution_new);

      amrex::Real residue = R(solution_new);
      amrex::Real update(0.);

      int iter(0);
      const amrex::Real update_rel_tol = rel_tol*std::abs(solution_new);
      //const amrex::Real residue_rel_tol = abs_tol*std::abs(residue);

      do {
        amrex::Real solution = solution_new;
        solution_new = solution - damping_factor * (residue / partial_R(solution));

        if (iter != 0) {
          if (std::abs(solution_new - solution_old) < std::abs(solution_new - solution)) {
            solution_new = (solution + solution_old) / 2.;
          }
          else if (std::abs(solution - solution_old) < std::abs(solution_new - solution)) {
            solution_new = (solution_new + solution_old) / 2.;
          }
        }

        residue = R(solution_new);
        update = solution_new - solution;
        solution_old = solution;

        ++iter;
        if(iter > max_iterations) {
          solution_new = original_solution;

          if (is_IOProc) {
#ifdef AMREX_USE_GPU
            AMREX_DEVICE_PRINTF("Damped-Newton total iterations = %d\n", iter);
            AMREX_DEVICE_PRINTF("Damped-Newton update norm = %e\n", std::abs(update));
            AMREX_DEVICE_PRINTF("Damped-Newton residue norm = %e\n", std::abs(residue));
#else
            std::printf("Damped-Newton total iterations = %d\n", iter);
            std::printf("Damped-Newton update norm = %e\n", std::abs(update));
            std::printf("Damped-Newton residue norm = %e\n", std::abs(residue));
#endif
          }

          return max_iterations;
        }
      } while(//(std::abs(residue) > residue_rel_tol) ||
              (std::abs(update) > update_rel_tol));

      return iter;
    }


    using NormMF = std::function<amrex::Real(const amrex::Vector<amrex::MultiFab*>&)>;
};

} // end namespace Solvers


#endif
