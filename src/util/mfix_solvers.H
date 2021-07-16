#ifndef MFIX_SOLVERS_H_
#define MFIX_SOLVERS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

#include <cmath>
#include <functional>
#include <exception>
#include <string>


namespace Solvers {


class DumpedNewton {

  public:

    // Constructor
    AMREX_GPU_HOST_DEVICE
    DumpedNewton() = default;


    // TODO description
    class NonConvergingExc : public std::exception {
      public:
        NonConvergingExc(const int iter,
                         const amrex::Real update_norm,
                         const amrex::Real residue_norm)
          : m_iter(iter)
          , m_update_norm(update_norm)
          , m_residue_norm(residue_norm)
        {}

        virtual const char* what() const throw()
        {
          std::string message =
            "DumpedNewton total iterations = " + std::to_string(m_iter) + "\n" +
            "DumpedNewton update norm = " + std::to_string(m_update_norm) + "\n" +
            "DumpedNewton residue norm = " + std::to_string(m_residue_norm) + "\n";

          return message.c_str();
        }

      private:
        const int m_iter;
        const amrex::Real m_update_norm;
        const amrex::Real m_residue_norm;
    };


    // TODO description
    class DumpingFactor {
      public:
        // Default Constructor
        AMREX_GPU_HOST_DEVICE
        DumpingFactor() : m_epg_exp(0), m_vfrac_exp(0) {};

        // Constructor
        AMREX_GPU_HOST_DEVICE
        DumpingFactor(const amrex::Real epg_exp,
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
              amrex::Real dumping_factor = 1.,
              const amrex::Real /*abs_tol*/ = 1.e-8,
              const amrex::Real rel_tol = 1.e-8,
              const int max_iterations = 100)
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
        solution_new = solution - dumping_factor * (residue / partial_R(solution));

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

          printf("DumpedNewton total iterations = %d\n", iter);
          printf("DumpedNewton update norm = %e\n", std::abs(update));
          printf("DumpedNewton residue norm = %e\n", std::abs(residue));

          return max_iterations;
        }
      } while(//(std::abs(residue) > residue_rel_tol) ||
              (std::abs(update) > update_rel_tol));

      return iter;
    }


    using ResidueMF = std::function<void(const amrex::Vector<amrex::MultiFab*>&,
                                         const amrex::Vector<amrex::MultiFab*>&)>;


    using GradientMF = std::function<void(const amrex::Vector<amrex::MultiFab*>&,
                                          const amrex::Vector<amrex::MultiFab*>&)>;


    using NormMF = std::function<amrex::Real(const amrex::Vector<amrex::MultiFab*>&)>;


    // Solve for MultiFabs
    static
    void solve(const amrex::Vector<amrex::MultiFab*>& solution,
               //const amrex::Vector<amrex::MultiFab*>& ep_g,
               ResidueMF& R,
               GradientMF& partial_R,
               NormMF& norm,
               DumpingFactor lean_rate = DumpingFactor(),
               const amrex::Real /*abs_tol*/ = 1.e-8,
               const amrex::Real rel_tol = 1.e-8,
               const int max_iterations = 100);

};

} // end namespace Solvers


#endif
