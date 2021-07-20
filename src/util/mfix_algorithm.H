#ifndef MFIX_ALGORITHM_H_
#define MFIX_ALGORITHM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

#include <cmath>


namespace aux {

template <class T>
class is_equal {
  public:
    AMREX_GPU_HOST_DEVICE
    is_equal (const T& value):
      m_value(value)
    {}

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    bool operator () (const T& other)
    { return other == m_value; }

  private:
    const T m_value;
};


template<class InputIterator, class UnaryPredicate>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool any_of (InputIterator first, InputIterator last, UnaryPredicate pred)
{
  while (first != last) {
    if (pred(*first))
      return true;
    
    ++first;
  }

  return false; 
}

} // end namespace aux


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int MFIXfind (const int* array, const int N, const int elem)
{
  for (int n(0); n < N; n++)
    if (array[n] == elem)
      return n;

  return -1; // -1 is the InvalidIdx value
}


namespace Solvers {

template<typename ResidueFunc, typename GradientFunc>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void NewtonRaphson(amrex::Real& solution_new,
                   const amrex::Real& initial_value,
                   ResidueFunc& R,
                   GradientFunc& partial_R,
                   const amrex::Real /*abs_tol*/ = 1.e-8,
                   const amrex::Real rel_tol = 1.e-8,
                   const int max_iterations = 100)
{
  amrex::Real solution_old = initial_value;
  amrex::Real residue = R(solution_old);
  amrex::Real update(0.);

  int iter(0);
  const amrex::Real update_rel_tol = rel_tol*std::abs(solution_old);
  //const amrex::Real residue_rel_tol = abs_tol*std::abs(residue);

  do {
    solution_new = solution_old - residue / partial_R(solution_old);
    residue = R(solution_new);
    update = std::abs(solution_new - solution_old);
    solution_old = solution_new;

    ++iter;
    if(iter > max_iterations) {
      amrex::Abort("Error: Solver did not converge");
    }
  } while(//(std::abs(residue) > residue_rel_tol) ||
          (std::abs(update) > update_rel_tol));
}

template<typename ResidueFunc, typename GradientFunc>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void NewtonStabilized(amrex::Real& solution_new,
                      const amrex::Real& initial_value,
                      ResidueFunc& R,
                      GradientFunc& partial_R,
                      const amrex::Real /*abs_tol*/ = 1.e-8,
                      const amrex::Real rel_tol = 1.e-8,
                      const int max_iterations = 100)
{
  amrex::Real solution_old = initial_value;
  solution_new = initial_value;
  amrex::Real residue = R(solution_new);
  amrex::Real update(0.);

  int iter(0);
  const amrex::Real update_rel_tol = rel_tol*std::abs(solution_new);
  //const amrex::Real residue_rel_tol = abs_tol*std::abs(residue);

  do {
    amrex::Real solution = solution_new;
    solution_new = solution - residue / partial_R(solution);

    if (iter != 0 and std::abs(solution_new - solution_old) < std::abs(solution_new - solution)) {
      if (std::abs(solution_new - solution) < std::abs(solution_old - solution)) {
        solution_new = (solution_new+solution) / 2.;
      } else {
        solution_new = (solution_old+solution) / 2.;
      }
    }

    residue = R(solution_new);
    update = solution_new - solution;
    solution_old = solution;

    ++iter;
    if(iter > max_iterations) {
      amrex::Abort("Error: Solver did not converge");
    }
  } while(//(std::abs(residue) > residue_rel_tol) ||
          (std::abs(update) > update_rel_tol));
}

} // end namespace Solvers


#endif
