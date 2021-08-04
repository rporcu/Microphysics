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


#endif
