#ifndef _MFIX_ALGORITHM_H_
#define _MFIX_ALGORITHM_H_

#include <AMReX_GpuUtility.H>

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

}

#endif
