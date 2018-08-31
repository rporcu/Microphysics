#ifndef MFIX_IF_LIST_
#define MFIX_IF_LIST_

#include <AMReX_Vector.H>
#include <AMReX_EB2.H>

#include <algorithm>

using namespace amrex;

template <class F>
class UnionListIF {

public:

    UnionListIF (const Vector<F> a_ifs)
        : m_ifs(a_ifs)
    {}

    ~UnionListIF () {}

    UnionListIF (const UnionListIF & rhs) = default;
    UnionListIF (UnionListIF && rhs)      = default;
    UnionListIF & operator= (UnionListIF & rhs)  = default;
    UnionListIF & operator= (UnionListIF && rhs) = default;

    Real operator() (const RealArray & p) const {
        F & f_max = std::max( m_ifs,
                      [&](const F & f1, const F & f2) {
                          return f1(p) < f2(p);
                      });
        return f_max(p);
    }

private:

    Vector<F> m_ifs;
};

#endif
