#ifndef LINPACK_H
#define LINPACK_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <ArrayUtilities.H>

namespace linpack {

void printA (RArray2D& a)
{
  for (int i(1); i <= 3; ++i) {
    for (int j(1); j <= 3; ++j) {
      if (j != 3) printf("%e,\t", a(i,j));
      else printf("%e\n", a(i,j));
    }
    if (i == 3) printf("\n");
  }
}

void printB (RArray1D& b)
{
  for (int j(1); j <= 3; ++j) {
    if (j != 3) printf("%e,\t", b(j));
    else printf("%e\n", b(j));
  }
  printf("\n");
}

void printP (IArray1D& b)
{
  for (int j(1); j <= 3; ++j) {
    if (j != 3) printf("%d,\t", b(j));
    else printf("%d\n", b(j));
  }
  printf("\n");
}

} // end namespace linpack

template <int num_eqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgesl (RArray2D& a, IArray1D& pivot, RArray1D& b)
{
    //printf("<<<<<<<<<<<<<<<<<<<<<<\n");

    printf("1)\n");
    linpack::printA(a);
    linpack::printB(b);
    linpack::printP(pivot);

    int nm1 = num_eqs - 1;
    printf("nm1 = %d\n", nm1);

    // solve a * x = b
    // first solve l * y = b
    if (nm1 >= 1) {
        for (int k = 1; k <= nm1; ++k) {
            printf("k = %d\n", k);
            int l = pivot(k);
            printf("l = %d\n", l);
            Real t = b(l);
            printf("t = %e\n", t);
            if (l != k) {
                b(l) = b(k);
                b(k) = t;
            }

            for (int j = k+1; j <= num_eqs; ++j) {
                printf("j = %d\n", j);
                b(j) += t * a(j,k);
            }
        }
    }

    printf("2)\n");
    linpack::printA(a);
    linpack::printB(b);
    linpack::printP(pivot);

    // now solve u * x = y
    for (int kb = 1; kb <= num_eqs; ++kb) {
        //printf("kb = %d\n", kb);

        int k = num_eqs + 1 - kb;
        //printf("k = %d\n", k);
        b(k) = b(k) / a(k,k);
        //printf("b = %e\n", b(k));
        Real t = -b(k);
        //printf("t = %e\n", t);
        for (int j = 1; j <= k-1; ++j) {
            //printf("j = %d\n", j);
            b(j) += t * a(j,k);
        }
    }

    printf("3)\n");
    linpack::printA(a);
    linpack::printB(b);
    linpack::printP(pivot);

}



template <int num_eqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgefa (RArray2D& a, IArray1D& pivot, int& info)
{
    //printf(">>>>>>>>>>>>>>>>>>>>>>\n");

    // dgefa factors a matrix by gaussian elimination.
    // a is returned in the form a = l * u where
    // l is a product of permutation and unit lower
    // triangular matrices and u is upper triangular.

    // gaussian elimination with partial pivoting

    printf("A)\n");
    linpack::printA(a);
    linpack::printP(pivot);

    info = 0;
    int nm1 = num_eqs - 1;
    //printf("nm1 = %d\n", nm1);

    Real t;

    if (nm1 >= 1) {

        for (int k = 1; k <= nm1; ++k) {
            printf("k = %d\n", k);

            // find l = pivot index
            int l = k;
            Real dmax = std::abs(a(k,k));
            for (int i = k+1; i <= num_eqs; ++i) {
                printf("i = %d\n", i);
                if (std::abs(a(i,k)) > dmax) {
                    l = i;
                    printf("l = %d\n", l);
                    dmax = std::abs(a(i,k));
                }
            }

            pivot(k) = l;

            // zero pivot implies this column already triangularized
            if (a(l,k) != 0.0e0_rt) {

                //printf("a = %e\n", a(l,k));

                // interchange if necessary
                if (l != k) {
                    t = a(l,k);
                    //printf("t = %e\n", t);
                    a(l,k) = a(k,k);
                    a(k,k) = t;
                }

                // compute multipliers
                t = -1.0e0_rt / a(k,k);
                //printf("t = %e\n", t);
                for (int j = k+1; j <= num_eqs; ++j) {
                    a(j,k) *= t;
                }

                // row elimination with column indexing
                for (int j = k+1; j <= num_eqs; ++j) {
                    //printf("j = %d\n", j);
                    t = a(l,j);
                    if (l != k) {
                        a(l,j) = a(k,j);
                        a(k,j) = t;
                        //printf("t = %e\n", t);
                    }
                    for (int i = k+1; i <= num_eqs; ++i) {
                        //printf("i = %d\n", i);
                        a(i,j) += t * a(i,k);
                    }
                }
            }
            else {

                info = k;
                //printf("k = %d\n", k);

            }

        }

    }

    printf("B)\n");
    linpack::printA(a);
    linpack::printP(pivot);

    pivot(num_eqs) = num_eqs;
    //printf("pivot = %d\n", pivot(num_eqs));

    if (a(num_eqs,num_eqs) == 0.0e0_rt) {
        info = num_eqs;
        //printf("a = %e\n", a(num_eqs,num_eqs));

    }

    printf("C)\n");
    linpack::printA(a);
    linpack::printP(pivot);

}

#endif
