#ifndef VODE_DVJAC_H
#define VODE_DVJAC_H

#include <vode_type.H>
#ifndef NEW_NETWORK_IMPLEMENTATION
#include <linpack.H>
#endif
#ifdef STRANG
#include <integrator_rhs_strang.H>
#endif
#ifdef SDC
#include <integrator_rhs_simplified_sdc.H>
#endif

namespace dvjac_aux {

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

} // end namespace dvjac_aux

template <typename IArray1D, typename BurnT, typename DvodeT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dvjac (IArray1D& pivot, int& IERPJ, BurnT& state, DvodeT& vstate)
{
    // dvjac is called by dvnlsd to compute and process the matrix
    // P = I - h*rl1*J , where J is an approximation to the Jacobian
    // that we obtain either through direct evaluation or caching from
    // a previous evaluation. P is then subjected to LU decomposition
    // in preparation for later solution of linear systems with P as
    // coefficient matrix. This is done by DGEFA.

    constexpr int int_neqs = integrator_neqs<BurnT>();

    IERPJ = 0;

#ifdef ALLOW_JACOBIAN_CACHING
    // See whether the Jacobian should be evaluated. Start by basing
    // the decision on whether we're caching the Jacobian.

    int evaluate_jacobian = 1;

    if (use_jacobian_caching) {
        evaluate_jacobian = 0;
    }

    if (use_jacobian_caching) {
        // Now evaluate the cases where we're caching the Jacobian but aren't
        // going to be using the cached Jacobian.

        // On the first step we don't have a cached Jacobian. Also, after enough
        // steps, we consider the cached Jacobian too old and will want to re-evaluate
        // it, so we look at whether the step of the last Jacobian evaluation (NSLJ)
        // is more than max_steps_between_jacobian_evals steps in the past.
        if (vstate.NST == 0 || vstate.NST > vstate.NSLJ + max_steps_between_jacobian_evals) {
            evaluate_jacobian = 1;
        }

        // See the non-linear solver for details on these conditions.
        if (vstate.ICF == 1 && vstate.DRC < CCMXJ) {
            evaluate_jacobian = 1;
        }

        if (vstate.ICF == 2) {
            evaluate_jacobian = 1;
        }

    }

    if (evaluate_jacobian == 1) {
#endif

        // We want to evaluate the Jacobian -- now the path depends on
        // whether we're using the numerical or analytic Jacobian.

        if (vstate.jacobian_type == 1) {

            // For the analytic Jacobian, call the user-supplied function.

            // Increment the Jacobian evaluation counter.
            vstate.NJE += 1;

            // Refresh the timestep marker for the last Jacobian evaluation.
            vstate.NSLJ = vstate.NST;

            // Indicate that the Jacobian is current for this solve.
            vstate.JCUR = 1;

            // Initialize the Jacobian to zero
            vstate.jac.zero();

            printf("11)\n");
            dvjac_aux::printA(vstate.jac);

            jac(vstate.tn, state, vstate, vstate.jac);

#ifdef ALLOW_JACOBIAN_CACHING
            // Store the Jacobian if we're caching.
            if (use_jacobian_caching == 1) {
                vstate.jac_save = vstate.jac;
            }
#endif

        }
        else {

            // For the numerical Jacobian, make N calls to the RHS to approximate it.

            // Increment the Jacobian evaluation counter.
            vstate.NJE += 1;

            // Refresh the timestep marker for the last Jacobian evaluation.
            vstate.NSLJ = vstate.NST;

            // Indicate that the Jacobian is current for this solve.
            vstate.JCUR = 1;

            Real fac = 0.0_rt;
            for (int i = 1; i <= int_neqs; ++i) {
                fac += (vstate.savf(i) * vstate.ewt(i)) * (vstate.savf(i) * vstate.ewt(i));
            }
            fac = std::sqrt(fac / int_neqs);

            Real R0 = 1000.0_rt * std::abs(vstate.H) * UROUND * int_neqs * fac;
            if (R0 == 0.0_rt) {
                R0 = 1.0_rt;
            }

            bool in_jacobian = true;
            for (int j = 1; j <= int_neqs; ++j) {
                Real yj = vstate.y(j);
                printf("yj = %e\n", yj);

                Real R = amrex::max(std::sqrt(UROUND) * std::abs(yj), R0 / vstate.ewt(j));
                printf("R = %e\n", R);
                vstate.y(j) += R;
                fac = 1.0_rt / R;
                printf("fac = %e\n", fac);

                rhs(vstate.tn, state, vstate, vstate.acor, in_jacobian);
                for (int i = 1; i <= int_neqs; ++i) {
                    vstate.jac.set(i, j, (vstate.acor(i) - vstate.savf(i)) * fac);
                }

                vstate.y(j) = yj;
            }

            printf("22)\n");
            dvjac_aux::printA(vstate.jac);

            // Increment the RHS evaluation counter by N.
            vstate.NFE += int_neqs;

#ifdef ALLOW_JACOBIAN_CACHING
            // Store the Jacobian if we're caching.
            if (use_jacobian_caching == 1) {
                vstate.jac_save = vstate.jac;
            }
#endif

        }

#ifdef ALLOW_JACOBIAN_CACHING
    }
    else {

        // Load the cached Jacobian.

        // Indicate the Jacobian is not current for this step.
        vstate.JCUR = 0;
        vstate.jac = vstate.jac_save;

        printf("33)\n");
        dvjac_aux::printA(vstate.jac);

    }
#endif

    // Multiply Jacobian by a scalar, add the identity matrix
    // (along the diagonal), and do LU decomposition.

    Real hrl1 = vstate.H * vstate.RL1;
    printf("vstate.H = %e\n", vstate.H);
    printf("vstate.RL1 = %e\n", vstate.RL1);
    printf("hrl1 = %e\n\n", hrl1);
    Real con = -hrl1;

    vstate.jac.mul(con);

    printf("44)\n");
    dvjac_aux::printA(vstate.jac);

    vstate.jac.add_identity();

    printf("55)\n");
    dvjac_aux::printA(vstate.jac);

    int IER;

#ifdef NEW_NETWORK_IMPLEMENTATION
    amrex::ignore_unused(pivot);
    RHS::dgefa(vstate.jac);
    IER = 0;
#else
    dgefa<int_neqs>(vstate.jac, pivot, IER);
#endif

    if (IER != 0) {
        IERPJ = 1;
    }
}

#endif
