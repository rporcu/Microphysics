#include <REAL.H>

extern "C"
{
  void mfix_MAIN(int* flag, Real* vol_surr, Real* A_m, Real* b_m);

  void mfix_add_argument(char* fname, int* nlen);
  void mfix_get_data(int* imax, int* jmax, int* kmax);

  void mfix_set_domain(int* flag);

  void mfix_init_fvars(Real* u, Real* v, Real* w, Real* ep,
                       Real* p, Real* ro, Real* rop, int* nlen);

  void mfix_set_constprop(Real* ro_g, Real* mu_g, Real* lambda_g, int* nlen);
}
