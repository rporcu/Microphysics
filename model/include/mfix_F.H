#include <REAL.H>

extern "C"
{
  void mfix_MAIN();

  void mfix_add_argument(char* fname, int* nlen);
  void mfix_get_data(int* imax, int* jmax, int* kmax);

  void mfix_set_domain(int* icbc_flag, int* flag, int* flag_e,
                       int* flag_n, int* flag_t,
                       int* nlen, const Real* dx);

  void mfix_init_fvars(Real* u, Real* v, Real* w, Real* ep,
                       Real* p, Real* ro, Real* rop, int* nlen);

  void mfix_set_constprop(Real* ro_g, Real* mu_g, Real* lambda_g, int* nlen);
}
