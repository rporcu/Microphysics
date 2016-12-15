#include <REAL.H>

extern "C"
{
  void mfix_MAIN(int* flag,     Real* vol_surr,
                 Real* A_m,     Real* b_m,
                 Real* ep_g,    Real* ep_go,
                 Real* p_g,     Real* p_go,
                 Real* ro_g,    Real* ro_go,
                 Real* rog_g,   Real* rog_go,
                 Real* u_g,     Real* u_go,
                 Real* v_g,     Real* v_go,
                 Real* w_g,     Real* w_go,
                 Real* pp_g,    Real* d_e,
                 Real* d_n,     Real* d_t,
                 Real* mu_g,    Real* lambda_g,
                 Real* trD_g,   Real* tau_u_g,
                 Real* tau_v_g, Real* tau_w_g,
                 Real* flux_gE, Real* flux_gN,
                 Real* flux_gT, Real* rop_gE,
                 Real* rop_gN,  Real* rop_gT,
                 Real* f_gds,   Real* drag_am,
                 Real* drag_bm  );

  void mfix_add_argument(char* fname, int* nlen);
  void mfix_get_data(int* imax, int* jmax, int* kmax);

  void mfix_set_domain(int* flag);

  void mfix_init_fvars(Real* u, Real* v, Real* w, Real* ep,
                       Real* p, Real* ro, Real* rop, int* nlen);

  void mfix_set_constprop(Real* ro_g, Real* mu_g, Real* lambda_g, int* nlen);
}
