#include <REAL.H>

extern "C"
{
  void mfix_add_argument(char* fname, int* nlen);

  void mfix_get_data(int* imax, int* jmax, int* kmax,
                     int* fluid, int* dem,
                     int* steady_state, Real* dt, Real* dt_min, Real* dt_max,
                     Real* tstop, Real* time, int* max_nit,
                     Real* normg, int* set_normg, int* call_udf,
                     int* cyclic_x, int* cyclic_y, int* cyclic_z, int* cyclic_mf,
                     Real* xlength, Real* ylength, Real* zlength, int* coord);

  void set_domain(int* slo, int* shi, const int* lo, const int* hi,
                  int* flag, Real* dx, Real* dy, Real* dz);

  void mfix_main1(int* slo, int* shi, const int* lo, const int* hi,
                 Real* time,    Real* dt,
                 Real* u_g,     Real* v_g,      Real* w_g,
                 Real* p_g,     Real* ep_g,
                 Real* ro_g,    Real* rop_g,
                 Real* d_e,     Real* d_n,      Real* d_t,
                 Real* flux_gE, Real* flux_gN,  Real* flux_gT,
                 Real* trD_g,   Real* lambda_g, Real* mu_g,
                 int* flag,
                 Real* dx,       Real* dy,       Real* dz);

  void init_fluid(int* slo, int* shi, const int* lo, const int* hi,
                 Real* ep_g,    Real* ro_g,
                 Real* rop_g,   Real* p_g,
                 Real* u_g,     Real* v_g,      Real* w_g,
                 Real* mu_g,    Real* lambda_g,
                 int* flag,
                 Real* dx,       Real* dy,       Real* dz);

  void mfix_make_arrays_des(const int* max_pip,
                 int* particle_state,
                 int* particle_phase,
                 Real* des_radius,
                 Real* ro_sol,
                 Real* pvol,
                 Real* pmass,
                 Real* omoi,
                 Real* des_pos_new,
                 Real* des_vel_new,
                 Real* des_usr_var,
                 Real* omega_new,
                 Real* fc,
                 Real* tow);

  void mfix_write_des_data(const int* max_pip,
                 int* particle_state,
                 Real* des_radius,
                 Real* des_pos_new,
                 Real* des_vel_new,
                 Real* des_usr_var);

  void comp_mean_fields(int* slo, int* shi, const int* lo, const int* hi, const int* max_pip,
         Real* ep_g, int* particle_state, Real* des_pos_new,
         Real* pvol, int* flag,
         Real* dx,       Real* dy,       Real* dz);

  void set_bc1(
         Real* time,     Real* dt,
         int* slo,       int* shi,
         Real* p_g,      Real* ep_g,
         Real* ro_g,     Real* rop_g,
         Real* u_g,      Real* v_g,      Real* w_g,
         Real* flux_gE,  Real* flux_gN,  Real* flux_gT,
         int* flag,
         Real* dx,       Real* dy,       Real* dz);

  void mfix_output_manager(int* slo, int* shi, const int* max_pip,
         Real* time, Real* dt, int* nstep,
         Real* ep_g,    Real* p_g,
         Real* ro_g,    Real* rop_g,
         Real* u_g,     Real* v_g,
         Real* w_g,
         int* particle_state,    Real* des_radius,
         Real* ro_sol,           Real* des_pos_new,
         Real* des_vel_new,      Real* des_usr_var,
         Real* omega_new,
         int* estatus, int* finsh);

  void calc_coeff_all(int* slo, int* shi, const int* lo, const int* hi, const int* max_pip,
         Real* ro_g, Real* p_g,
         Real* ep_g, Real* rop_g,
         Real* u_g,  Real* v_g,   Real* w_g,
         Real* mu_g, Real* f_gds, Real* drag_bm,
         int* particle_phase,  int* particle_state,
         Real* pvol, Real* des_pos_new,
         Real* des_vel_new, Real* des_radius,
         int* flag,
         Real* dx,       Real* dy,       Real* dz);

  void calc_coeff(int* slo, int* shi, const int* lo, const int* hi, const int* max_pip,
         int* flag, int* level,
         Real* ro_g, Real* p_g,
         Real* ep_g, Real* rop_g,
         Real* u_g,  Real* v_g,   Real* w_g,
         Real* mu_g, Real* f_gds, Real* drag_bm,
         int* particle_phase,  int* particle_state,
         Real* pvol, Real* des_pos_new,
         Real* des_vel_new, Real* des_radius,
         Real* dx,       Real* dy,       Real* dz);

  void calc_trd_and_tau(int* slo, int* shi, const int* lo, const int* hi,
         Real* tau_u_g,  Real* tau_v_g, Real* tau_w_g,
         Real* trd_g,    Real* ep_g,
         Real* u_g,      Real* v_g,     Real* w_g,
         Real* lambda_g, Real* mu_g,
         int* flag,
         Real* dx,       Real* dy,       Real* dz);

  void conv_rop(int* slo, int* shi, const int* lo, const int* hi,
         Real* u_g,            Real* v_g,            Real* w_g,
         Real* rop_g,
         Real* rop_ge,         Real* rop_gn,         Real* rop_gt,
         int* flag,            Real* dt,
         Real* dx,             Real* dy,             Real* dz);

  void calc_mflux(int* slo, int* shi, const int* lo, const int* hi,
         Real* u_g,            Real* v_g,            Real* w_g,
         Real* rop_ge,         Real* rop_gn,         Real* rop_gt,
         Real* flux_ge,        Real* flux_gn,        Real* flux_gt,
         int* flag,
         Real* dx,             Real* dy,             Real* dz);

  void physical_prop(int* slo, int* shi,
         int* level, Real* ro_g, Real* p_g,
         Real* ep_g, Real* rop_g, int* flag);

  void solve_u_g_star(int* slo, int* shi, const int* lo, const int* hi,
          Real* u_g,      Real* v_g,      Real* w_g,
          Real* u_go,     Real* p_g,      Real* ro_g,
          Real* rop_g,    Real* rop_go,   Real* ep_g,
          Real* tau_u_g,  Real* d_e,
          Real* flux_gE,  Real* flux_gN,  Real* flux_gT,
          Real* mu_g,     Real* f_gds,
          Real* A_m,      Real* b_m,      Real* drag_bm,
          int* flag,      Real* dt,
          Real* dx,       Real* dy,       Real* dz);

  void solve_v_g_star(int* slo, int* shi, const int* lo, const int* hi,
          Real* u_g,      Real* v_g,      Real* w_g,
          Real* v_go,     Real* p_g,      Real* ro_g,
          Real* rop_g,    Real* rop_go,   Real* ep_g,
          Real* tau_v_g,  Real* d_n,
          Real* flux_gE,  Real* flux_gN,  Real* flux_gT,
          Real* mu_g,     Real* f_gds,
          Real* A_m,      Real* b_m,      Real* drag_bm,
          int* flag,      Real* dt,
          Real* dx,       Real* dy,       Real* dz);

  void solve_w_g_star(int* slo, int* shi, const int* lo, const int* hi,
          Real* u_g,      Real* v_g,      Real* w_g,
          Real* w_go,     Real* p_g,      Real* ro_g,
          Real* rop_g,    Real* rop_go,   Real* ep_g,
          Real* tau_w_g,  Real* d_t,
          Real* flux_gE,  Real* flux_gN,  Real* flux_gT,
          Real* mu_g,     Real* f_gds,
          Real* A_m,      Real* b_m,      Real* drag_bm,
          int*  flag,     Real* dt,
          Real* dx,       Real* dy,       Real* dz);

  void solve_pp_g(int* slo, int* shi, const int* lo, const int* hi,
          Real* u_g,      Real* v_g,      Real* w_g,
          Real* p_g,      Real* ep_g,
          Real* rop_g,    Real* rop_go,
          Real* ro_g,
          Real* rop_gE,   Real* rop_gN,   Real* rop_gT,
          Real* d_e,      Real* d_n,      Real* d_t,
          Real* A_m,      Real* b_m,
          int*  flag,     Real* dt,
          Real* lnormg,   Real* resg,
          Real* dx,       Real* dy,       Real* dz);

  void correct0(int* slo, int* shi, const int* lo, const int* hi,
          Real* p_g,      Real* pp_g,
          Real* u_g,      Real* v_g,      Real* w_g,
          Real* d_e,      Real* d_n,      Real* d_t,
          int*  flag);

  void set_wall_bc(int* slo, int* shi, const int* lo, const int* hi,
          Real* u_g,      Real* v_g,      Real* w_g,
          int*  flag);

  void mfix_iterate(
         Real* u_g,            Real* v_g,            Real* w_g,
         Real* u_go,           Real* v_go,           Real* w_go,
         Real* p_g,            Real* pp_g,
         Real* ep_g,           Real* ro_g,
         Real* rop_g,          Real* rop_go,
         Real* rop_ge,         Real* rop_gn,         Real* rop_gt,
         Real* d_e,            Real* d_n,            Real* d_t,
         Real* flux_ge,        Real* flux_gn,        Real* flux_gt,
         Real* mu_g,           Real* f_gds,
         Real* A_m,            Real* b_m,            Real* drag_bm,
         Real* tau_u_g,        Real* tau_v_g,        Real* tau_w_g,
         int* particle_phase,  int* particle_state,
         Real* pvol,           Real* des_radius,
         Real* des_pos_new,    Real* des_vel_new,
         int* flag,
         Real* time,           Real* dt,             int* nstep,
         int* ier,             int* nit,
         int* mustit, int* gsmf,
         Real* delP_MF, Real* lMFlux, Real* resg,
         Real* lnormg, int* lset_normg);

  void display_resid(int* nit);

  int check_convergence(int* nit);

  int goal_seek_mFlux(const int* slo, const int* shi,
              int*nit,        int* gsmf,      Real* delP_MF,  Real* lMFlux,
              Real* flux_gE,  Real* flux_gN,  Real* flux_gT,  int* flag,
              Real* dx     ,  Real* dy,       Real* dz);

  int mfix_adjustdt(int* ier, int* nit, Real* dt);

  void  mfix_des_time_march(
          int* slo, int* shi, const int* lo, const int* hi, const int* max_pip,
          Real* ep_g,           Real* p_g,
          Real* u_g,            Real* v_g,            Real* w_g,
          Real* ro_g,           Real* rop_g,          Real* mu_g,
          int* particle_state,  int* particle_phase,
          Real* des_radius,     Real* ro_sol,
          Real* pvol,           Real* pmass,
          Real* omoi,           Real* des_usr_var,
          Real* des_pos_new,    Real* des_vel_new,    Real* omega_new,
          Real* des_acc_old,    Real* rot_acc_old,
          Real* drag_fc,        Real* fc,             Real* tow,
          int* pairs,           int* pair_count,
          int* flag,
          Real* time,           Real* dt,
          Real* dx  ,           Real* dy,             Real* dz,
          int * nstep);

  void mfix_usr0();

  void mfix_usr1();

  void mfix_usr2();

  void mfix_usr3(int* slo, int* shi,
                 Real* u_g,     Real* v_g,
                 Real* w_g,     Real* p_g,
                 Real* dx ,     Real* dy,   Real* dz);

  void mfix_finl_err_msg();

  void mfix_solve_lin_eq(int* vno, Real* var, Real* A_m, Real* b_m,
                         int* sweep_type, Real* tol, int* precond_type,
                         int* max_it, int* ier,
                         int* slo, int* shi, const int* lo, const int* hi);

  void leq_residual(const Real*, const Real*, const Real*, Real*,
                    int* slo, int* shi, const int* lo, const int* hi);

  void leq_matvec(const Real*, const Real*, Real*,
                    int* slo, int* shi, const int* lo, const int* hi);

  void leq_scale(const Real*, const Real*,
                 int* slo, int* shi, const int* lo, const int* hi);

  void leq_msolve1(const int* slo, const int* shi, const Real*, const Real*, Real*,
                     const int* sweep_type);

  void get_solver_params (int* eq_id, int* sweep, int* pc, int* maxiter, Real* tol);
}
