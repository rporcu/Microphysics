#include <REAL.H>

extern "C"
{
  void mfix_add_argument(char* fname, int* nlen);

  void mfix_get_data(int* imax, int* jmax, int* kmax, int* fluid, int* dem,
                     int* steady_state, Real* dt, Real* dt_min, Real* dt_max,
                     Real* tstop, Real* time, int* max_nit,
                     Real* normg, int* set_normg);

  void mfix_set_domain(int* flag);

  void mfix_MAIN(
                 Real* time,    Real* dt,       int* nstep,
                 Real* u_g,     Real* v_g,      Real* w_g,
                 Real* u_go,    Real* v_go,     Real* w_go,
                 Real* p_g,     Real* p_go,     Real* pp_g,
                 Real* ep_g,    Real* ep_go,
                 Real* ro_g,    Real* ro_go,
                 Real* rog_g,   Real* rog_go,
                 Real* rop_gE,  Real* rop_gN,   Real* rop_gT,
                 Real* d_e,     Real* d_n,      Real* d_t,
                 Real* tau_u_g, Real* tau_v_g,  Real* tau_w_g,
                 Real* flux_gE, Real* flux_gN,  Real* flux_gT,
                 Real* trD_g,   Real* lambda_g, Real* mu_g,
                 Real* f_gds,   Real* A_m,      Real* b_m,
                 Real* drag_bm,
                 int* flag,
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
                 Real* des_acc_old,
                 Real* rot_acc_old,
                 Real* drag_fc,
                 Real* fc,
                 Real* tow);

  void set_bc1(
         Real* time,     Real* dt,
         Real* p_g,      Real* ep_g,
         Real* ro_g,     Real* rop_g,
         Real* u_g,      Real* v_g,      Real* w_g,
         Real* flux_gE,  Real* flux_gN,  Real* flux_gT,
         int* flag);

  void mfix_output_manager(
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

  void calc_coeff_all(
         Real* ro_g, Real* p_g,
         Real* ep_g, Real* rop_g,
         Real* u_g,  Real* v_g,   Real* w_g,
         Real* mu_g, Real* f_gds, Real* drag_bm,
         int* particle_phase,  int* particle_state,
         Real* pvol, Real* des_pos_new,
         Real* des_vel_new, Real* des_radius,
         int* flag);

  void calc_coeff_all(
         Real* ro_g, Real* p_g,
         Real* ep_g, Real* rop_g,
         Real* u_g,  Real* v_g,   Real* w_g,
         Real* mu_g, Real* f_gds, Real* drag_bm,
         int* particle_phase,  int* particle_state,
         Real* pvol, Real* des_pos_new,
         Real* des_vel_new, Real* des_radius,
         int* flag);

  void calc_coeff(
         int* flag, int* level,
         Real* ro_g, Real* p_g,
         Real* ep_g, Real* rop_g,
         Real* u_g,  Real* v_g,   Real* w_g,
         Real* mu_g, Real* f_gds, Real* drag_bm,
         int* particle_phase,  int* particle_state,
         Real* pvol, Real* des_pos_new,
         Real* des_vel_new, Real* des_radius);

  void calc_trd_and_tau(
         Real* tau_u_g,  Real* tau_v_g, Real* tau_w_g,
         Real* trd_g,    Real* ep_g,
         Real* u_g,      Real* v_g,     Real* w_g,
         Real* lambda_g, Real* mu_g,
         int* flag);

  void conv_rop(const int* lo, const int* hi,
         Real* u_g,            Real* v_g,            Real* w_g,
         Real* rop_g,
         Real* rop_ge,         Real* rop_gn,         Real* rop_gt,
         int* flag,            Real* dt);

  void calc_mflux(
         Real* u_g,            Real* v_g,            Real* w_g,
         Real* rop_ge,         Real* rop_gn,         Real* rop_gt,
         Real* flux_ge,        Real* flux_gn,        Real* flux_gt,
         int* flag);

  void calc_coeff(
         int* flag, int* level,
         Real* ro_g, Real* p_g,
         Real* ep_g, Real* rop_g,
         Real* u_g,  Real* v_g,   Real* w_g,
         Real* mu_g, Real* f_gds, Real* drag_bm,
         int* particle_phase,  int* particle_state,
         Real* pvol, Real* des_pos_new,
         Real* des_vel_new, Real* des_radius);

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


  int mfix_adjustdt(int* ier, int* nit, Real* dt);

  void  mfix_des_time_march(
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
          Real* time,           Real* dt,             int* nstep);

  void mfix_usr1();

  void mfix_usr2();

  void mfix_usr3(Real* u_g,     Real* v_g,
                 Real* w_g,     Real* p_g);

}
