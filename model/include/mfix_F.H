#include <REAL.H>

extern "C"
{
  void mfix_add_argument(char* fname, int* nlen);

  void mfix_get_data(int* imax, int* jmax, int* kmax, int* dem_solids);

  void mfix_set_domain(int* flag);

  void mfix_MAIN(
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
                 Real* f_gds,   Real* A_m,     Real* b_m,
                 Real* drag_am, Real* drag_bm,
                 int* flag,     Real* vol_surr,
                 int* pijk,
                 int* dg_pijk,
                 int* dg_pijkprv,
                 int* iglobal_id,
                 int* particle_state,
                 int* particle_phase,
                 Real* des_radius,
                 Real* ro_sol,
                 Real* pvol,
                 Real* pmass,
                 Real* omoi,
                 Real* ppos,
                 Real* des_pos_new,
                 Real* des_vel_new,
                 Real* des_usr_var,
                 Real* omega_new,
                 Real* des_acc_old,
                 Real* rot_acc_old,
                 Real* fc,
                 Real* tow,
                 Real* wall_collision_pft);

  void mfix_time_march(
                 Real* u_g,     Real* v_g,      Real* w_g,
                 Real* u_go,    Real* v_go,     Real* w_go,
                 Real* p_g,     Real* p_go,     Real* pp_g,    
                 Real* ep_g,    Real* ep_go,   
                 Real* ro_g,    Real* ro_go,   
                 Real* rop_g,   Real* rop_go,  
                 Real* rop_gE,  Real* rop_gN,   Real* rop_gT,
                 Real* d_e,     Real* d_n,      Real* d_t,
                 Real* tau_u_g, Real* tau_v_g,  Real* tau_w_g,
                 Real* flux_gE, Real* flux_gN,  Real* flux_gT, 
                 Real* trD_g,   Real* lambda_g, Real* mu_g,
                 Real* f_gds,   Real* A_m,      Real* b_m,
                 Real* drag_am, Real* drag_bm,
                 int* flag,     Real* vol_surr,
                 int* pijk,
                 int* dg_pijk,
                 int* dg_pijkprv,
                 int* iglobal_id,
                 int* particle_state,
                 int* particle_phase,
                 Real* des_radius,
                 Real* ro_sol,
                 Real* pvol,
                 Real* pmass,
                 Real* omoi,
                 Real* ppos,
                 Real* des_pos_new,
                 Real* des_vel_new,
                 Real* des_usr_var,
                 Real* omega_new,
                 Real* des_acc_old,
                 Real* rot_acc_old,
                 Real* fc,
                 Real* tow,
                 Real* wall_collision_pft);

  void mfix_usr3(Real* u_g,     Real* v_g,
                 Real* w_g,     Real* p_g);

}
