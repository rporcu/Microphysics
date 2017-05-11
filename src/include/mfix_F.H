#include <AMReX_REAL.H>

using namespace amrex;

extern "C"
{
    void mfix_add_argument(char* fname, int* nlen);

    void set_resid_p(Real& val);

    void write_out0(Real* time, Real* dt, Real* dx, Real* dy, Real* dz,
		    Real* xlength, Real* ylength, Real* zlength,
		    const int* domlo, const int* domhi);

    void write_usr0();

    void set_ps(Real* dx, Real* dy, Real* dz, int* err, int* is_ioproc);

    void mfix_get_data( int* fluid, int* dem,
			int* steady_state, Real* dt, Real* dt_min, Real* dt_max,
			Real* tstop, int* max_nit,
			Real* normg, int* set_normg, int* call_udf,
			int* cyclic_x, int* cyclic_y, int* cyclic_z, int* cyclic_mf,
			Real* xlength, Real* ylength, Real* zlength, int* coord);

    void check_inputs(Real* dt);

    void check_domain(Real* dx, Real* dy, Real* dz,
		      Real* xlength, Real* ylength, Real* zlength,
		      const int* domlo, const int* domhi);

    void   set_bc_flow(Real* xlength, Real* ylength, Real* zlength,
		       Real* dx,      Real* dy,      Real* dz  );
    void check_bc_flow();

    void set_bc_type(const int* slo, const int* shi,
		     int* bc_ilo, int* bc_ihi,
		     int* bc_jlo, int* bc_jhi,
		     int* bc_klo, int* bc_khi,
		     const int* domlo , const int* domhi,
		     Real* dx, Real* dy, Real* dz,
		     Real* xlen, Real* ylen, Real* zlen);

    void set_bc0(const int* slo, const int* shi,
		 const int* ulo, const int* uhi,
		 const int* vlo, const int* vhi,
		 const int* wlo, const int* whi,
		 Real* u_g,     Real* v_g,      Real* w_g,
		 Real* p_g,     Real* ep_g,
		 int* bc_ilo, int* bc_ihi,
		 int* bc_jlo, int* bc_jhi,
		 int* bc_klo, int* bc_khi,
		 const int* domlo , const int* domhi);

    void init_fluid(const int* slo, const int* shi,
		    const int* ulo, const int* uhi,
		    const int* vlo, const int* vhi,
		    const int* wlo, const int* whi,
		    const int*  lo, const int*  hi,
		    const int* domlo , const int* domhi,
		    Real* ep_g,    Real* ro_g,
		    Real* rop_g,   Real* p_g,
		    Real* u_g,     Real* v_g,      Real* w_g,
		    Real* mu_g,    Real* lambda_g,
		    Real* dx,       Real* dy,       Real* dz,
		    Real* xlength, Real* ylength, Real* zlength);

    void mfix_set_particle_properties(const int* pstate, Real* pradius, Real* pdensity,
				      Real* pvol, Real* pmass, Real* omoi);

    void mfix_init_collision();

   
   
    void mfix_output_manager ( const int* np,  Real* time, Real* dt,
			       Real* xlength, Real* ylength, Real* zlength,
			       int* nstep, void* particles, int* finsh);


    void calc_solids_volume ( const int* slo, const int* shi, const int* np,
			      void* particles, Real* dx, Real* dy, Real* dz,
			      Real* ep_g );

    void calc_volume_fraction(
	const int* slo, const int* shi,
	const int*  lo, const int*  hi,
	Real* rop_g, Real* ep_g, Real* ro_g,
	Real* dx,    Real* dy,   Real* dz);
    



    void calc_drag_fluid(
	const int* slo, const int* shi,
	const int* ulo, const int* uhi,
	const int* vlo, const int* vhi,
	const int* wlo, const int* whi,
	const int* np,
	Real* ep_g, Real* ro_g,
	Real* u_g,  Real* v_g,   Real* w_g,
	Real* mu_g, Real* f_gds, Real* drag_bm,
	void* particles,
	Real* dx,       Real* dy,       Real* dz);

    
    void calc_drag_particle(
	const int* slo, const int* shi,
	const int* ulo, const int* uhi,
	const int* vlo, const int* vhi,
	const int* wlo, const int* whi,
	const int* np,
	Real* p_g, Real* u_g,  Real* v_g,   Real* w_g,
	void* particles, 
	Real* dx, Real* dy, Real* dz, Real* xlen, Real* ylen, Real* zlen);


    void calc_coeff(
	const int* slo, const int* shi, const int* lo,  const int* hi, int* level,
	Real* ro_g, Real* p_g, Real* ep_g, Real* rop_g);

    void calc_trd_and_tau(
	const int* slo, const int* shi,
	const int* ulo, const int* uhi,
	const int* vlo, const int* vhi,
	const int* wlo, const int* whi,
	const int* lo,  const int* hi,
	Real* tau_u_g,  Real* tau_v_g, Real* tau_w_g,
	Real* trd_g,
	Real* u_g,      Real* v_g,     Real* w_g,
	Real* lambda_g, Real* mu_g,
	Real* dx,       Real* dy,       Real* dz);

    void conv_rop(
	const int* slo,       const int* shi,
	const int* ulo, const int* uhi,
	const int* vlo, const int* vhi,
	const int* wlo, const int* whi,
	const int*  lo,       const int*  hi,
	Real* u_g,            Real* v_g,            Real* w_g,
	Real* rop_g,
	Real* ropX,         Real* ropY,         Real* ropZ);

    void calc_mflux(
	const int* ulo, const int* uhi,
	const int* vlo, const int* vhi,
	const int* wlo, const int* whi,
	Real* u_g,            Real* v_g,            Real* w_g,
	Real* ropX,         Real* ropY,         Real* ropZ,
	Real* fluxX,        Real* fluxY,        Real* fluxZ,
	Real* dx,             Real* dy,             Real* dz);

    void physical_prop(const int* slo, const int* shi,
		       const int* lo, const int* hi,
		       int* level, Real* ro_g, Real* p_g,
		       Real* ep_g, Real* rop_g);

    void fill_bc0(Real* s, const int* slo, const int* shi,
		  int* bc_ilo, int* bc_ihi,
		  int* bc_jlo, int* bc_jhi,
		  int* bc_klo, int* bc_khi,
		  const int* domlo , const int* domhi);

    void solve_u_g_star(
	const int* slo, const int* shi,
	const int* ulo, const int* uhi,
	const int* vlo, const int* vhi,
	const int* wlo, const int* whi,
	const int* alo, const int* ahi,
	const int*  lo, const int*  hi,
	Real* u_g,      Real* v_g,      Real* w_g,
	Real* u_go,     Real* p_g,      Real* ro_g,
	Real* rop_g,    Real* rop_go,   Real* ep_g,
	Real* tau_u_g,  Real* d_e,
	Real* fluxX,  Real* fluxY,  Real* fluxZ,
	Real* mu_g,     Real* f_gds,
	Real* A_m,      Real* b_m,      Real* drag_bm,
	int* bc_ilo,    int* bc_ihi,
	int* bc_jlo,    int* bc_jhi,
	int* bc_klo,    int* bc_khi,
	const int* domlo , const int* domhi,
	Real* dt,
	Real* dx,       Real* dy,       Real* dz, Real (&residuals)[16]);

    void solve_v_g_star(
	const int* slo, const int* shi,
	const int* ulo, const int* uhi,
	const int* vlo, const int* vhi,
	const int* wlo, const int* whi,
	const int* alo, const int* ahi,
	const int* lo,  const int* hi,
	Real* u_g,      Real* v_g,      Real* w_g,
	Real* v_go,     Real* p_g,      Real* ro_g,
	Real* rop_g,    Real* rop_go,   Real* ep_g,
	Real* tau_v_g,  Real* d_n,
	Real* fluxX,  Real* fluxY,  Real* fluxZ,
	Real* mu_g,     Real* f_gds,
	Real* A_m,      Real* b_m,      Real* drag_bm,
	int* bc_ilo,    int* bc_ihi,
	int* bc_jlo,    int* bc_jhi,
	int* bc_klo,    int* bc_khi,
	const int* domlo , const int* domhi,
	Real* dt,
	Real* dx,       Real* dy,       Real* dz, Real (&residuals)[16]);

    void solve_w_g_star(
	const int* slo, const int* shi,
	const int* ulo,  const int* uhi,
	const int* vlo,  const int* vhi,
	const int* wlo,  const int* whi,
	const int* alo, const int* ahi,
	const int* lo,  const int* hi,
	Real* u_g,      Real* v_g,      Real* w_g,
	Real* w_go,     Real* p_g,      Real* ro_g,
	Real* rop_g,    Real* rop_go,   Real* ep_g,
	Real* tau_w_g,  Real* d_t,
	Real* fluxX,  Real* fluxY,  Real* fluxZ,
	Real* mu_g,     Real* f_gds,
	Real* A_m,      Real* b_m,      Real* drag_bm,
	int* bc_ilo,    int* bc_ihi,
	int* bc_jlo,    int* bc_jhi,
	int* bc_klo,    int* bc_khi,
	const int* domlo , const int* domhi,
	Real* dt,
	Real* dx,       Real* dy,       Real* dz, Real (&residuals)[16]);

    void solve_pp_g(
	const int* slo, const int* shi,
	const int* ulo,  const int* uhi,
	const int* vlo,  const int* vhi,
	const int* wlo,  const int* whi,
	const int* alo,  const int* ahi,
	const int* lo,  const int* hi,
	Real* u_g,      Real* v_g,      Real* w_g,
	Real* p_g,      Real* ep_g,
	Real* rop_g,    Real* rop_go,
	Real* ro_g,
	Real* ropX,   Real* ropY,   Real* ropZ,
	Real* d_e,      Real* d_n,      Real* d_t,
	Real* A_m,      Real* b_m,      Real* b_mmax,
	Real* dt,       Real* dx, Real* dy, Real* dz,
	const int* domlo , const int* domhi,
	Real (&residuals)[16]);

    void correct_0(
	const int* slo, const int* shi,
	const int* ulo,  const int* uhi,
	const int* vlo,  const int* vhi,
	const int* wlo,  const int* whi,
	const int* lo,  const int* hi,
	const int* domlo,  const int* domhi,
	Real* p_g,      Real* pp_g,
	Real* u_g,      Real* v_g,      Real* w_g,
	Real* d_e,      Real* d_n,      Real* d_t);

    void display_resid(int* nit, Real (&residuals)[16]);

    int check_convergence(int* nit, Real (&residuals)[16]);

    int goal_seek_mflux(
	const int* slo, const int* shi,
	int *nit,       int* gsmf,      Real* delP_MF,  Real* lMFlux,
	Real* fluxX,  Real* fluxY,  Real* fluxZ,
	Real* dx     ,  Real* dy,       Real* dz);

    int mfix_adjustdt(int* ier, int* nit, Real* dt);


    void mfix_des_set_dt( Real* time, Real* dt, int * nsteps, Real* dts, Real* dts_tmp);

    void mfix_call_usr3_des( const int* np, void* particles );

    void mfix_call_usr2_des( const int* np, void* particles );
    
    void mfix_des_check_dt( int* quit, Real* dt, Real* time, Real* tmp_dts );

    void mfix_des_init_time_loop ( const  Real* tstart, const Real* dt, int* nsubsteps, Real* subdt );
  
    void mfix_des_finalize_time_loop( const Real* dt, Real* dts, Real* dts_tmp);

    void mfix_des_time_loop_ops( const int* np,  void* particles, Real* subdt, Real* time, Real* dx,
				 Real* dy, Real* dz, Real* xlength, Real* ylength, Real* zlength,
				 int * nstep );

    int mfix_des_continuum_coupled ();
    
    void mfix_usr0();

    void mfix_usr1();

    void mfix_usr2();

    void mfix_usr3(Real* u_g, const int*, const int*,
		   Real* v_g, const int*, const int*,
		   Real* w_g, const int*, const int*,
		   Real* p_g, const int*, const int*,
		   Real* dx ,     Real* dy,   Real* dz);

    void mfix_finl_err_msg();

    void mfix_solve_lin_eq(int* vno, Real* var, Real* A_m, Real* b_m,
			   int* sweep_type, Real* tol, int* precond_type,
			   int* max_it, int* ier,
			   int* slo, int* shi, const int* lo, const int* hi);

    void leq_residual(const Real*, const int*, const int*,
		      const Real*, const int*, const int*,
		      const Real*, const int*, const int*,
		      const Real*, const int*, const int*);

    void leq_matvec(const Real*, const int*, const int*,
		    const Real*, const int*, const int*,
		    const Real*, const int*, const int*);

    void leq_scale(const Real*, const int*, const int*,
		   const Real*, const int*, const int*);

    void leq_msolve1(const Real*, const int*, const int*,
		     const Real*, const int*, const int*,
		     const Real*, const int*, const int*);

    void get_solver_params (int* eq_id, int* sweep, int* pc, int* maxiter, Real* tol);
}
