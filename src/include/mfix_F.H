#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include "AMReX_RealVect.H"

#ifdef BL_SINGLE_PRECISION_PARTICLES
typedef float amrex_particle_real;
#else
typedef double amrex_particle_real;
#endif

#ifdef __cplusplus
extern "C"
{
#endif
    void mfix_add_argument(char* fname, int* nlen);

    void set_resid_p(amrex::Real& val);

    void write_out0(amrex::Real* time, amrex::Real* dt, amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
        amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
        const int* domlo, const int* domhi);

    void write_usr0();

    void set_ps(amrex::Real* dx, amrex::Real* dy, amrex::Real* dz, int* err, int* is_ioproc);

    void mfix_get_data( int* fluid, int* dem,
      int* steady_state, amrex::Real* dt, amrex::Real* dt_min, amrex::Real* dt_max,
      amrex::Real* tstop, int* max_nit,
      amrex::Real* normg, int* set_normg, int* call_udf);

    void check_inputs(amrex::Real* dt);

    void mfix_set_cyclic(int* cyc_x, int* cyc_y, int* cyc_z);

    void check_initial_conditions(amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
          const int* domlo, const int* domhi);

    void check_boundary_conditions(amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
           amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
           const int* domlo, const int* domhi);

    void check_point_sources(amrex::Real* dx, amrex::Real* dy, amrex::Real* dz);

    void set_bc_flow(amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
         amrex::Real* dx,      amrex::Real* dy,      amrex::Real* dz  );

    void check_bc_flow();

    void set_bc_type(
        int* bc_ilo, int* bc_ihi,
        int* bc_jlo, int* bc_jhi,
        int* bc_klo, int* bc_khi,
        const int* domlo , const int* domhi,
        amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
        amrex::Real* xlen, amrex::Real* ylen, amrex::Real* zlen,
        const int* nghost_bc);

    void set_bc0(const int* slo, const int* shi,
     const int* ulo, const int* uhi,
     const int* vlo, const int* vhi,
     const int* wlo, const int* whi,
     amrex::Real* u_g,     amrex::Real* v_g,      amrex::Real* w_g,
     amrex::Real* p_g,     amrex::Real* ep_g,
     amrex::Real* ro_g,    amrex::Real* rop_g,
     amrex::Real* mu_g,    amrex::Real* lambda_g,
     int* bc_ilo, int* bc_ihi,
     int* bc_jlo, int* bc_jhi,
     int* bc_klo, int* bc_khi,
     const int* domlo , const int* domhi);

    void set_bc1(const int* slo, const int* shi,
     const int* ulo, const int* uhi,
     const int* vlo, const int* vhi,
     const int* wlo, const int* whi,
     amrex::Real* u_g, amrex::Real* v_g,      amrex::Real* w_g,
     amrex::Real* p_g,     amrex::Real* ep_g,
     amrex::Real* ro_g,    amrex::Real* rop_g,
     amrex::Real* mu_g,    amrex::Real* lambda_g,
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
		    amrex::Real* ep_g,    amrex::Real* ro_g,
		    amrex::Real* rop_g,   amrex::Real* p_g,      amrex::Real* p0_g,
		    amrex::Real* u_g,     amrex::Real* v_g,      amrex::Real* w_g,
		    amrex::Real* mu_g,    amrex::Real* lambda_g,
		    amrex::Real* dx,       amrex::Real* dy,       amrex::Real* dz,
		    amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength);

    
    void init_fluid_restart(const int* slo, const int* shi,
          const int*  lo, const int*  hi,
          amrex::Real* mu_g,    amrex::Real* lambda_g);

    void flip_particle_vol(const int* slo, const int* shi,
         amrex::Real* ep_g,
         int* bc_ilo, int* bc_ihi,
         int* bc_jlo, int* bc_jhi,
         int* bc_klo, int* bc_khi,
         const int* domlo , const int* domhi);

    void set_drag_bcs(
         const int* ulo, const int* uhi,
         const int* vlo, const int* vhi,
         const int* wlo, const int* whi,
         amrex::Real* drag_u, amrex::Real* drag_v, amrex::Real* drag_w,
         int* bc_ilo, int* bc_ihi,
         int* bc_jlo, int* bc_jhi,
         int* bc_klo, int* bc_khi,
         const int* domlo , const int* domhi);

    void zero_wall_norm_vel(const int* slo, const int* shi,
          const int* ulo, const int* uhi,
          const int* vlo, const int* vhi,
          const int* wlo, const int* whi,
          amrex::Real* u_g,     amrex::Real* v_g,      amrex::Real* w_g,
          int* bc_ilo, int* bc_ihi,
          int* bc_jlo, int* bc_jhi,
          int* bc_klo, int* bc_khi,
          const int* domlo , const int* domhi);

    void set_particle_properties(const int* pstate, amrex::Real* pradius, amrex::Real* pdensity,
         amrex::Real* pvol, amrex::Real* pmass, amrex::Real* omoi, amrex::Real* omega);

    void init_collision(amrex::Real *sum_dp, amrex::Real *sum_ro);

    void sum_particle_props( const int* nrp,  void* particles,
           amrex::Real *sum_np, amrex::Real *sum_dp,
           amrex::Real *sum_ro );


    void output_manager ( const int* np,  amrex::Real* time, amrex::Real* dt,
             amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
             int* nstep, void* particles, int* finsh);

    void calc_volume_fraction (
        const int*  lo, const int*  hi,
        const int* slo, const int* shi,
        const int*  np, void* particles,
        amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
        amrex::Real* ep_g );

    void mfix_deposit_cic(const amrex_particle_real*, int ns, int np, int nc,
                          amrex_real* mf, const int* lo, const int* hi,
                          const amrex_real* plo, const amrex_real* dx,
                          const int* fortran_comp);

    void mfix_multi_deposit_cic(const amrex_particle_real*, int ns, int np,
                                amrex_real* mf_bx, amrex_real* mf_by, amrex_real* mf_bz,
                                amrex_real* mf_bu, amrex_real* mf_bv, amrex_real* mf_bw,
                                const int* lo_x, const int* hi_x,
                                const int* lo_y, const int* hi_y,
                                const int* lo_z, const int* hi_z,
                                const amrex_real* plo, const amrex_real* dx,
                                const int* beta_comp, const int* beta_vel_comp);

    void mfix_interpolate_cic(const amrex_particle_real*, int ns, int np,
                              const amrex::Real* acc, const int* lo, const int* hi, int ncomp,
                              const amrex::Real* plo, const amrex::Real* dx);

    void calc_particle_beta(
       const int* slo, const int* shi,
       const int* ulo, const int* uhi,
       const int* vlo, const int* vhi,
       const int* wlo, const int* whi,
       const int* np,
       amrex::Real* ep_g, amrex::Real* ro_g,
       amrex::Real* u_g,  amrex::Real* v_g,   amrex::Real* w_g,
       amrex::Real* mu_g,
       void* particles,
       amrex::Real* dx,       amrex::Real* dy,       amrex::Real* dz);

    void calc_drag_particle(
	const int* slo, const int* shi,
	const int* ulo, const int* uhi,
	const int* vlo, const int* vhi,
	const int* wlo, const int* whi,
	const int* np,
	amrex::Real* p_g, amrex::Real* p0_g,
	amrex::Real* u_g,  amrex::Real* v_g,   amrex::Real* w_g,
	void* particles,
	amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
	amrex::Real* xlen, amrex::Real* ylen, amrex::Real* zlen);
    
    void calc_coeff(
	const int* slo, const int* shi, const int* lo,  const int* hi,int* level,
	amrex::Real* ro_g, amrex::Real* p_g, amrex::Real* p0_g, amrex::Real* ep_g,
	amrex::Real* rop_g);

    void calc_trd_g(const int* slo, const int* shi,
        const int* ulo, const int* uhi,
        const int* vlo, const int* vhi,
        const int* wlo, const int* whi,
        const int* lo,  const int* hi,
        amrex::Real* trd_g,
        amrex::Real* u_g, amrex::Real* v_g, amrex::Real* w_g,
        amrex::Real* dx,  amrex::Real* dy,  amrex::Real* dz);

    void calc_tau_g(
       const int* slo, const int* shi,
       const int* ulo, const int* uhi,
       const int* vlo, const int* vhi,
       const int* wlo, const int* whi,
       const int* lo,  const int* hi,
       amrex::Real* tau_u_g,  amrex::Real* tau_v_g,  amrex::Real* tau_w_g,
       amrex::Real* u_g,      amrex::Real* v_g,      amrex::Real* w_g,
       amrex::Real* trd_g,    amrex::Real* lambda_g, amrex::Real* mu_g,
       amrex::Real* dx,       amrex::Real* dy,       amrex::Real* dz);

    void conv_rop(
       const int*  lo,       const int*  hi,
       amrex::Real* rop_g, const int*  slo, const int*  shi,
       amrex::Real*   u_g, const int*  ulo, const int*  uhi,
       amrex::Real*   v_g, const int*  vlo, const int*  vhi,
       amrex::Real*   w_g, const int*  wlo, const int*  whi,
       amrex::Real* ropX , const int* rxlo, const int* rxhi,
       amrex::Real* ropY , const int* rylo, const int* ryhi,
       amrex::Real* ropZ , const int* rzlo, const int* rzhi);

    void calc_mflux(
       const int* ulo, const int* uhi,
       const int* vlo, const int* vhi,
       const int* wlo, const int* whi,
       amrex::Real* u_g,            amrex::Real* v_g,            amrex::Real* w_g,
       amrex::Real* ropX,         amrex::Real* ropY,         amrex::Real* ropZ,
       amrex::Real* fluxX,        amrex::Real* fluxY,        amrex::Real* fluxZ,
       amrex::Real* dx,             amrex::Real* dy,             amrex::Real* dz);

    void physical_prop(const int* slo, const int* shi,
           const int* lo, const int* hi,
           int* level, amrex::Real* ro_g, amrex::Real* p_g,
           amrex::Real* ep_g, amrex::Real* rop_g);

    void fill_bc0(amrex::Real* s, const int* slo, const int* shi,
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
	const int* dlo, const int* dhi,
	const int*  lo, const int*  hi,
	amrex::Real* u_g,      amrex::Real* v_g,     amrex::Real* w_g,
	amrex::Real* u_go,     amrex::Real* p_g,     amrex::Real* p0_g,     amrex::Real* ro_g,
	amrex::Real* rop_g,    amrex::Real* rop_go,  amrex::Real* ep_g,
	amrex::Real* tau_u_g,  amrex::Real* d_e,
	amrex::Real* fluxX,  amrex::Real* fluxY,     amrex::Real* fluxZ,
	amrex::Real* mu_g,     amrex::Real* f_gds_u, amrex::Real* drag_u,
	amrex::Real* A_m,      amrex::Real* b_m,     amrex::Real* mask,
	int* bc_ilo,    int* bc_ihi,
	int* bc_jlo,    int* bc_jhi,
	int* bc_klo,    int* bc_khi,
	const int* domlo , const int* domhi,
	amrex::Real* dt,
	amrex::Real* dx,       amrex::Real* dy,       amrex::Real* dz,
	amrex::Real *num_u, amrex::Real *denom_u);



    void solve_v_g_star(
	const int* slo, const int* shi,
	const int* ulo, const int* uhi,
	const int* vlo, const int* vhi,
	const int* wlo, const int* whi,
	const int* alo, const int* ahi,
	const int* dlo, const int* dhi,
	const int* lo,  const int* hi,
	amrex::Real* u_g,      amrex::Real* v_g,     amrex::Real* w_g,
	amrex::Real* v_go,     amrex::Real* p_g,     amrex::Real* p0_g,     amrex::Real* ro_g,
	amrex::Real* rop_g,    amrex::Real* rop_go,  amrex::Real* ep_g,
	amrex::Real* tau_v_g,  amrex::Real* d_n,
	amrex::Real* fluxX,    amrex::Real* fluxY,   amrex::Real* fluxZ,
	amrex::Real* mu_g,     amrex::Real* f_gds_v, amrex::Real* drag_v,
	amrex::Real* A_m,      amrex::Real* b_m,     amrex::Real* mask,
	int* bc_ilo,    int* bc_ihi,
	int* bc_jlo,    int* bc_jhi,
	int* bc_klo,    int* bc_khi,
	const int* domlo , const int* domhi,
	amrex::Real* dt,
	amrex::Real* dx,       amrex::Real* dy,       amrex::Real* dz,
	amrex::Real *num_v, amrex::Real *denom_v);


    void solve_w_g_star(
	const int* slo, const int* shi,
	const int* ulo,  const int* uhi,
	const int* vlo,  const int* vhi,
	const int* wlo,  const int* whi,
	const int* alo, const int* ahi,
	const int* dlo, const int* dhi,
	const int* lo,  const int* hi,
	amrex::Real* u_g,      amrex::Real* v_g,     amrex::Real* w_g,
	amrex::Real* w_go,     amrex::Real* p_g,     amrex::Real* p0_g,     amrex::Real* ro_g,
	amrex::Real* rop_g,    amrex::Real* rop_go,  amrex::Real* ep_g,
	amrex::Real* tau_w_g,  amrex::Real* d_t,
	amrex::Real* fluxX,    amrex::Real* fluxY,   amrex::Real* fluxZ,
	amrex::Real* mu_g,     amrex::Real* f_gds_w, amrex::Real* drag_w,
	amrex::Real* A_m,      amrex::Real* b_m,     amrex::Real* mask,
	int* bc_ilo,    int* bc_ihi,
	int* bc_jlo,    int* bc_jhi,
	int* bc_klo,    int* bc_khi,
	const int* domlo , const int* domhi,
	amrex::Real* dt,
	amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
	amrex::Real *num_w, amrex::Real *denom_w);
    
    
    void solve_pp_g(
       const int* slo, const int* shi,
       const int* ulo,  const int* uhi,
       const int* vlo,  const int* vhi,
       const int* wlo,  const int* whi,
       const int* alo,  const int* ahi,
       const int* lo,  const int* hi,
       amrex::Real* u_g,      amrex::Real* v_g,      amrex::Real* w_g,
       amrex::Real* p_g,      amrex::Real* ep_g,
       amrex::Real* rop_g,    amrex::Real* rop_go,
       amrex::Real* ro_g,
       amrex::Real* ropX,   amrex::Real* ropY,   amrex::Real* ropZ,
       amrex::Real* d_e,      amrex::Real* d_n,      amrex::Real* d_t,
       amrex::Real* A_m,      amrex::Real* b_m,      amrex::Real* b_mmax,
       int* bc_ilo, int* bc_ihi,
       int* bc_jlo, int* bc_jhi,
       int* bc_klo, int* bc_khi,
       amrex::Real* dt,       amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
       const int* domlo , const int* domhi,
       amrex::Real *num_p, amrex::Real *denom_p);

    void correct_p_0(
       const int*  lo,  const int*  hi,
       const int* slo,  const int* shi,
       amrex::Real* p_g,      amrex::Real* pp_g);

    void correct_u_0(
       const int*  lo,  const int*  hi,
       const int* ulo,  const int* uhi,
       const int* slo,  const int* shi,
       amrex::Real* pp_g,      amrex::Real* u_g,      amrex::Real* d_e);

    void correct_v_0(
       const int*  lo,  const int*  hi,
       const int* vlo,  const int* vhi,
       const int* slo,  const int* shi,
       amrex::Real* pp_g,      amrex::Real* v_g,      amrex::Real* d_n);

    void correct_w_0(
       const int*  lo,  const int*  hi,
       const int* wlo,  const int* whi,
       const int* slo,  const int* shi,
       amrex::Real* pp_g,      amrex::Real* w_g,      amrex::Real* d_t);

    void display_resid(amrex::Real* time, amrex::Real* dt, int* nit, amrex::Real *residuals);

    int check_convergence(int* nit, amrex::Real *residuals);

    int goal_seek_mflux(
       const int* slo, const int* shi,
       int *nit,       int* gsmf,      amrex::Real* delP_MF,  amrex::Real* lMFlux,
       amrex::Real* fluxX,  amrex::Real* fluxY,  amrex::Real* fluxZ,
       amrex::Real* dx     ,  amrex::Real* dy,       amrex::Real* dz);

    int adjustdt(int* ier, int* nit, amrex::Real* dt);

    void call_usr3_des( const int* np, void* particles );

    void call_usr2_des( const int* np, void* particles );

    void des_init_time_loop ( const  amrex::Real* tstart, const amrex::Real* dt, int* nsubsteps, amrex::Real* subdt, 
            const int* subdt_io );

    void des_time_loop ( const int* nrp,  void* rparticles,
                         const int* nf, amrex::Real* tow, amrex::Real* fc,
                         amrex::Real* subdt, amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
                         amrex::Real* stime, int* nstep);

    void particle_get_position (void* particles, const int* nrp, amrex::Real* x, amrex::Real* y, amrex::Real* z);

    void calc_particle_collisions ( void* rparticles, const int* nrp,
                                    void* gparticles, const int* ngp,
                                    const int* nbor_list, const int* size_nl,
                                    amrex::Real* tow, amrex::Real* fc,
                                    amrex::Real* subdt, int* ncoll);

    void calc_wall_collisions ( void* rparticles, const int* ntot, const int* nrp,
             amrex::Real* tow, amrex::Real* fc, amrex::Real* subdt,
             const void* flag, const int* fglo, const int* fghi,
             const amrex_real* normal, const int* nlo, const int* nhi,
             const amrex_real* bcent, const int* blo, const int* bhi,
             const amrex::Real* dx);

    void compute_normals (
             const int* lo , const int* hi,
             const void* flag, const int* fglo, const int* fghi,
             amrex_real* normal, const int* nlo, const int* nhi,
             const amrex_real* ax, const int* axlo, const int* axhi,
             const amrex_real* ay, const int* aylo, const int* ayhi,
             const amrex_real* az, const int* azlo, const int* azhi);

    int des_continuum_coupled ();

    void mfix_particle_write( const int* nrp,  void* particles);

    void mfix_particle_generator(int* pcount, const int* domlo , const int* domhi,
                                 amrex::Real* dx, amrex::Real* dy, amrex::Real* dz);


    void rm_wall_collisions ( void* rparticles, const int* nrp,
                              const void* flag, const int* fglo, const int* fghi,
                              const amrex_real* bcent, const int* blo, const int* bhi,
                              const amrex_real* ax, const int* axlo, const int* axhi,
                              const amrex_real* ay, const int* aylo, const int* ayhi,
                              const amrex_real* az, const int* azlo, const int* azhi,
                              const amrex::Real* dx);



    void mfix_particle_generator_prop( const int* nrp,  void* particles);

  void mfix_get_walls(int* bcv, int* exists,
                      amrex::RealVect* normal, amrex::RealVect* center);

    void mfix_usr0();

    void mfix_usr1();

    void mfix_usr2();

    void mfix_usr3(amrex::Real* u_g, const int*, const int*,
       amrex::Real* v_g, const int*, const int*,
       amrex::Real* w_g, const int*, const int*,
       amrex::Real* p_g, const int*, const int*,
       amrex::Real* dx ,     amrex::Real* dy,   amrex::Real* dz);

    void mfix_finl_err_msg();

    void leq_residual(const int* lo, const int* hi,
          const amrex::Real*, const int*, const int*,
          const amrex::Real*, const int*, const int*,
          const amrex::Real*, const int*, const int*,
          const amrex::Real*, const int*, const int*);

    void out_array(const amrex::Real*, const int*, const int*);
    void out_matrix(const amrex::Real*, const int*, const int*);

    void leq_matvec(const int* lo, const int* hi,
        const amrex::Real*, const int*, const int*,
        const amrex::Real*, const int*, const int*,
        const amrex::Real*, const int*, const int*);

    void leq_scale(const int* lo, const int* hi,
       const amrex::Real*, const int*, const int*,
       const amrex::Real*, const int*, const int*);

    void leq_msolve1(const int* lo, const int* hi,
         const amrex::Real*, const int*, const int*,
         const amrex::Real*, const int*, const int*,
         const amrex::Real*, const int*, const int*);

    void get_solver_params (int* eq_id, int* sweep, int* pc, int* maxiter, amrex::Real* tol);

    
    void compute_divu (
	const int* lo,  const int*  hi,
	amrex::Real* divu,  const int* slo, const int* shi,
	amrex::Real* u_g,   const int* ulo, const int* uhi,
	amrex::Real* v_g,   const int* vlo, const int* vhi,
	amrex::Real* w_g,   const int* wlo, const int* whi,
	const amrex::Real* dx );

    
    void compute_diveu (
	const int* lo,  const int*  hi,
	      amrex::Real* diveu,  const int* slo, const int* shi,
	const amrex::Real* ep_g,
	const amrex::Real* u_g,    const int* ulo, const int* uhi,
	const amrex::Real* v_g,    const int* vlo, const int* vhi,
	const amrex::Real* w_g,    const int* wlo, const int* whi,
	const amrex::Real* dx );

    
    void compute_vort (
	const int* lo,  const int*  hi,
	amrex::Real* vort,  const int* slo, const int* shi,
	amrex::Real* u_g,   const int* ulo, const int* uhi,
	amrex::Real* v_g,   const int* vlo, const int* vhi,
	amrex::Real* w_g,   const int* wlo, const int* whi,
	const amrex::Real* dx );


    void compute_bcoeff (
	const int* lo,  const int*  hi,
	      amrex::Real* bcoeff, const int* alo, const int* ahi,
	const amrex::Real* rop_g, const int* slo, const int* shi,
	const amrex::Real* ep_g, const int* dir );
    
    void compute_new_dt ( const amrex::Real* umax,
			  const amrex::Real* vmax,
			  const amrex::Real* wmax,
			  const amrex::Real* romin,
			  const amrex::Real* mumax,
			  const amrex::Real* dx,
			  const amrex::Real* cfl,
			  amrex::Real*  dt ); 

    void init_periodic_vortices (
	const int* utlo, const int* uthi,
	amrex::Real* u,  const int* ulo, const int* uhi,
        const int* vtlo, const int* vthi,
	amrex::Real* v,  const int* vlo, const int* vhi,
	const int* wtlo, const int* wthi,
	amrex::Real* w,  const int* wlo, const int* whi,
	const amrex::Real* dx,
	const int*  domlo, const int* plane );

  void add_grad_p (
	const int* lo, const int* hi,
	amrex::Real* u_i, const int* ulo, const int* uhi,
        const amrex::Real* ro_g, const int* slo, const int* shi,
	const amrex::Real* p, const amrex::Real* p0,
	const amrex::Real* dx, const amrex::Real* c,
	const int* dir );

  void add_grad_phi (
	const int* lo, const int* hi,
	amrex::Real* u_i, const int* ulo, const int* uhi,
	const amrex::Real* ro_g,
        const amrex::Real* phi, const int* slo, const int* shi,
	const amrex::Real* dx, const amrex::Real* c,
	const int* dir );    

    void compute_fluid_acceleration (
	const int* lo, const int* hi,
	      amrex::Real* rhs, const int* rlo, const int* rhi, 
	const amrex::Real* sl,
	const amrex::Real* u,   const int* ulo, const int* uhi, 
	const amrex::Real* v,   const int* vlo, const int* vhi, 
	const amrex::Real* w,   const int* wlo, const int* whi, 
	const amrex::Real* mu,  const int* slo, const int* shi, 
	const amrex::Real* rop, const amrex::Real* dx,
	const int* dir );

    void compute_u_slopes (
	const int* lo, const int* hi,
	const amrex::Real* u,   const int* ulo, const int* uhi, 
	amrex::Real* slopes,
	const int* domlo, const int* domhi,
	const int* bc_ilo_type, const int* bc_ihi_type );

    void compute_v_slopes (
	const int* lo, const int* hi,
	const amrex::Real* v,   const int* vlo, const int* vhi, 
	amrex::Real* slopes,
	const int* domlo, const int* domhi,
	const int* bc_jlo_type, const int* bc_jhi_type );
    
    void compute_w_slopes (
	const int* lo, const int* hi,
	const amrex::Real* w, const int* wlo, const int* whi, 
	amrex::Real* slopes,
	const int* domlo, const int* domhi,
	const int* bc_klo_type, const int* bc_khi_type );

    void set_ppe_bc (
	int* bc_lo, int* bc_hi,
	const int* domlo, const int* domhi,
	const int* bct_ilo, const int* bct_ihi, 
	const int* bct_jlo, const int* bct_jhi,
	const int* bct_klo, const int* bct_khi,
	int* singular );
    
    void set_poisson_solver_bcs (
	const int* lo, const int* hi,
	amrex::Real* rhs, const int* slo, const int* shi,
	amrex::Real* phi, 
        const amrex::Real* oro_x, const int* ulo, const int* uhi,
        const amrex::Real* oro_y, const int* vlo, const int* vhi,
	const amrex::Real* oro_z, const int* wlo, const int* whi,
	const int* bc_lo, const int* bc_hi,
	const amrex::Real* dx, const int* domlo, const int* domhi );
    
    void add_forcing (
	const int* lo, const int* hi,
	      amrex::Real* u_i,  const int* ulo, const int* uhi,
	const amrex::Real* drag_i, const int* dlo, const int* dhi,
	const amrex::Real* ro_g, const int* slo, const int* shi,
	const amrex::Real* rop_g,
	const int* domlo, const int* domhi,
	const amrex::Real* dx, const amrex::Real* dt, 
	const int* dir);

    void compute_intermediate_velocity (
	const int* lo, const int* hi,
	      amrex::Real* u_i,  const int* ulo, const int* uhi,
	const amrex::Real* f_gds_i, const int* flo, const int* fhi,
	const amrex::Real* rop, const int* slo, const int* shi,
	const int* dir, const amrex::Real* dt );


    void extrap_pressure_to_ghost_cells (
	amrex::Real* p_g, const int* slo, const int* shi,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
	const int* domlo , const int* domhi);

    void set_projection_bcs (
	amrex::Real* u_g, const int* ulo, const int* uhi,
	amrex::Real* v_g, const int* vlo, const int* vhi,
	amrex::Real* w_g, const int* wlo, const int* whi,
	amrex::Real* ep_g, const int* slo, const int* shi,
	amrex::Real* ro_g,    amrex::Real* rop_g,
	amrex::Real* mu_g,    amrex::Real* lambda_g,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
	const int* domlo , const int* domhi);
    
    void mfix_eb_to_polygon ( const amrex::Real* dx,
			      const int* slo, const int* shi,
			      const void* flag, const int* fglo, const int* fghi,
			      const amrex_real* bcent, const int* blo, const int* bhi,
			      const amrex_real* ax, const int* axlo, const int* axhi,
			      const amrex_real* ay, const int* aylo, const int* ayhi,
			      const amrex_real* az, const int* azlo, const int* azhi);

    void mfix_write_eb_vtp ( int* myID);
    void mfix_write_pvtp ( int* nProcs);

    void init_levelset(const int * lo,   const int * hi,
                       amrex_real * phi, const int * phlo, const int * phhi);

    void fill_levelset_eb(const int * lo,             const int * hi,
                          const amrex_real * eb_list, const int * l_eb,
                          int * valid,                const int * vlo,     const int * vhi,
                          amrex_real * phi,           const int * phlo,    const int * phhi, 
                          const amrex_real * dx,      const amrex_real * dx_eb);

    void update_levelset(const int * lo,           const int * hi,
                         const amrex_real * ls_in, const int * lslo,    const int * lshi,
                         int * valid,              const int * vlo,     const int * vhi,
                         amrex_real * phi,         const int * phlo,    const int * phhi,
                         const amrex_real * dx,    const int * n_pad);

    void calc_wall_collisions_ls(void * rparticles, const int * ntot, const int * nrp,
                                 amrex_real * tow, amrex_real * fc, amrex_real * subdt,
                                 const int * valid, const int * vlo, const int * vhi,
                                 const amrex_real * phi, const int * phlo, const int * phhi,
                                 const amrex_real * dx, const int * n_refine);

    void count_eb_facets(const int * lo,    const int * hi,
                         const void * flag, const int * flo, const int * fhi,
                         int * n_facets);

    void eb_as_list(const int * lo,           const int * hi,  int * c_facets,
                    const void * flag,        const int * flo, const int * fhi,
                    const amrex_real * norm,  const int * nlo, const int * nhi,
                    const amrex_real * bcent, const int * blo, const int * bhi,
                    amrex_real * list_out,    const int * lsize,
                    const amrex_real * dx);

    void validate_levelset(const int * lo,          const int * hi,   const int * n_pad,
                           const amrex_real * impf, const int * imlo, const int * imhi,
                           const int * valid,       const int * vlo,  const int * vhi,
                           amrex_real * phi,        const int * phlo, const int * phhi);


#ifdef __cplusplus
}
#endif
