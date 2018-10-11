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

    void set_delp_dir(int* delp_dir);

    void mfix_get_data( 
        int* fluid, int* dem,
        int* steady_state, amrex::Real* dt, amrex::Real* dt_min, amrex::Real* dt_max,
        amrex::Real* tstop, int* call_udf);

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
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
        amrex::Real* xlen, amrex::Real* ylen, amrex::Real* zlen,
        const int* ng
    );

    void set_bc0(
        const int* slo, const int* shi,
        amrex::Real* ep_g,
        amrex::Real* ro_g,  amrex::Real* rop_g,
        amrex::Real* mu_g,  amrex::Real* lambda_g,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng    , const int* nodal_pressure
    );

    void init_fluid(const int* slo, const int* shi,
        const int*  lo, const int*  hi,
        const int* domlo , const int* domhi,
        amrex::Real* ep_g,  amrex::Real* ro_g,
        amrex::Real* rop_g, amrex::Real* p_g,
        amrex::Real* vel_g,
        amrex::Real* mu_g,  amrex::Real* lambda_g,
        amrex::Real* dx,    amrex::Real* dy,       amrex::Real* dz,
        amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength
       );

    void init_fluid_restart( 
          const int* slo, const int* shi,
          const int*  lo, const int*  hi,
          amrex::Real* mu_g,    amrex::Real* lambda_g
       );

    void set_p0(const int*  lo   , const int*  hi,
	        const int* domlo , const int* domhi,
                amrex::Real* p0_g , const int* rlo, const int* rhi,
                amrex::Real* gp0_g, const int* glo, const int* ghi,
	        amrex::Real* dx,       amrex::Real* dy,       amrex::Real* dz,
	        amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
	        int* delp_dir,
	        const int* bc_ilo, const int* bc_ihi,
        	const int* bc_jlo, const int* bc_jhi,
        	const int* bc_klo, const int* bc_khi,
                const int* ng, const int* nodal_pressure
       );

    void set_velocity_bcs (
        amrex::Real* vel_g, const int* ulo, const int* uhi,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
	const int* domlo , const int* domhi,
        const int* ng, const int* extrap_dir_bcs
       );

    void set_vec_bcs (
        amrex::Real* vec, const int* ulo, const int* uhi,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
	const int* domlo , const int* domhi,
        const int* ng
       );

    void flip_particle_vol(const int* slo, const int* shi,
        amrex::Real* ep_g,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng
        );

  void mfix_cap_eps(const int* slo, const int* shi, amrex::Real* ep_g);


    void set_drag_bcs(
        amrex::Real* drag, const int* slo, const int* shi,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng
        );

    void zero_wall_norm_vel(
        const int* slo, const int* shi,
        amrex::Real* vel_g,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng
        );

    void set_particle_properties(
            const int* pstate, amrex::Real* pradius, amrex::Real* pdensity,
            amrex::Real* pvol, amrex::Real* pmass, amrex::Real* omoi, amrex::Real* omega
        );

    void init_collision(amrex::Real *sum_dp, amrex::Real *sum_roi
        );

    void sum_particle_props(
            const int* nrp, void* particles,
            amrex::Real *sum_np, amrex::Real *sum_dp, amrex::Real *sum_ro
        );

    void output_manager ( const int* np,  amrex::Real* time, amrex::Real* dt,
             amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
             int* nstep, void* particles, int* finsh
        );

    void calc_volume_fraction (
        const int*  lo, const int*  hi,
        const int* slo, const int* shi,
        const int*  np, void* particles,
        amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
        amrex::Real* ep_g 
        );

    void mfix_deposit_cic(const amrex_particle_real*, int ns, int np, int nc,
                          amrex_real* mf, const int* lo, const int* hi,
                          const amrex_real* plo, const amrex_real* dx,
                          const int* fortran_comp
        );

    void mfix_multi_deposit_cic(const amrex_particle_real*, int ns, int np,
                                amrex_real* mf_bx, amrex_real* mf_bu,
                                const int* lo, const int* hi,
                                const amrex_real* plo, const amrex_real* dx,
                                const int* beta_comp, const int* beta_vel_comp
        );

    void calc_particle_beta(
       const int* slo, const int* shi,
       amrex::Real*  ep_g, amrex::Real* ro_g,
       amrex::Real* vel_g, amrex::Real* mu_g,
       const int* np, void* particles,
       amrex::Real* dx, amrex::Real* dy, amrex::Real* dz
        );

    void calc_drag_particle(  
            const amrex::Real* gp , const int* gplo , const int* gphi,
            const amrex::Real* gp0, const int* gp0lo, const int* gp0hi,
            const amrex::Real* vel, const int* ulo, const int* uhi,
            const amrex::Real* xsl, const amrex::Real* ysl, 
            const amrex::Real* zsl, const int* slo, const int* shi,
            const int* np, void* particles, const amrex::Real* dx, 
            const amrex::Real* x0, int interp_type
        );

    void physical_prop(const int* slo, const int* shi,
           const int* lo, const int* hi,
           int* level, amrex::Real* ro_g, amrex::Real* p_g,
           amrex::Real* ep_g, amrex::Real* rop_g
        );

    void fill_bc0(amrex::Real* s, const int* slo, const int* shi,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng
        );

    void call_usr3_des( const int* np, void* particles );

    void call_usr2_des( const int* np, void* particles );

    void des_init_time_loop ( const  amrex::Real* tstart, const amrex::Real* dt, int* nsubsteps, amrex::Real* subdt, int* subdt_io );

    void des_time_loop ( const int* nrp,  void* rparticles,
                         const int* nf, amrex::Real* tow, amrex::Real* fc,
                         amrex::Real* subdt, amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
                         amrex::Real* stime, int* nstep
        );

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

    int des_continuum_coupled ();

    void mfix_particle_write( const int* nrp,  void* particles);

    void mfix_particle_generator(int* pcount, const int* domlo , const int* domhi,
                                 amrex::Real* dx, amrex::Real* dy, amrex::Real* dz);


    void mfix_particle_generator_prop( const int* nrp,  void* particles);

    void mfix_usr0();

    void mfix_usr1(amrex::Real* time);

    void mfix_usr2();

    void mfix_usr3(
            amrex::Real* u_g, const int*, const int*,
            amrex::Real* v_g, const int*, const int*,
            amrex::Real* w_g, const int*, const int*,
            amrex::Real* p_g, const int*, const int*,
            amrex::Real* dx, amrex::Real* dy, amrex::Real* dz);

    void mfix_finl_err_msg();

    void construct_gradp( 
	    const int*  lo,  const int*  hi,
            amrex::Real* p_g, 
            amrex::Real* p0_g, const int* rlo, const int* rhi,
            amrex::Real* gpx, const int* xlo, const int* xhi,
            amrex::Real* gpy, const int* ylo, const int* yhi,
            amrex::Real* gpz, const int* zlo, const int* zhi,
	    amrex::Real* dx, amrex::Real* dy, amrex::Real* dz, 
	    const int* bc_ilo, const int* bc_ihi,
    	    const int* bc_jlo, const int* bc_jhi,
	    const int* bc_klo, const int* bc_khi,
            const int* domlo , const int* domhi,
            const int* ng, const int* nodal_pressure
         );

    void set_gradp_bcs (
             const int* slo, const int* shi,
             amrex::Real* gp, const int* glo, const int* ghi,
	     const int* bc_ilo, const int* bc_ihi,
	     const int* bc_jlo, const int* bc_jhi,
	     const int* bc_klo, const int* bc_khi,
	     const int* domlo , const int* domhi,
             const int* ng
       );

    void init_periodic_vortices (
	const int* lo,  const int*  hi,
	amrex::Real* vel,  const int* ulo, const int* uhi,
	const amrex::Real* dx,
	const int*  domlo, const int* plane
       );

    void mfix_collect_fluid(
            const int*  lo, const int*  hi,
            const int* domlo , const int* domhi,
            amrex::Real* p   , const int* p_lo, const int* p_hi,
            amrex::Real* ep_g, const int* s_lo, const int* s_hi,
            amrex::Real* dx,   amrex::Real* dy,  amrex::Real* dz,
            amrex::Real (&accumulator)[256]
            );

    void mfix_write_fluid(const int* domlo , const int* domhi,
            amrex::Real* dx,   amrex::Real* dy,  amrex::Real* dz,
            amrex::Real* time, amrex::Real* dt,
            amrex::Real (&accumulator)[256]
            );
   
#ifdef __cplusplus
}
#endif
