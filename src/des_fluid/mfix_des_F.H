#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef BL_SINGLE_PRECISION_PARTICLES
typedef float amrex_particle_real;
#else
typedef double amrex_particle_real;
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    void get_collision_coefficients(int* nphase, amrex::Real* kt, amrex::Real* kt_w,
                                    amrex::Real* kn, amrex::Real* kn_w, amrex::Real* etan,
                                    amrex::Real* etan_w, amrex::Real* etat,amrex::Real* etat_w);
    
    int des_continuum_coupled ();

    void get_gravity (
        amrex::Real* grav
        );

    void calc_volume_fraction(
        const int*  lo, const int*  hi,
        const int* slo, const int* shi,
        const int*  np, void* particles,
        amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
        amrex::Real* ep_g
        );


    void mfix_deposit_cic_eb(
        const amrex_particle_real*, int ns, int np,
        amrex_real* vol,          const int* vlo, const int* vhi,
        const amrex_real* vratio, const int* slo, const int* shi,
        const void*        flags, const int* flo, const int* fhi,
        const amrex_real* plo, const amrex_real* dx,
        const int* fortran_comp
        );


    void mfix_multi_deposit_cic_eb(
        const amrex_particle_real*, int ns, int np,
        amrex_real* mf_bx, amrex_real* mf_bu,
        const int* lo, const int* hi,
        const amrex_real* vratio, const int* vlo, const int* vhi,
        const void*        flags, const int* flo, const int* fhi,
        const amrex_real* plo, const amrex_real* dx,
        const int* beta_comp, const int* beta_vel_comp
        );

    void mfix_interpolate_cic(
        const amrex_particle_real*, int ns, int np,
        const amrex::Real* acc, const int* lo, const int* hi, int ncomp,
        const amrex::Real* plo, const amrex::Real* dx
        );

    void calc_particle_beta(
        const amrex::Real*  ep_g, const int* slo, const int* shi,
        const amrex::Real* ro_g,  const amrex::Real* mu_g,
        const amrex::Real* vel_g, const int* ulo, const int* uhi,
        const int* np, void* particles,
        const amrex::Real* x0, const amrex::Real* dx
        );


    void calc_drag_particle(
        const amrex::Real* gp , const int* gplo , const int* gphi,
        const amrex::Real* gp0,
        const amrex::Real* vel, const int* ulo, const int* uhi,
        const int* np, void* particles, const amrex::Real* dx,
        const amrex::Real* x0
        );


    void calc_drag_particle_eb(
        const amrex::Real* gp , const int* gplo , const int* gphi,
        const amrex::Real* gp0,
        const amrex::Real* vel, const int* ulo,   const int* uhi,
        const void*  flags    , const int* flo,   const int* fhi,
        const int* np, void* particles, const amrex::Real* dx,
        const amrex::Real* x0
        );

    void calc_particle_collisions(
        void* rparticles, const int* nrp,
        void* gparticles, const int* ngp,
        const int* nbor_list, const int* size_nl,
        amrex::Real* tow, amrex::Real* fc,
        amrex::Real* subdt, int* ncoll
        );

    void flip_particle_vol(
        const int* slo, const int* shi,
        amrex::Real* ep_g,
        const int* bc_ilo, const int* bc_ihi,
        const int* bc_jlo, const int* bc_jhi,
        const int* bc_klo, const int* bc_khi,
        const int* domlo, const int* domhi,
        const int* nghost
        );

    void set_particle_properties(
        const int* pstate, amrex::Real* pradius, amrex::Real* pdensity,
        amrex::Real* pvol, amrex::Real* pmass, amrex::Real* omoi, amrex::Real* omega
        );

    void mfix_particle_write(const int* nrp, void* particles);

    void mfix_particle_generator(
        int* pcount, const int* domlo , const int* domhi,
        amrex::Real* dx, amrex::Real* dy, amrex::Real* dz
        );

    void mfix_particle_generator_prop(const int* nrp, void* particles);

    void particle_get_position (void* particles, const int* nrp, amrex::Real* x, amrex::Real* y, amrex::Real* z);

    void des_init_time_loop(
        const  amrex::Real* tstart, const amrex::Real* dt, int* nsubsteps, amrex::Real* subdt
        );

    void des_time_loop(
        const int* nrp,  void* rparticles,
        const int* nf, amrex::Real* tow, amrex::Real* fc,
        amrex::Real* subdt,
        amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
        amrex::Real* stime, int* nstep
        );

    void call_usr3_des(const int* np, void* particles);

    void call_usr2_des(const int* np, void* particles);

#ifdef __cplusplus
}
#endif
