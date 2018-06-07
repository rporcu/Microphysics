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
    void calc_volume_fraction(
            const int*  lo, const int*  hi,
            const int* slo, const int* shi,
            const int*  np, void* particles,
            amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
            amrex::Real* ep_g
        );

    void mfix_deposit_cic(
            const amrex_particle_real*, int ns, int np, int nc,
            amrex_real* mf, const int* lo, const int* hi,
            const amrex_real* plo, const amrex_real* dx,
            const int* fortran_comp
        );

    void mfix_multi_deposit_cic(
            const amrex_particle_real*, int ns, int np,
            amrex_real* mf_bx, amrex_real* mf_by, amrex_real* mf_bz,
            amrex_real* mf_bu, amrex_real* mf_bv, amrex_real* mf_bw,
            const int* lo_x, const int* hi_x,
            const int* lo_y, const int* hi_y,
            const int* lo_z, const int* hi_z,
            const amrex_real* plo, const amrex_real* dx,
            const int* beta_comp, const int* beta_vel_comp
        );

    void mfix_interpolate_cic(
            const amrex_particle_real*, int ns, int np,
            const amrex::Real* acc, const int* lo, const int* hi, int ncomp,
            const amrex::Real* plo, const amrex::Real* dx
        );

    void calc_wall_collisions(
            void* rparticles, const int* ntot, const int* nrp,
            amrex::Real* tow, amrex::Real* fc, amrex::Real* subdt,
            const void* flag, const int* fglo, const int* fghi,
            const amrex_real* normal, const int* nlo, const int* nhi,
            const amrex_real* bcent,  const int* blo, const int* bhi,
            const amrex::Real* dx
        );

    void rm_wall_collisions(
            void* rparticles, const int* nrp,
            const void* flag, const int* fglo, const int* fghi,
            const amrex_real* bcent, const int* blo, const int* bhi,
            const amrex_real* ax, const int* axlo, const int* axhi,
            const amrex_real* ay, const int* aylo, const int* ayhi,
            const amrex_real* az, const int* azlo, const int* azhi,
            const amrex::Real* dx
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

    void set_drag_bcs(
            const int* ulo, const int* uhi,
            const int* vlo, const int* vhi,
            const int* wlo, const int* whi,
            amrex::Real* drag_u, amrex::Real* drag_v, amrex::Real* drag_w,
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

    void des_init_time_loop(
            const  amrex::Real* tstart, const amrex::Real* dt, int* nsubsteps, amrex::Real* subdt,
            const int* subdt_io
        );

    void des_time_loop(
            const int* nrp,  void* rparticles,
            const int* nf, amrex::Real* tow, amrex::Real* fc,
            amrex::Real* subdt,
            amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
            amrex::Real* stime, int* nstep
        );

    void output_manager(
            const int* np, amrex::Real* time, amrex::Real* dt,
            amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
            int* nstep, void* particles, int* finsh
        );

    void call_usr3_des(const int* np, void* particles);

    void call_usr2_des(const int* np, void* particles);

#ifdef __cplusplus
}
#endif
