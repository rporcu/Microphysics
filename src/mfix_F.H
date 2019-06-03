#ifndef _MFIX_F_H_
#define _MFIX_F_H_

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
    void mfix_add_argument(char* fname, int* nlen);

    void write_out0(amrex::Real* time, 
             amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
             amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
             const int* domlo, const int* domhi);

    void write_usr0();

    void set_ps(amrex::Real* dx, amrex::Real* dy, amrex::Real* dz, int* err, int* is_ioproc);

    void set_delp_dir(int* delp_dir);

    void mfix_get_data(
             int* fluid, int* dem,
             int* call_udf, int* name_len, const char *cmfix_dat);

    void check_inputs();

    void mfix_set_cyclic(int* cyc_x, int* cyc_y, int* cyc_z);

    void check_initial_conditions(amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
             const int* domlo, const int* domhi);

    void check_boundary_conditions(amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
             amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
             const int* domlo, const int* domhi);

    void check_point_sources(amrex::Real* dx, amrex::Real* dy, amrex::Real* dz);

    void check_bc_flow();

    void mfix_set_bc_mod(const int* pID, const int* pType,
             const amrex::Real* pLo, const amrex::Real* pHi,
             amrex::Real* pLoc,
             amrex::Real* pPg,
             amrex::Real* pVel);

  void mfix_set_bc_mod_add_mi(const int* pPlane,
                       amrex::Real* xLo, amrex::Real* yLo, amrex::Real* zLo,
                       amrex::Real* xHi, amrex::Real* yHi, amrex::Real* zHi,
                       amrex::Real* pPg, amrex::Real* pVel
                              );

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
             amrex::Real* ro_g,
             amrex::Real* mu_g,
             const int* bc_ilo, const int* bc_ihi,
             const int* bc_jlo, const int* bc_jhi,
             const int* bc_klo, const int* bc_khi,
             const int* domlo , const int* domhi,
             const int* ng
    );

    void set_ls_inflow(
        amrex::Real* ls_phi , const int* phi_lo, const int* phi_hi,
        const int* bc_ilo, const int* bc_ihi,
        const int* bc_jlo, const int* bc_jhi,
        const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng    , const int* n     ,
        const amrex::Real* dx
    );

    void set_p0(const int*  lo   , const int*  hi,
             const int* domlo , const int* domhi,
             amrex::Real* p0_g , const int* rlo, const int* rhi,
             amrex::Real* gp0,
             amrex::Real* dx,       amrex::Real* dy,       amrex::Real* dz,
             amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
             int* delp_dir,
             const int* bc_ilo, const int* bc_ihi,
             const int* bc_jlo, const int* bc_jhi,
             const int* bc_klo, const int* bc_khi,
             const int* ng
       );

    void mfix_cap_eps(const int* slo, const int* shi, amrex::Real* ep_g);

    void zero_wall_norm_vel(
             const int* slo, const int* shi,
             amrex::Real* vel_g,
             const int* bc_ilo, const int* bc_ihi,
             const int* bc_jlo, const int* bc_jhi,
             const int* bc_klo, const int* bc_khi,
             const int* domlo , const int* domhi,
             const int* ng
        );

    void init_collision(amrex::Real *sum_dp, amrex::Real *sum_roi
        );

    void fill_bc0(
             amrex::Real* s, const int* slo, const int* shi,
             const int* bc_ilo, const int* bc_ihi,
             const int* bc_jlo, const int* bc_jhi,
             const int* bc_klo, const int* bc_khi,
             const int* domlo , const int* domhi,
             const int* ng
        );

    void mfix_usr0();

    void mfix_usr1(amrex::Real* time);

    void mfix_usr2();

    void mfix_usr3(
             amrex::Real* vel_g, const int* ulo, const int* uhi,
             amrex::Real* p_g,   const int* slo, const int* shi,
             amrex::Real* dx, amrex::Real* dy, amrex::Real* dz);

    void mfix_finl_err_msg();

    void mfix_collect_fluid(
            const int*    lo,  const int*   hi,
            const int* domlo , const int* domhi,
            amrex::Real* ep_g,        const int* slo, const int* shi,
            amrex::Real* p_g,         const int* plo, const int* phi,
            amrex::Real* vel_g,       const int* vlo, const int* vhi,
            const amrex_real* vratio, const int* rlo, const int* rhi,
            const amrex::Real* x_w, const amrex::Real* x_e,
            const amrex::Real* y_s, const amrex::Real* y_n,
            const amrex::Real* z_b, const amrex::Real* z_t, const amrex::Real* dx,
            amrex::Real* sum_ep_g, amrex::Real* sum_p_g,  amrex::Real* sum_vol,
            amrex::Real* sum_velx, amrex::Real* sum_vely, amrex::Real* sum_velz
            );

    void mfix_write_fluid(const int* domlo , const int* domhi,
            amrex::Real* dx,   amrex::Real* dy,  amrex::Real* dz,
            amrex::Real* time, amrex::Real* dt,
            amrex::Real (&accumulator)[256]
            );

#ifdef __cplusplus
}
#endif

#endif
