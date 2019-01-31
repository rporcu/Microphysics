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

    void mfix_sum_mf ( const int* lo, const int* hi,
                       const BL_FORT_FAB_ARG_3D(rho),
                       const amrex::Real* dx,
                       const BL_FORT_FAB_ARG_3D(vol),
                       amrex::Real* mass );

    void compute_vort (
        const int* lo,  const int*  hi,
        amrex::Real* vort , const int* slo, const int* shi,
        amrex::Real* vel_g, const int* ulo, const int* uhi,
        const amrex::Real* dx );

    void extrap_pressure_to_ghost_cells (
       amrex::Real* p_g, const int* slo, const int* shi,
       const int* bc_ilo, const int* bc_ihi,
       const int* bc_jlo, const int* bc_jhi,
       const int* bc_klo, const int* bc_khi,
       const int* domlo , const int* domhi,
       const int* nghost
       );

   void set_scalar_bcs (
      amrex::Real* ep_g, const int* slo, const int* shi,
      amrex::Real* ro_g,    amrex::Real* rop_g,
      amrex::Real* mu_g,    amrex::Real* lambda_g,
      const int* bc_ilo, const int* bc_ihi,
      const int* bc_jlo, const int* bc_jhi,
      const int* bc_klo, const int* bc_khi,
      const int* domlo , const int* domhi,
      const int* nghost
      );

   void set_pressure_bcs (
      amrex::Real* p_g, const int* slo, const int* shi,
      const int* bc_ilo, const int* bc_ihi,
      const int* bc_jlo, const int* bc_jhi,
      const int* bc_klo, const int* bc_khi,
      const int* domlo , const int* domhi);

    void set_velocity_bcs (
        amrex::Real* time,
        amrex::Real* vel_g, const int* ulo, const int* uhi,
        const int* bc_ilo, const int* bc_ihi,
        const int* bc_jlo, const int* bc_jhi,
        const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng, const int* extrap_dir_bcs
       );


   void reconstruct_velocity (
       amrex::Real* vel_out,        const int* volo, const int* vohi,
       const amrex::Real* vel_in,   const int* vilo, const int* vihi,
       const amrex::Real* phi,      const int* phlo, const int* phhi,
       const int n_refine,
       const void *     flags,      const int*  flo, const int*  fhi,
       const amrex::Real* x0, const amrex::Real* dx );

    void phifill (BL_FORT_FAB_ARG_3D(state),
                  const int* dlo, const int* dhi,
                  const amrex_real* dx, const amrex_real* glo,
                  const amrex_real* time, const int* bc);
#ifdef __cplusplus
}
#endif
