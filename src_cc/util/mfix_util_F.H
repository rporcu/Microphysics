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

    void average_cc_to_fc ( const int* lo, const int* hi,
                            amrex::Real* fcx, const int* fcxlo, const int* fcxhi,
                            amrex::Real* fcy, const int* fcylo, const int* fcyhi,
                            amrex::Real* fcz, const int* fczlo, const int* fczhi,
                            const amrex::Real* cc, const int* slo, const int* shi
      );

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
        amrex::Real* vel_g, const int* ulo, const int* uhi,
        const int* bc_ilo, const int* bc_ihi,
        const int* bc_jlo, const int* bc_jhi,
        const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng, const int* extrap_dir_bcs
       );

#ifdef __cplusplus
}
#endif
