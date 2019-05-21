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

   void set_scalar_bcs (
      amrex::Real* ep_g, const int* slo, const int* shi,
      amrex::Real* ro_g, 
      amrex::Real* mu_g,
      const int* bc_ilo, const int* bc_ihi,
      const int* bc_jlo, const int* bc_jhi,
      const int* bc_klo, const int* bc_khi,
      const int* domlo , const int* domhi,
      const int* nghost
      );

    void set_velocity_bcs (
        amrex::Real* time,
        amrex::Real* vel_g, const int* ulo, const int* uhi,
        const int* bc_ilo, const int* bc_ihi,
        const int* bc_jlo, const int* bc_jhi,
        const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng, const int* extrap_dir_bcs
       );

    void phifill (BL_FORT_FAB_ARG_3D(state),
                  const int* dlo, const int* dhi,
                  const amrex_real* dx, const amrex_real* glo,
                  const amrex_real* time, const int* bc);
#ifdef __cplusplus
}
#endif
