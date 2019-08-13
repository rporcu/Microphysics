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

    void mfix_get_walls(
            int* bcv, int* exists,
            void* normal, void* center
        );

    void mfix_get_real_walls(
            int* bcv, int* exists,
            void* normal, void* center
        );

    void ls_has_walls(int * has_wall,
                      const amrex_real * phi, const int * phlo, const int * phhi,
                      const amrex_real * tol
        );
    
    void calc_wall_collisions(
        void * rparticles, const int * ntot, const int * nrp,
        amrex_real * tow, amrex_real * fc, amrex_real * subdt,
        const int * valid,      const int * vlo,  const int * vhi,
        const amrex_real * phi, const int * phlo, const int * phhi,
        const amrex_real * dx, const int * n_refine
        );

#ifdef __cplusplus
}
#endif
