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

    void mfix_eb_to_polygon(
            const amrex::Real* dx,
            const int* slo, const int* shi,
            const void* flag,        const int* fglo, const int* fghi,
            const amrex_real* bcent, const int* blo,  const int* bhi,
            const amrex_real* ax,    const int* axlo, const int* axhi,
            const amrex_real* ay,    const int* aylo, const int* ayhi,
            const amrex_real* az,    const int* azlo, const int* azhi
        );

    void mfix_write_eb_vtp(int* myID);
    void mfix_write_pvtp(int* nProcs);

    void mfix_eb_grid_coverage(
            int* myID, const amrex::Real* dx,
            const int* slo, const int* shi,
            const void* flag, const int* fglo, const int* fghi
        );

    void calc_wall_collisions_ls(
            void * rparticles, const int * ntot, const int * nrp,
            amrex_real * tow, amrex_real * fc, amrex_real * subdt,
            const int * valid,      const int * vlo,  const int * vhi,
            const amrex_real * phi, const int * phlo, const int * phhi,
            const amrex_real * dx, const int * n_refine
        );

    void rm_wall_collisions(
            void* rparticle,  const int* nrp,
            const int * valid,      const int * vlo,  const int * vhi,
            const amrex_real * phi, const int * phlo, const int * phhi,
            const amrex_real * dx,  const int * n_refine
        );
    
    void rm_wall_collisions_eb(
            void* rparticle,  const int* nrp,
            const int*        valid, const int*  vlo, const int*  vhi,
            const amrex_real*   phi, const int* phlo, const int* phhi,
            const void*       flags, const int*  flo, const int*  fhi, 
            const amrex_real* plo, const amrex_real * dx,
            const int * n_refine
        );
    

    // According to Weiqun `fix_up_geometry` is deprecated.
    //void mfix_fixup_eb_geom (const int* lo, const int* hi,
    //                         const void* flag, const int* fglo, const int* fghi,
    //                         const amrex_real* volfrac, const int* vlo, const int* vhi,
    //                         const int* domlo, const int* domhi
    //    );

#ifdef __cplusplus
}
#endif
