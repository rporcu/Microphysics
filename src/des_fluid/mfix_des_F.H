#ifndef _MFIX_DES_F_H_
#define _MFIX_DES_F_H_

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

// TODO move this to another file
// Boundary conditions list auxiliary structure
typedef struct BcList {
  public:
    BcList(): minf(0), pinf(0), pout(0) {};
    int minf;
    int pinf;
    int pout;
} BcList;

#ifdef BL_SINGLE_PRECISION_PARTICLES
typedef float amrex_particle_real;
#else
typedef double amrex_particle_real;
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    void get_collision_model(int* coll_model_flag);
    
    void get_lsd_collision_coefficients(int* nphase,
                                        amrex::Real* kt, amrex::Real* kt_w,
                                        amrex::Real* kn, amrex::Real* kn_w,
                                        amrex::Real* mew, amrex::Real* mew_w,
                                        amrex::Real* etan, amrex::Real* etan_w,
                                        amrex::Real* etat, amrex::Real* etat_w);
    
    int des_continuum_coupled ();

    void get_gravity (
        amrex::Real* grav
        );

    // TODO move this to another file
    void get_bc_list (
        BcList& bc_list
        );

    void get_domain_bc (
        int* grav
        );

    void des_drag_gp(
        const int* p_id, 
        const amrex::Real* particle_vel ,
        const amrex::Real*    fluid_vel ,
        const amrex::Real* ep_g, const amrex::Real* ro, 
        const amrex::Real* mu, const amrex::Real* f_gp, 
        const int* i, const int* j, const int* k, 
        const amrex::Real* rad, const amrex::Real* vol, 
        const amrex::Real* ro_s);

    void calc_particle_collisions(
        void* rparticles, const int* nrp,
        void* gparticles, const int* ngp,
        const int* nbor_list, const int* size_nl,
        amrex::Real* tow, amrex::Real* fc,
        amrex::Real* subdt, int* ncoll
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

    void call_usr3_des(const int* np, void* particles);

    void call_usr2_des(const int* np, void* particles);

#ifdef __cplusplus
}
#endif

#endif
