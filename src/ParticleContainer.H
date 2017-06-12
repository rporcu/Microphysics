#ifndef _ParticleContainer_H_
#define _ParticleContainer_H_

#include <memory>
#include <map>

#include <AMReX_Particles.H>
#include <AMReX_AmrParGDB.H>

using namespace amrex;



struct PIdx
{
    enum { // Particle data stored in p.m_data
  particle_phase = 0,
  particle_state,
        des_radius,
        ro_sol,
        pvol, pmass, omoi,
        des_pos_new,
        des_vel_new,
        des_usr_var,
        omega_new,
        des_acc_old,
        des_rot_old,
        drag_fc,
        fc,
        tow,
        pairs,
  nattribs
    };
};

class MyParticleContainer
    : public amrex::ParticleContainer<PIdx::nattribs,0> //std::vector<Particle<PIdx::nattribs,0> > >
{

public:

    MyParticleContainer (AmrCore* amr_core);

    virtual ~MyParticleContainer() {}

    void InitData ();

    void Evolve (int lev, Real dt);

private:

    static void ReadStaticParameters ();
};
#endif /*_ParticleContainer_H_*/
