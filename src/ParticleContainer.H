
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


class MFIXParIter
    : public amrex::ParIter<0,0,PIdx::nattribs>
{
public:
    using amrex::ParIter<0,0,PIdx::nattribs>::ParIter;

    const SoA& GetAttribs () const { return GetStructOfArrays(); }
    SoA&       GetAttribs ()       { return GetStructOfArrays(); }
    const amrex::Array<amrex::Real>& GetAttribs (int comp) const { return GetStructOfArrays()[comp]; }
    amrex::Array      <amrex::Real>& GetAttribs (int comp)       { return GetStructOfArrays()[comp]; }
};



class MFIXParticleContainer
    : public amrex::ParticleContainer<0,0,PIdx::nattribs> 
{

public:

    MFIXParticleContainer (AmrCore* amr_core);

    virtual ~MFIXParticleContainer() {}

    void InitData ();

    void Evolve (int lev, Real dt);

private:

    static void ReadStaticParameters ();

    // runtime parameters
    static int do_tiling;
    static IntVect tile_size;
};
#endif /*_ParticleContainer_H_*/
