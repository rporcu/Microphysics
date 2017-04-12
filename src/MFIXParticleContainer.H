#ifndef _MFIXParticleContainer_H_
#define _MFIXParticleContainer_H_


#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include "AMReX_Particles.H"


using namespace amrex;


struct PIdx
{
    enum { // Particle data stored in p.m_data
	phase = 0, 
	state,
        radius,
	volume,
	mass,
        density, 
        oneOverI,
	velx,
	vely,
	velz,
        omegax,
	omegay,
	omegaz,
	accx,
	accy,
	accz,
        alphax,  //angular acceleration
	alphay, 
	alphaz,
        dragx,
	dragy,
	dragz,
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

    virtual ~MFIXParticleContainer() {};

    void InitData ();

    void Evolve (int lev, Real dt);

    void InitParticlesAscii (const std::string& file);

    void printParticles ();

    void AllocData ();

    const int  GetNParticles () { return numberOfParticles; };

    void GetParticlesPosition ( Array<Real>& des_pos_new );

    void GetParticlesAttributes (
	Array<int>& particle_state, Array<int>& particle_phase,
	Array<Real>& des_radius,  Array<Real>& ro_sol,
	Array<Real>& pvol,   Array<Real>& pmass,
	Array<Real>& omoi,   Array<Real>& des_vel_new,
	Array<Real>& omega_new,   Array<Real>& des_acc_old,
	Array<Real>& rot_acc_old,   Array<Real>& drag_fc);    

private:

    void Pack3DArrays( Array<Real>& vec, Array<Real>& comp1,
		       Array<Real>& comp2, Array<Real>& comp3 );
    
    static void ReadStaticParameters ();

    int  numberOfParticles = 0; 
    
    // runtime parameters
    static int do_tiling;
    static IntVect tile_size;
};
#endif 
