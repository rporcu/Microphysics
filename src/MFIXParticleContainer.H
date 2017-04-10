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
        cforcex,  //Collision force for linear momenta
	cforcey,
	cforcez,
        acforcex,  //Collision force for angular momenta
	acforcey,
	acforcez,
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

    virtual ~MFIXParticleContainer() {};

    void InitData ();

    void Evolve (int lev, Real dt);

    void InitParticlesAscii (const std::string& file);

    void printParticles ();

    void AllocData ();

    const int  GetNParticles () { return numberOfParticles; };

    void GetParticlesPosition ( Array<Real>& des_pos_new );

    

private:

    static void ReadStaticParameters ();

    int  numberOfParticles = 0; 
    
    // runtime parameters
    static int do_tiling;
    static IntVect tile_size;
};
#endif 
