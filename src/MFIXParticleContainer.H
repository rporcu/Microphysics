#ifndef _MFIXParticleContainer_H_
#define _MFIXParticleContainer_H_


#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include "AMReX_Particles.H"


using namespace amrex;
using namespace std;


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


    ///
    /// Define some convenient wrappers for accessing particle data
    ///
    ParticleType::RealType& x(int i) { return GetArrayOfStructs()[i].pos(0); }
    ParticleType::RealType& y(int i) { return GetArrayOfStructs()[i].pos(1); }
    ParticleType::RealType& z(int i) { return GetArrayOfStructs()[i].pos(2); }

    ParticleType::RealType& velx(int i) { return GetStructOfArrays()[PIdx::velx][i]; }
    ParticleType::RealType& vely(int i) { return GetStructOfArrays()[PIdx::vely][i]; }
    ParticleType::RealType& velz(int i) { return GetStructOfArrays()[PIdx::velz][i]; }

    int phase(int i) { return int(GetStructOfArrays()[PIdx::phase][i]); }
    int state(int i) { return int(GetStructOfArrays()[PIdx::state][i]); }


    // Getter and Setter for attributes
    const ParticleType::RealType& vel();
    const ParticleType::RealType& vel(int grid, int tile);


};








class MFIXParticleContainer
    : public amrex::ParticleContainer<0,0,PIdx::nattribs> 
{

public:

    // Get number of particles

    const int NumberOfParticles(MFIXParIter& pti) { return pti.GetArrayOfStructs().size(); } 

    MFIXParticleContainer (AmrCore* amr_core);

    virtual ~MFIXParticleContainer() {};

    void InitData ();

    void InitParticlesAscii (const std::string& file);

    void printParticles ();

    void AllocData ();

    const int  GetNParticles () { return numberOfParticles; };

    void GetParticlesPosition ( Array<Real>& des_pos_new );
    void GetPositionsArray( Real** posPtr, MFIXParIter& pti );
    void RestorePositionsArray( Real** posPtr, MFIXParIter& pti );
    void GetVelocitiesArray( Real** velPtr, MFIXParIter& pti );
    void RestoreVelocitiesArray( Real** velPtr, MFIXParIter& pti );

    void GetParticlesAttributes (
	Array<int>& particle_state, Array<int>& particle_phase,
	Array<Real>& des_radius,  Array<Real>& ro_sol,
	Array<Real>& pvol,   Array<Real>& pmass,
	Array<Real>& omoi,   Array<Real>& des_vel_new,
	Array<Real>& omega_new,   Array<Real>& des_acc_old,
	Array<Real>& rot_acc_old,   Array<Real>& drag_fc);    

    void GetParticlesAttributes ( int** pstate, int** pphase,
				  Real** pradius,  Real** pdensity,
				  Real** pvol, Real** pmass,
				  Real** pomoi, Real** ppos,  Real** pvel,
				  Real** pomega, Real** pacc,
				  Real** palpha, Real** pdrag,
				  MFIXParIter& pti); 

    void RestoreParticlesAttributes ( int** pstate, int** pphase,
				  Real** pradius,  Real** pdensity,
				  Real** pvol, Real** pmass,
				  Real** pomoi, Real** ppos,  Real** pvel,
				  Real** pomega, Real** pacc,
				  Real** palpha, Real** pdrag,
				  MFIXParIter& pti); 

    void EvolveParticles(Array< unique_ptr<MultiFab> >& ep_g,
			 const Array< unique_ptr<MultiFab> >& u_g,
			 const Array< unique_ptr<MultiFab> >& v_g,
			 const Array< unique_ptr<MultiFab> >& w_g,
			 const Array< unique_ptr<MultiFab> >& p_g,
			 const Array< unique_ptr<MultiFab> >& ro_g,
			 const Array< unique_ptr<MultiFab> >& mu_g,
			 int lev, int nstep, Real dt, Real time);

    void output(int lev, int estatus, int finish, int nstep, Real dt, Real time);

private:

    void Pack3DArrays( Array<Real>& vec, const Array<Real>& comp1,
		       const Array<Real>& comp2, const Array<Real>& comp3 );

    void Unpack3DArrays( Array<Real>& comp1,  Array<Real>& comp2, Array<Real>& comp3,
			 const Array<Real>& vec );
    
    static void ReadStaticParameters ();

    int  numberOfParticles = 0; 
    
    // runtime parameters
    static int do_tiling;
    static IntVect tile_size;

    // These are working arrays used to convert from Real to int to conform to AMReX 
    // way of doing things
    Array<int>     phase;
    Array<int>     state;

    // These are working arrays used to pack 3D components into a 3D vector
    Array<Real>      pos;
    Array<Real>      vel;
    Array<Real>      acc;
    Array<Real>    omega;
    Array<Real>    alpha;
    Array<Real>     drag;


};
#endif 
