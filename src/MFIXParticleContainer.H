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

    const int NumberOfParticles(int& lev, MFIter& mfi)
	{
	    const int gridIndex = mfi.index();
	    const int tileIndex = mfi.LocalTileIndex();
	    auto&     particles = GetParticles(lev)[std::make_pair(gridIndex,tileIndex)];
	    auto&     structs   = particles.GetArrayOfStructs();
	    auto&     attribs   = particles.GetStructOfArrays();
	    const int np        = structs.size();

	    return structs.size();
	} 

    MFIXParticleContainer (AmrCore* amr_core);

    virtual ~MFIXParticleContainer() {};

    void InitData ();

    void InitParticlesAscii (const std::string& file);

    void printParticles ();

    void AllocData ();
   

    void GetParticlesAttributes ( MFIXParIter& pti,
				  int** pstate, int** pphase,
				  Real** pradius,  Real** pdensity,
				  Real** pvol, Real** pmass,
				  Real** pomoi, Real** ppos,  Real** pvel,
				  Real** pomega, Real** pacc,
				  Real** palpha, Real** pdrag ); 

    void GetParticlesAttributes ( const int& lev, const MFIter& mfi,
				  int** pstate, int** pphase,
				  Real** pradius,  Real** pdensity,
				  Real** pvol, Real** pmass,
				  Real** pomoi, Real** ppos,  Real** pvel,
				  Real** pomega, Real** pacc,
				  Real** palpha, Real** pdrag );


    void RestoreParticlesAttributes ( MFIXParIter& pti,
				      int** pstate, int** pphase,
				      Real** pradius,  Real** pdensity,
				      Real** pvol, Real** pmass,
				      Real** pomoi, Real** ppos,  Real** pvel,
				      Real** pomega, Real** pacc,
				      Real** palpha, Real** pdrag ); 
    
    void RestoreParticlesAttributes ( const int& lev, const MFIter& mfi,
				      int** pstate, int** pphase,
				      Real** pradius,  Real** pdensity,
				      Real** pvol, Real** pmass,
				      Real** pomoi, Real** ppos,  Real** pvel,
				      Real** pomega, Real** pacc,
				      Real** palpha, Real** pdrag );


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
