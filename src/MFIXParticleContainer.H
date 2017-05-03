#ifndef _MFIXParticleContainer_H_
#define _MFIXParticleContainer_H_


#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include "AMReX_Particles.H"


using namespace amrex;
using namespace std;

struct realData
{
   enum { // Particle data stored in p.m_data
      radius = 0,
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
      count

   };
};

struct intData
{
   enum { // Particle data stored in p.m_data
      phase = 0,
      state,
      count
   };
};



class MFIXParIter
   : public amrex::ParIter<0,0,realData::count,intData::count>
{
public:
   using amrex::ParIter<0,0,realData::count,intData::count>::ParIter;
};

class MFIXParticleContainer
   : public amrex::ParticleContainer<0,0,realData::count,intData::count>
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
	 const int np        = structs.size();

	 return structs.size();
      }

   MFIXParticleContainer (AmrCore* amr_core);

   virtual ~MFIXParticleContainer() {};

   void InitData ();

   void InitParticlesAscii (const std::string& file);

   void printParticles ();

   void AllocData ();

   void GetIntData ( MFIXParIter& pti, int dataIdx, int** dataPtr );

   void GetIntData ( const int& lev, const MFIter& mfi, int dataIdx, int** dataPtr );

   void RestoreIntData ( int** dataPtr ) { *dataPtr = NULL;};

   void GetRealData ( MFIXParIter& pti, int dataIdx, Real** dataPtr );

   void GetRealData ( const int& lev, const MFIter& mfi, int dataIdx, Real** dataPtr );

   void RestoreRealData ( Real** dataPtr ) { *dataPtr = NULL;};

   void GetPosition ( MFIXParIter& pti, Real** dataPtr );

   void GetPosition ( const int& lev, const MFIter& mfi, Real** dataPtr );

   void RestorePosition ( MFIXParIter& pti, Real** dataPtr );

   void RestorePosition ( const int& lev, const MFIter& mfi, Real** dataPtr );

   void GetVectorData ( MFIXParIter& pti, int firstCompIdx, Real** dataPtr );

   void GetVectorData ( const int& lev, const MFIter& mfi, int firstCompIdx, Real** dataPtr );

   void RestoreVectorData ( MFIXParIter& pti, int firstCompIdx, Real** dataPtr);

   void RestoreVectorData ( const int& lev, const MFIter& mfi, int firstCompIdx, Real** dataPtr );

   void EvolveParticles(Array< unique_ptr<MultiFab> >& ep_g,
			const Array< unique_ptr<MultiFab> >& u_g,
			const Array< unique_ptr<MultiFab> >& v_g,
			const Array< unique_ptr<MultiFab> >& w_g,
			const Array< unique_ptr<MultiFab> >& p_g,
			const Array< unique_ptr<MultiFab> >& ro_g,
			const Array< unique_ptr<MultiFab> >& mu_g,
			int lev, int nstep, Real dt, Real time);

   void output(int lev, int estatus, int finish, int nstep, Real dt, Real time);

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel(int lev = 0);

private:

   void Pack3DArrays( Array<Real>& vec, const Array<Real>& comp1,
		      const Array<Real>& comp2, const Array<Real>& comp3 );

   void Unpack3DArrays( Array<Real>& comp1,  Array<Real>& comp2, Array<Real>& comp3,
			const Array<Real>& vec );

   static void ReadStaticParameters ();



   // runtime parameters
   static int do_tiling;
   static IntVect tile_size;

   // These are working arrays used to pack 3D components into a 3D vector
   Array<Real>      pos;
   Array<Real>      vel;
   Array<Real>      acc;
   Array<Real>    omega;
   Array<Real>    alpha;
   Array<Real>     drag;


};
#endif
