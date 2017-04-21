#include <AMReX.H>
#include "AMReX_Particles.H"
#include <iostream>
#include <MFIXParticleContainer.H>

#include<math.h>

#include "mfix_F.H"


using namespace amrex;
using namespace std;




void
MFIXParticleContainer::Pack3DArrays( Array<Real>& vec, const Array<Real>& comp1,
             const Array<Real>& comp2, const Array<Real>& comp3 )
{

    const int  np = comp1.size();

    vec.resize(3*np);

    for ( int i = 0; i < np; i++ )
    {
      vec[i]      = comp1[i];
      vec[i+np]   = comp2[i];
      vec[i+2*np] = comp3[i];
    }

}


void
MFIXParticleContainer::Unpack3DArrays( Array<Real>& comp1,  Array<Real>& comp2, Array<Real>& comp3,
				       const Array<Real>& vec )
{

    const int  np = comp1.size();

    for ( int i = 0; i < np; i++ )
    {
	comp1[i] = vec[i];
	comp2[i] = vec[i+np];
	comp3[i] = vec[i+2*np];
    }

}



void MFIXParticleContainer:: GetParticlesAttributes ( MFIXParIter& pti,
						      int** pstate, int** pphase,
						      Real** pradius,  Real** pdensity,
						      Real** pvol, Real** pmass,
						      Real** pomoi,  Real** ppos, Real** pvel,
						      Real** pomega, Real** pacc,
						      Real** palpha, Real** pdrag ) 
{ 


    auto& structs = pti.GetArrayOfStructs();
    auto& attribs = pti.GetStructOfArrays();

    if (pstate)
	*pstate = attribs.GetIntData(intData::state).dataPtr();
    if (pphase)
	*pphase = attribs.GetIntData(intData::phase).dataPtr();
    if (pradius)
	*pradius = attribs.GetRealData(realData::radius).dataPtr();
    if (pdensity)
	*pdensity = attribs.GetRealData(realData::density).dataPtr();
    if (pvol)
	*pvol = attribs.GetRealData(realData::volume).dataPtr();
    if (pmass)
	*pmass = attribs.GetRealData(realData::mass).dataPtr();
    if (pomoi)
	*pomoi = attribs.GetRealData(realData::oneOverI).dataPtr();

    if (ppos) {
	const int np  = structs.size();
	pos.resize(3*np);
	for (int i = 0; i < np; ++i ) {
	    pos[i]      = structs[i].pos(0);
	    pos[i+np]   = structs[i].pos(1);
	    pos[i+2*np] = structs[i].pos(2);
	}
	*ppos = pos.dataPtr();	
    }

    if (pvel) {
	Pack3DArrays( vel, attribs.GetRealData(realData::velx), 
		      attribs.GetRealData(realData::vely),
		      attribs.GetRealData(realData::velz) );
	*pvel = vel.dataPtr();	
    }

    if (pomega) {
	Pack3DArrays( omega, attribs.GetRealData(realData::omegax), 
		      attribs.GetRealData(realData::omegay),
		      attribs.GetRealData(realData::omegaz) );
	*pomega = omega.dataPtr();	
    }

    if (pacc) {
	Pack3DArrays( acc, attribs.GetRealData(realData::accx), 
		      attribs.GetRealData(realData::accy),
		      attribs.GetRealData(realData::accz) );
	*pacc = acc.dataPtr();	
    }

    if (palpha) {
	Pack3DArrays( alpha, attribs.GetRealData(realData::alphax), 
		      attribs.GetRealData(realData::alphay),
		      attribs.GetRealData(realData::alphaz) );
	*palpha= alpha.dataPtr();	
    }

    if (pdrag) {
	Pack3DArrays( drag, attribs.GetRealData(realData::dragx), 
		      attribs.GetRealData(realData::dragy),
		      attribs.GetRealData(realData::dragz) );
	*pdrag= drag.dataPtr();	
    }

}


void MFIXParticleContainer::GetParticlesAttributes ( const int& lev, const MFIter& mfi,
						     int** pstate, int** pphase,
						     Real** pradius,  Real** pdensity,
						     Real** pvol, Real** pmass,
						     Real** pomoi, Real** ppos,  Real** pvel,
						     Real** pomega, Real** pacc,
						     Real** palpha, Real** pdrag )
{ 

    const int gridIndex = mfi.index();
    const int tileIndex = mfi.LocalTileIndex();
    auto&     particles = GetParticles(lev)[std::make_pair(gridIndex,tileIndex)];
    auto&     structs   = particles.GetArrayOfStructs();
    auto&     attribs   = particles.GetStructOfArrays();

    if (pstate)
	*pstate = attribs.GetIntData(intData::state).dataPtr();
    if (pphase)
	*pphase = attribs.GetIntData(intData::phase).dataPtr();
    if (pradius)
	*pradius = attribs.GetRealData(realData::radius).dataPtr();
    if (pdensity)
	*pdensity = attribs.GetRealData(realData::density).dataPtr();
    if (pvol)
	*pvol = attribs.GetRealData(realData::volume).dataPtr();
    if (pmass)
	*pmass = attribs.GetRealData(realData::mass).dataPtr();
    if (pomoi)
	*pomoi = attribs.GetRealData(realData::oneOverI).dataPtr();

    if (ppos) {
	const int np  = structs.size();
	pos.resize(3*np);
	for (int i = 0; i < np; ++i ) {
	    pos[i]      = structs[i].pos(0);
	    pos[i+np]   = structs[i].pos(1);
	    pos[i+2*np] = structs[i].pos(2);
	}
	*ppos = pos.dataPtr();	
    }

    if (pvel) {
	Pack3DArrays( vel, attribs.GetRealData(realData::velx), 
		      attribs.GetRealData(realData::vely),
		      attribs.GetRealData(realData::velz) );
	*pvel = vel.dataPtr();	
    }

    if (pomega) {
	Pack3DArrays( omega, attribs.GetRealData(realData::omegax), 
		      attribs.GetRealData(realData::omegay),
		      attribs.GetRealData(realData::omegaz) );
	*pomega = omega.dataPtr();	
    }

    if (pacc) {
	Pack3DArrays( acc, attribs.GetRealData(realData::accx), 
		      attribs.GetRealData(realData::accy),
		      attribs.GetRealData(realData::accz) );
	*pacc = acc.dataPtr();	
    }

    if (palpha) {
	Pack3DArrays( alpha, attribs.GetRealData(realData::alphax), 
		      attribs.GetRealData(realData::alphay),
		      attribs.GetRealData(realData::alphaz) );
	*palpha= alpha.dataPtr();	
    }

    if (pdrag) {
	Pack3DArrays( drag, attribs.GetRealData(realData::dragx), 
		      attribs.GetRealData(realData::dragy),
		      attribs.GetRealData(realData::dragz) );
	*pdrag= drag.dataPtr();	
    }
}




void MFIXParticleContainer:: RestoreParticlesAttributes ( MFIXParIter& pti,
							  int** pstate, int** pphase,
							  Real** pradius,  Real** pdensity,
							  Real** pvol, Real** pmass,
							  Real** pomoi,  Real** ppos, Real** pvel,
							  Real** pomega, Real** pacc,
							  Real** palpha, Real** pdrag  ) 
{ 

    auto& structs = pti.GetArrayOfStructs();
    auto& attribs = pti.GetStructOfArrays();
    const int np  = structs.size();

    if (pstate)
	*pstate = NULL;
    if (pphase)
	*pphase = NULL;
    if (pradius)
	*pradius = NULL;
    if (pdensity)
	*pdensity = NULL;
    if (pvol)
	*pvol = NULL;
    if (pmass)
	*pmass = NULL;
    if (pomoi)
	*pomoi = NULL;

    if (ppos) {
	const int np  = structs.size();
	for (int i = 0; i < np; ++i ) {
	    structs[i].pos(0) = pos[i];
	    structs[i].pos(1) = pos[i+np];
	    structs[i].pos(2) = pos[i+2*np];
	}
	*ppos = NULL;
    }

    if (pvel) {
	Unpack3DArrays( attribs.GetRealData(realData::velx), 
		      attribs.GetRealData(realData::vely),
			attribs.GetRealData(realData::velz), vel );
	*pvel = NULL;	
    }

    if (pomega) {
	Unpack3DArrays( attribs.GetRealData(realData::omegax), 
		      attribs.GetRealData(realData::omegay),
			attribs.GetRealData(realData::omegaz), omega );
	*pomega = NULL;	
    }

    if (pacc) {
	Unpack3DArrays( attribs.GetRealData(realData::accx), 
		      attribs.GetRealData(realData::accy),
			attribs.GetRealData(realData::accz), acc );
	*pacc = NULL;	
    }

    if (palpha) {
	Unpack3DArrays( attribs.GetRealData(realData::alphax), 
		      attribs.GetRealData(realData::alphay),
			attribs.GetRealData(realData::alphaz), alpha );
	*palpha= NULL;	
    }

    if (pdrag) {
	Unpack3DArrays( attribs.GetRealData(realData::dragx), 
		      attribs.GetRealData(realData::dragy),
			attribs.GetRealData(realData::dragz), drag );
	*pdrag= NULL;	
    }

}


void MFIXParticleContainer:: RestoreParticlesAttributes ( const int& lev, const MFIter& mfi,
							  int** pstate, int** pphase,
							  Real** pradius,  Real** pdensity,
							  Real** pvol, Real** pmass,
							  Real** pomoi,  Real** ppos, Real** pvel,
							  Real** pomega, Real** pacc,
							  Real** palpha, Real** pdrag ) 
{ 


    const int gridIndex = mfi.index();
    const int tileIndex = mfi.LocalTileIndex();
    auto&     particles = GetParticles(lev)[std::make_pair(gridIndex,tileIndex)];
    auto&     structs   = particles.GetArrayOfStructs();
    auto&     attribs   = particles.GetStructOfArrays();
    const int np        = structs.size();

    if (pstate)
	*pstate = NULL;
    if (pphase)
	*pphase = NULL;
    if (pradius)
	*pradius = NULL;
    if (pdensity)
	*pdensity = NULL;
    if (pvol)
	*pvol = NULL;
    if (pmass)
	*pmass = NULL;
    if (pomoi)
	*pomoi = NULL;

    if (ppos) {
	const int np  = structs.size();
	for (int i = 0; i < np; ++i ) {
	    structs[i].pos(0) = pos[i];
	    structs[i].pos(1) = pos[i+np];
	    structs[i].pos(2) = pos[i+2*np];
	}
	*ppos = NULL; 
    }

    if (pvel) {
	Unpack3DArrays( attribs.GetRealData(realData::velx), 
		      attribs.GetRealData(realData::vely),
			attribs.GetRealData(realData::velz), vel );
	*pvel = NULL;	
    }

    if (pomega) {
	Unpack3DArrays( attribs.GetRealData(realData::omegax), 
		      attribs.GetRealData(realData::omegay),
			attribs.GetRealData(realData::omegaz), omega );
	*pomega = NULL;	
    }

    if (pacc) {
	Unpack3DArrays( attribs.GetRealData(realData::accx), 
		      attribs.GetRealData(realData::accy),
			attribs.GetRealData(realData::accz), acc );
	*pacc = NULL;	
    }

    if (palpha) {
	Unpack3DArrays( attribs.GetRealData(realData::alphax), 
		      attribs.GetRealData(realData::alphay),
			attribs.GetRealData(realData::alphaz), alpha );
	*palpha= NULL;	
    }

    if (pdrag) {
	Unpack3DArrays( attribs.GetRealData(realData::dragx), 
		      attribs.GetRealData(realData::dragy),
			attribs.GetRealData(realData::dragz), drag );
	*pdrag= NULL;	
    }


}




