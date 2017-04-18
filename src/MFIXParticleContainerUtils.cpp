#include <AMReX.H>
#include "AMReX_Particles.H"
#include <iostream>
#include <MFIXParticleContainer.H>

#include<math.h>

#include "mfix_F.H"


using namespace amrex;
using namespace std;







void
MFIXParticleContainer::GetVelocitiesArray( Real** velPtr, MFIXParIter& pti ) {


    auto& structs = pti.GetArrayOfStructs();
    const int np  = structs.size();
    auto& attribs = pti.GetStructOfArrays();

    vel.resize(3*np);

    for (int i = 0; i < np; ++i )
    {
	vel[i]      = attribs[PIdx::velx][i];
	vel[i+np]   = attribs[PIdx::vely][i];
	vel[i+2*np] = attribs[PIdx::velz][i];
    }

    *velPtr = vel.dataPtr();
}


void
MFIXParticleContainer::RestoreVelocitiesArray( Real** velPtr, MFIXParIter& pti ) {


    auto& structs = pti.GetArrayOfStructs();
    const int np  = structs.size();
    auto& attribs = pti.GetStructOfArrays();

    vel.resize(3*np);

    for (int i = 0; i < np; ++i )
    {
	attribs[PIdx::velx][i] = vel[i];
	attribs[PIdx::vely][i] = vel[i+np];
	attribs[PIdx::velz][i] = vel[i+2*np];
    }

    *velPtr = NULL;
}




void
MFIXParticleContainer::GetPositionsArray( Real** posPtr, MFIXParIter& pti ) {


    auto& structs = pti.GetArrayOfStructs();
    const int np  = structs.size();

    pos.resize(3*np);

    for (int i = 0; i < np; ++i )
    {
	pos[i]      = structs[i].pos(0);
	pos[i+np]   = structs[i].pos(1);
	pos[i+2*np] = structs[i].pos(2);
    }

    *posPtr = pos.dataPtr();
}

void
MFIXParticleContainer::RestorePositionsArray( Real** posPtr, MFIXParIter& pti ) {


    auto& structs = pti.GetArrayOfStructs();
    const int np  = structs.size();

    for (int i = 0; i < np; ++i )
    {
	structs[i].pos(0) = pos[i];
	structs[i].pos(1) = pos[i+np];
	structs[i].pos(2) = pos[i+2*np];
    }

    *posPtr = NULL;
}



void
MFIXParticleContainer:: GetParticlesPosition( Array<Real>& des_pos_new ) {



    Array<Real>  xp, yp, zp;

    xp.resize(numberOfParticles);
    yp.resize(numberOfParticles);
    zp.resize(numberOfParticles);

    int lev = 0;
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
	pti.GetPosition(xp, yp, zp);
    }

    Pack3DArrays(des_pos_new, xp, yp, zp);

}






void
MFIXParticleContainer::Pack3DArrays( Array<Real>& vec, const Array<Real>& comp1,
             const Array<Real>& comp2, const Array<Real>& comp3 )
{

    const int  np = comp1.size();

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



void MFIXParticleContainer:: GetParticlesAttributes ( int** pstate, int** pphase,
						      Real** pradius,  Real** pdensity,
						      Real** pvol, Real** pmass,
						      Real** pomoi,  Real** ppos, Real** pvel,
						      Real** pomega, Real** pacc,
						      Real** palpha, Real** pdrag,
						      MFIXParIter& pti) 
{ 


    auto& structs = pti.GetArrayOfStructs();
    auto& attribs = pti.GetStructOfArrays();
    const int np  = structs.size();

    state.resize(np);
    phase.resize(np);
    pos.resize(3*np);
    vel.resize(3*np);
    omega.resize(3*np);
    acc.resize(3*np);
    alpha.resize(3*np);
    drag.resize(3*np);

    for (int i = 0; i < np; ++i )
    {

	state[i] = int( attribs[PIdx::state][i] );
	phase[i] = int( attribs[PIdx::phase][i] );

	pos[i]      = structs[i].pos(0);
	pos[i+np]   = structs[i].pos(1);
	pos[i+2*np] = structs[i].pos(2);	

	vel[i]      = attribs[PIdx::velx][i];
	vel[i+np]   = attribs[PIdx::vely][i];
	vel[i+2*np] = attribs[PIdx::velz][i];

	omega[i]      = attribs[PIdx::omegax][i];
	omega[i+np]   = attribs[PIdx::omegay][i];
	omega[i+2*np] = attribs[PIdx::omegaz][i];

	acc[i]      = attribs[PIdx::accx][i];
	acc[i+np]   = attribs[PIdx::accy][i];
	acc[i+2*np] = attribs[PIdx::accz][i];

	alpha[i]      = attribs[PIdx::alphax][i];
	alpha[i+np]   = attribs[PIdx::alphay][i];
	alpha[i+2*np] = attribs[PIdx::alphaz][i];

	drag[i]      = attribs[PIdx::dragx][i];
	drag[i+np]   = attribs[PIdx::dragy][i];
	drag[i+2*np] = attribs[PIdx::dragz][i];
	
    }

    if (pstate)
	*pstate = state.dataPtr();
    if (pphase)
	*pphase = phase.dataPtr();
    if (pradius)
	*pradius = attribs[PIdx::radius].dataPtr();
    if (pdensity)
	*pdensity = attribs[PIdx::density].dataPtr();
    if (pvol)
	*pvol = attribs[PIdx::volume].dataPtr();
    if (pmass)
	*pmass = attribs[PIdx::mass].dataPtr();
    if (pomoi)
	*pomoi = attribs[PIdx::oneOverI].dataPtr();
    if (ppos)
	*ppos = pos.dataPtr();
    if (pvel)
	*pvel = vel.dataPtr();
    if (pomega)
	*pomega = omega.dataPtr();
    if (pacc)
	*pacc = acc.dataPtr();
    if (palpha)
	*palpha = alpha.dataPtr();
    if (pdrag)
	*pdrag = drag.dataPtr();
}


void MFIXParticleContainer:: RestoreParticlesAttributes ( int** pstate, int** pphase,
							  Real** pradius,  Real** pdensity,
							  Real** pvol, Real** pmass,
							  Real** pomoi,  Real** ppos, Real** pvel,
							  Real** pomega, Real** pacc,
							  Real** palpha, Real** pdrag,
							  MFIXParIter& pti ) 
{ 


    auto& structs = pti.GetArrayOfStructs();
    auto& attribs = pti.GetStructOfArrays();
    const int np  = structs.size();



    for (int i = 0; i < np; ++i )
    {

	attribs[PIdx::state][i] = state[i];
	attribs[PIdx::phase][i] = phase[i];

	structs[i].pos(0) = pos[i];
	structs[i].pos(1) = pos[i+np];
	structs[i].pos(2) = pos[i+2*np];

	attribs[PIdx::velx][i] = vel[i];
	attribs[PIdx::vely][i] = vel[i+np];
	attribs[PIdx::velz][i] = vel[i+2*np];

	attribs[PIdx::omegax][i] = omega[i];
	attribs[PIdx::omegay][i] = omega[i+np];
	attribs[PIdx::omegaz][i] = omega[i+2*np];

	attribs[PIdx::accx][i] = acc[i];
	attribs[PIdx::accy][i] = acc[i+np];
	attribs[PIdx::accz][i] = acc[i+2*np];

	attribs[PIdx::alphax][i] = alpha[i];
	attribs[PIdx::alphay][i] = alpha[i+np];
	attribs[PIdx::alphaz][i] = alpha[i+2*np];

	attribs[PIdx::dragx][i] = drag[i];
	attribs[PIdx::dragx][i] = drag[i+np];
	attribs[PIdx::dragz][i] = drag[i+2*np];
	
    }


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
    if (ppos)
	*ppos = NULL;
    if (pvel)
	*pvel = NULL;
    if (pomega)
	*pomega = NULL;
    if (pacc)
	*pacc = NULL;
    if (palpha)
	*palpha = NULL;
    if (pdrag)
	*pdrag = NULL;

}






