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


void MFIXParticleContainer:: GetIntData ( MFIXParIter& pti, int dataIdx, int** dataPtr)
{
    BL_ASSERT( dataIdx < intData::count );
    
    *dataPtr = pti.GetStructOfArrays().GetIntData(dataIdx).dataPtr();

}


void MFIXParticleContainer:: GetIntData( const int& lev, const MFIter& mfi, int dataIdx, int** dataPtr)
{ 

    BL_ASSERT( dataIdx < intData::count );

    const int gridIndex = mfi.index();
    const int tileIndex = mfi.LocalTileIndex();
    auto&     particles = GetParticles(lev)[std::make_pair(gridIndex,tileIndex)];

    *dataPtr = particles.GetStructOfArrays().GetIntData(dataIdx).dataPtr(); 
}


void MFIXParticleContainer:: GetRealData ( MFIXParIter& pti, int dataIdx, Real** dataPtr)
{
    BL_ASSERT( dataIdx < realData::count );
    
    *dataPtr = pti.GetStructOfArrays().GetRealData(dataIdx).dataPtr();

}


void MFIXParticleContainer:: GetRealData( const int& lev, const MFIter& mfi, int dataIdx, Real** dataPtr)
{ 

    BL_ASSERT( dataIdx < realData::count );

    const int gridIndex = mfi.index();
    const int tileIndex = mfi.LocalTileIndex();
    auto&     particles = GetParticles(lev)[std::make_pair(gridIndex,tileIndex)];

    *dataPtr = particles.GetStructOfArrays().GetRealData(dataIdx).dataPtr(); 
}


void MFIXParticleContainer::GetPosition ( MFIXParIter& pti, Real** dataPtr)
{

    auto& structs = pti.GetArrayOfStructs();
    const int np  = structs.size();

    pos.resize(3*np);

    for (int i = 0; i < np; ++i ) {
	pos[i]      = structs[i].pos(0);
	pos[i+np]   = structs[i].pos(1);
	pos[i+2*np] = structs[i].pos(2);
    }

    *dataPtr = pos.dataPtr();
}



void MFIXParticleContainer::GetPosition ( const int& lev, const MFIter& mfi, Real** dataPtr)
{
    const int gridIndex = mfi.index();
    const int tileIndex = mfi.LocalTileIndex();

    auto& structs = GetParticles(lev)[std::make_pair(gridIndex,tileIndex)].GetArrayOfStructs();
    const int np  = structs.size();

    pos.resize(3*np);

    for (int i = 0; i < np; ++i ) {
	pos[i]      = structs[i].pos(0);
	pos[i+np]   = structs[i].pos(1);
	pos[i+2*np] = structs[i].pos(2);
    }

    *dataPtr = pos.dataPtr();
}


void MFIXParticleContainer::RestorePosition ( const int& lev, const MFIter& mfi, Real** dataPtr)
{
    const int gridIndex = mfi.index();
    const int tileIndex = mfi.LocalTileIndex();

    auto& structs = GetParticles(lev)[std::make_pair(gridIndex,tileIndex)].GetArrayOfStructs();
    const int np  = structs.size();

    for (int i = 0; i < np; ++i ) {
	structs[i].pos(0) = pos[i];
	structs[i].pos(1) = pos[i+np];
	structs[i].pos(2) = pos[i+2*np];
    }

    *dataPtr = NULL;
}


void MFIXParticleContainer::RestorePosition ( MFIXParIter& pti, Real** dataPtr)
{

    auto& structs = pti.GetArrayOfStructs();
    const int np  = structs.size();

    for (int i = 0; i < np; ++i ) {
	structs[i].pos(0) = pos[i];
	structs[i].pos(1) = pos[i+np];
	structs[i].pos(2) = pos[i+2*np];
    }

    *dataPtr = NULL;
}


void MFIXParticleContainer::GetVectorData ( MFIXParIter& pti, int firstCompIdx, Real** dataPtr)
{
     BL_ASSERT( firstCompIdx < realData::count );

     auto& attribs = pti.GetStructOfArrays();    
   
     if ( firstCompIdx == realData::velx ) {
	 Pack3DArrays( vel, attribs.GetRealData(realData::velx), 
		       attribs.GetRealData(realData::vely),
		       attribs.GetRealData(realData::velz) );
	 *dataPtr = vel.dataPtr();	
     }
     else if ( firstCompIdx == realData::omegax ) {
	 Pack3DArrays( omega, attribs.GetRealData(realData::omegax), 
		       attribs.GetRealData(realData::omegay),
		       attribs.GetRealData(realData::omegaz) );
	 *dataPtr = omega.dataPtr();	
     }
     else if ( firstCompIdx == realData::accx ) {
	 Pack3DArrays( acc, attribs.GetRealData(realData::accx), 
		       attribs.GetRealData(realData::accy),
		       attribs.GetRealData(realData::accz) );
	 *dataPtr = acc.dataPtr();	
     }     
     else if ( firstCompIdx == realData::alphax ) {
	 Pack3DArrays( alpha, attribs.GetRealData(realData::alphax), 
		       attribs.GetRealData(realData::alphay),
		       attribs.GetRealData(realData::alphaz) );
	 *dataPtr = alpha.dataPtr();	
     }
     else if ( firstCompIdx == realData::dragx ) {
	 Pack3DArrays( drag, attribs.GetRealData(realData::dragx), 
		       attribs.GetRealData(realData::dragy),
		       attribs.GetRealData(realData::dragz) );
	 *dataPtr = drag.dataPtr();	
     }
     else
     {
	 *dataPtr = NULL;
     }
}

void MFIXParticleContainer::GetVectorData (  const int& lev, const MFIter& mfi, int firstCompIdx, Real** dataPtr)
{
     BL_ASSERT( firstCompIdx < realData::count );


    const int gridIndex = mfi.index();
    const int tileIndex = mfi.LocalTileIndex();

    auto& attribs = GetParticles(lev)[std::make_pair(gridIndex,tileIndex)].GetStructOfArrays(); 

   
     if ( firstCompIdx == realData::velx ) {
	 Pack3DArrays( vel, attribs.GetRealData(realData::velx), 
		       attribs.GetRealData(realData::vely),
		       attribs.GetRealData(realData::velz) );
	 *dataPtr = vel.dataPtr();	
     }
     else if ( firstCompIdx == realData::omegax ) {
	 Pack3DArrays( omega, attribs.GetRealData(realData::omegax), 
		       attribs.GetRealData(realData::omegay),
		       attribs.GetRealData(realData::omegaz) );
	 *dataPtr = omega.dataPtr();	
     }
     else if ( firstCompIdx == realData::accx ) {
	 Pack3DArrays( acc, attribs.GetRealData(realData::accx), 
		       attribs.GetRealData(realData::accy),
		       attribs.GetRealData(realData::accz) );
	 *dataPtr = acc.dataPtr();	
     }     
     else if ( firstCompIdx == realData::alphax ) {
	 Pack3DArrays( alpha, attribs.GetRealData(realData::alphax), 
		       attribs.GetRealData(realData::alphay),
		       attribs.GetRealData(realData::alphaz) );
	 *dataPtr = alpha.dataPtr();	
     }
     else if ( firstCompIdx == realData::dragx ) {
	 Pack3DArrays( drag, attribs.GetRealData(realData::dragx), 
		       attribs.GetRealData(realData::dragy),
		       attribs.GetRealData(realData::dragz) );
	 *dataPtr = drag.dataPtr();	
     }
     else
     {
	 *dataPtr = NULL;
     }
}


void MFIXParticleContainer::RestoreVectorData ( MFIXParIter& pti, int firstCompIdx, Real** dataPtr)
{
     BL_ASSERT( firstCompIdx < realData::count );

     auto& attribs = pti.GetStructOfArrays();    
   
     if ( firstCompIdx == realData::velx ) {
	 Unpack3DArrays( attribs.GetRealData(realData::velx), 
			 attribs.GetRealData(realData::vely),
			 attribs.GetRealData(realData::velz), vel );
	 *dataPtr = NULL;	
     }
     else if ( firstCompIdx == realData::omegax ) {
	 Unpack3DArrays( attribs.GetRealData(realData::omegax), 
			 attribs.GetRealData(realData::omegay),
			 attribs.GetRealData(realData::omegaz), omega );
	 *dataPtr = NULL;	
     }
     else if ( firstCompIdx == realData::accx ) {
	 Unpack3DArrays( attribs.GetRealData(realData::accx), 
			 attribs.GetRealData(realData::accy),
			 attribs.GetRealData(realData::accz), acc );
	 *dataPtr = NULL;	
     }
     else if ( firstCompIdx == realData::alphax ) {
	 Unpack3DArrays( attribs.GetRealData(realData::alphax), 
			 attribs.GetRealData(realData::alphay),
			 attribs.GetRealData(realData::alphaz), alpha );
	 *dataPtr = NULL;
     }
     else if ( firstCompIdx == realData::dragx ) {
	 Unpack3DArrays( attribs.GetRealData(realData::dragx), 
			 attribs.GetRealData(realData::dragy),
			 attribs.GetRealData(realData::dragz), drag );
	 *dataPtr = NULL;	
     }
     else
     {
	 *dataPtr = NULL;
     }
}


void MFIXParticleContainer::RestoreVectorData (  const int& lev, const MFIter& mfi, int firstCompIdx, Real** dataPtr)
{
     BL_ASSERT( firstCompIdx < realData::count );


    const int gridIndex = mfi.index();
    const int tileIndex = mfi.LocalTileIndex();

    auto& attribs = GetParticles(lev)[std::make_pair(gridIndex,tileIndex)].GetStructOfArrays(); 

   
     if ( firstCompIdx == realData::velx ) {
	 Unpack3DArrays( attribs.GetRealData(realData::velx), 
			 attribs.GetRealData(realData::vely),
			 attribs.GetRealData(realData::velz), vel );
	 *dataPtr = NULL;	
     }
     else if ( firstCompIdx == realData::omegax ) {
	 Unpack3DArrays( attribs.GetRealData(realData::omegax), 
			 attribs.GetRealData(realData::omegay),
			 attribs.GetRealData(realData::omegaz), omega );
	 *dataPtr = NULL;	
     }
     else if ( firstCompIdx == realData::accx ) {
	 Unpack3DArrays( attribs.GetRealData(realData::accx), 
			 attribs.GetRealData(realData::accy),
			 attribs.GetRealData(realData::accz), acc );
	 *dataPtr = NULL;	
     }
     else if ( firstCompIdx == realData::alphax ) {
	 Unpack3DArrays( attribs.GetRealData(realData::alphax), 
			 attribs.GetRealData(realData::alphay),
			 attribs.GetRealData(realData::alphaz), alpha );
	 *dataPtr = NULL;
     }
     else if ( firstCompIdx == realData::dragx ) {
	 Unpack3DArrays( attribs.GetRealData(realData::dragx), 
			 attribs.GetRealData(realData::dragy),
			 attribs.GetRealData(realData::dragz), drag );
	 *dataPtr = NULL;	
     }
     else
     {
	 *dataPtr = NULL;
     }
}
