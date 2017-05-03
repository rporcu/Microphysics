#include <AMReX.H>
#include "AMReX_Particles.H"
#include <iostream>
#include <MFIXParticleContainer.H>

#include<math.h>

#include "mfix_F.H"


using namespace amrex;
using namespace std;


int     MFIXParticleContainer::do_tiling = 0;
IntVect MFIXParticleContainer::tile_size   { D_DECL(1024000,8,8) };

MFIXParticleContainer::MFIXParticleContainer (AmrCore* amr_core)
    : ParticleContainer<0,0,realData::count,intData::count>
    (amr_core->GetParGDB())
{
    ReadStaticParameters();
  
    this->SetVerbose(0);
}

void MFIXParticleContainer::AllocData ()
{
    reserveData();
    resizeData();
}



void MFIXParticleContainer::InitParticlesAscii(const std::string& file) {

    // only read the file on the IO proc
    if (ParallelDescriptor::MyProc() ==  ParallelDescriptor::IOProcessorNumber()) {
	std::ifstream ifs;
	ifs.open(file.c_str(), std::ios::in);

	if (!ifs.good())
	    amrex::FileOpenFailed(file);

	int np = -1;
	ifs >> np >> std::ws;

	// Issue an error if nparticles = 0 is specified
	if ( np == -1 ){
	    Abort("\nCannot read number of particles from particle_input.dat: file is corrupt.\
\nPerhaps you forgot to specify the number of particles on the first line??? ");
	}

	// we add all the particles to grid 0 and tile 0 and let
	// Redistribute() put them in the right places.
	const int lev  = 0;
	const int grid = 0;
	const int tile = 0;

	ParticleType p;
	int        pstate, pphase;
	Real       pradius, pdensity, pvolume, pomoi, pmass;

	std::array<Real, realData::count> realAttribs;
	std::array<int,   intData::count>  intAttribs;

	pstate = 1;


	for (int i = 0; i < np; i++) {

	    // Init to 0 all attributes
	    realAttribs.fill(0.0);
	    intAttribs.fill(0);

	    // Read from input file
	    ifs >> pphase;
	    ifs >> p.pos(0);
	    ifs >> p.pos(1);
	    ifs >> p.pos(2);
	    ifs >> pradius;
	    ifs >> pdensity;
	    ifs >> realAttribs[realData::velx];
	    ifs >> realAttribs[realData::vely];
	    ifs >> realAttribs[realData::velz];

	    // Set id and cpu for this particle
	    p.id()  = ParticleType::NextID();
	    p.cpu() = ParallelDescriptor::MyProc();

	    // Compute other particle properties
	    mfix_set_particle_properties( &pstate, &pradius, &pdensity, &pvolume, &pmass, &pomoi);

	    // Set other particle properties
	    intAttribs[intData::phase]      = pphase;
	    intAttribs[intData::state]      = pstate;
	    realAttribs[realData::volume]   = pvolume;
	    realAttribs[realData::density]  = pdensity;
	    realAttribs[realData::mass]     = pmass;
	    realAttribs[realData::oneOverI] = pomoi;
	    realAttribs[realData::radius]   = pradius;

	    // Add everything to the data structure
	    auto& particle_tile = GetParticles(lev)[std::make_pair(grid,tile)];
	    particle_tile.push_back(p);
	    particle_tile.push_back_real(realAttribs);
	    particle_tile.push_back_int(intAttribs);
	}
    }
    Redistribute();

}


void MFIXParticleContainer:: printParticles() {
    const int lev = 0;
    const auto& plevel = GetParticles(lev);
    for (const auto& kv : plevel) {
	const auto& particle_structs = kv.second.GetArrayOfStructs();
	const auto& attribs = kv.second.GetStructOfArrays();

	for (unsigned i = 0; i < particle_structs.numParticles(); ++i) {
	    std::cout << "Particle ID  = " << i << " " << std::endl;
	    for (int j = 0; j < realData::count; j++) {
		std::cout << attribs.GetRealData(j)[i] << " " << std::endl;
	    }
	    std::cout << std::endl;
	}
    }
}


void MFIXParticleContainer::ReadStaticParameters ()
{
    static bool initialized = false;

    if (!initialized)
    {
	ParmParse pp("particles");

	pp.query("do_tiling",  do_tiling);

	Array<int> ts(BL_SPACEDIM);

	if (pp.queryarr("tile_size", ts)) {
	    tile_size = IntVect(ts);
	}

	initialized = true;
    }
}

void
MFIXParticleContainer::InitData()
{
}



// void MFIXParticleContainer::EvolveParticles(Array< unique_ptr<MultiFab> >& ep_g,
// 					    const Array< unique_ptr<MultiFab> >& u_g,
// 					    const Array< unique_ptr<MultiFab> >& v_g,
// 					    const Array< unique_ptr<MultiFab> >& w_g,
// 					    const Array< unique_ptr<MultiFab> >& p_g,
// 					    const Array< unique_ptr<MultiFab> >& ro_g,
// 					    const Array< unique_ptr<MultiFab> >& mu_g,
// 					    int lev, int nstep, Real dt, Real time )
// {

//     Box domain(Geom(lev).Domain());

//     Real dx = Geom(lev).CellSize(0);
//     Real dy = Geom(lev).CellSize(1);
//     Real dz = Geom(lev).CellSize(2);

//     Real xlen = Geom(lev).ProbHi(0) - Geom(lev).ProbLo(0);
//     Real ylen = Geom(lev).ProbHi(1) - Geom(lev).ProbLo(1);
//     Real zlen = Geom(lev).ProbHi(2) - Geom(lev).ProbLo(2);


//     for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
//     {

// 	//number of particles
// 	const int np =  NumberOfParticles(pti);

// 	// Temporary arrays
// 	Real         *pradius, *pdensity, *pvol, *pmass, *pomoi;
// 	Real         *ppos, *pvel, *pomega, *pacc, *palpha, *pdrag;
// 	int          *pstate, *pphase;
// 	int          nsubsteps;


// 	GetIntData( pti, intData::state, &pstate );
// 	GetIntData( pti, intData::phase, &pphase );

// 	GetRealData( pti, realData::radius, &pradius );
// 	GetRealData( pti, realData::volume, &pvol );
// 	GetRealData( pti, realData::mass, &pmass );
// 	GetRealData( pti, realData::oneOverI, &pomoi );

// 	GetPosition( pti, &ppos );

// 	GetVectorData( pti, realData::velx, &pvel );
// 	GetVectorData( pti, realData::omegax, &pomega );
// 	GetVectorData( pti, realData::accx, &pacc );
// 	GetVectorData( pti, realData::alphax, &palpha );
// 	GetVectorData( pti, realData::dragx, &pdrag );


// 	mfix_des_init_time_loop( &np,
// 				 pstate, pphase, pradius, pvol, ppos,  pvel, pomega,
// 				 pdrag, &time, &dt, &dx, &dy, &dz,
// 				 &xlen, &ylen, &zlen, &nstep, &nsubsteps);

// 	int quit;

// 	for ( int n = 0; n < nsubsteps; ++n ) {
// 	    mfix_des_time_loop_ops( &np,
// 				    pstate, pphase, pradius, pvol, pmass, pomoi, ppos,  pvel, pomega,
// 				    pacc,  palpha, pdrag, &time, &dt, &dx, &dy, &dz,
// 				    &xlen, &ylen, &zlen, &nstep, &quit);

// 	    if ( quit ) break;
// 	}


// 	mfix_des_finalize_time_loop( &np, &dt, ppos, pvel, pomega );


// 	RestoreIntData( &pstate );
// 	RestoreIntData( &pphase );

// 	RestoreRealData( &pradius );
// 	RestoreRealData( &pvol );
// 	RestoreRealData( &pmass );
// 	RestoreRealData( &pomoi );

// 	RestorePosition( pti, &ppos );

// 	RestoreVectorData( pti, realData::velx, &pvel );
// 	RestoreVectorData( pti, realData::omegax, &pomega );
// 	RestoreVectorData( pti, realData::accx, &pacc );
// 	RestoreVectorData( pti, realData::alphax, &palpha );
// 	RestoreVectorData( pti, realData::dragx, &pdrag );
//     }

//     Redistribute();

// }



void MFIXParticleContainer::EvolveParticles(Array< unique_ptr<MultiFab> >& ep_g,
					    const Array< unique_ptr<MultiFab> >& u_g,
					    const Array< unique_ptr<MultiFab> >& v_g,
					    const Array< unique_ptr<MultiFab> >& w_g,
					    const Array< unique_ptr<MultiFab> >& p_g,
					    const Array< unique_ptr<MultiFab> >& ro_g,
					    const Array< unique_ptr<MultiFab> >& mu_g,
					    int lev, int nstep, Real dt, Real time )
{

    Box domain(Geom(lev).Domain());

    Real dx = Geom(lev).CellSize(0);
    Real dy = Geom(lev).CellSize(1);
    Real dz = Geom(lev).CellSize(2);

    Real xlen = Geom(lev).ProbHi(0) - Geom(lev).ProbLo(0);
    Real ylen = Geom(lev).ProbHi(1) - Geom(lev).ProbLo(1);
    Real zlen = Geom(lev).ProbHi(2) - Geom(lev).ProbLo(2);


    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {

	//number of particles
	const int np  =  NumberOfParticles(pti);
	auto& structs = pti.GetArrayOfStructs();
	auto& attribs = pti.GetStructOfArrays();
	
	// Temporary arrays
	Real         *pradius, *pdensity, *pvol, *pmass, *pomoi;
	Real         *ppos, *pvel, *pomega, *pacc, *palpha, *pdrag;
	int          *pstate, *pphase;
	int          nsubsteps;

	Array< Particle<realData::count, intData::count> > particles;
	
	particles.resize(np);
	
	
	// Fill in particle data
	for ( int i = 0; i < np; ++i ) {

	    particles[i].pos(0)   = structs[i].pos(0);  
	    particles[i].pos(1)   = structs[i].pos(1);
	    particles[i].pos(2)   = structs[i].pos(2);
	    particles[i].id()     = structs[i].id();
	    particles[i].cpu()    = structs[i].cpu();
	    particles[i].idata(0) = attribs.GetIntData(intData::phase)[i];
	    particles[i].idata(1) = attribs.GetIntData(intData::state)[i];
	    particles[i].rdata(0) = attribs.GetRealData(realData::radius)[i];
	    particles[i].rdata(1) = attribs.GetRealData(realData::volume)[i];
	    particles[i].rdata(2) = attribs.GetRealData(realData::mass)[i];
	    particles[i].rdata(3) = attribs.GetRealData(realData::density)[i];
	    particles[i].rdata(4) = attribs.GetRealData(realData::oneOverI)[i];
	    particles[i].rdata(5) = attribs.GetRealData(realData::velx)[i];
	    particles[i].rdata(6) = attribs.GetRealData(realData::vely)[i];
	    particles[i].rdata(7) = attribs.GetRealData(realData::velz)[i];
	    particles[i].rdata(8)  = attribs.GetRealData(realData::omegax)[i];
	    particles[i].rdata(9)  = attribs.GetRealData(realData::omegay)[i];
	    particles[i].rdata(10) = attribs.GetRealData(realData::omegaz)[i];
	    particles[i].rdata(11)  = attribs.GetRealData(realData::accx)[i];
	    particles[i].rdata(12)  = attribs.GetRealData(realData::accy)[i];
	    particles[i].rdata(13) = attribs.GetRealData(realData::accz)[i];
	    particles[i].rdata(14)  = attribs.GetRealData(realData::alphax)[i];
	    particles[i].rdata(15)  = attribs.GetRealData(realData::alphay)[i];
	    particles[i].rdata(16) = attribs.GetRealData(realData::alphaz)[i];
	    particles[i].rdata(17)  = attribs.GetRealData(realData::dragx)[i];
	    particles[i].rdata(18)  = attribs.GetRealData(realData::dragy)[i];
	    particles[i].rdata(19) = attribs.GetRealData(realData::dragz)[i];
		
	}

	mfix_des_init_time_loop_aos( &np, particles.dataPtr(), &time, &dt, &dx, &dy, &dz,
				     &xlen, &ylen, &zlen, &nstep, &nsubsteps);

	int quit;

	for ( int n = 0; n < nsubsteps; ++n ) {
	    mfix_des_time_loop_ops_aos( &np, particles.dataPtr(), &time, &dt, &dx, &dy, &dz,
					&xlen, &ylen, &zlen, &nstep, &quit);
		
	    if ( quit ) break;
	}


	mfix_des_finalize_time_loop_aos( &np, &dt, particles.dataPtr() );


	// Copy back particle data
	for ( int i = 0; i < np; ++i ) {

	    structs[i].pos(0) = particles[i].pos(0);  
	    structs[i].pos(1) = particles[i].pos(1);
	    structs[i].pos(2) = particles[i].pos(2);
	    structs[i].id()   = particles[i].id();
	    structs[i].cpu()  = particles[i].cpu();
	    attribs.GetIntData(intData::phase)[i] = particles[i].idata(0);
	    attribs.GetIntData(intData::state)[i] = particles[i].idata(1);

	    for ( int j = 0; j < realData::count; ++j )  {
		attribs.GetRealData(j)[i] = particles[i].rdata(j);
	    }
	}



	
    }

    Redistribute();

}







void MFIXParticleContainer::output(int lev, int estatus, int finish, int nstep, Real dt, Real time)
{

    Real xlen = Geom(lev).ProbHi(0) - Geom(lev).ProbLo(0);
    Real ylen = Geom(lev).ProbHi(1) - Geom(lev).ProbLo(1);
    Real zlen = Geom(lev).ProbHi(2) - Geom(lev).ProbLo(2);


    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {

	//number of particles
	const int np = NumberOfParticles(pti);

	// Temporary arrays
	Real         *ppos, *pvel, *pomega, *pradius;
	int          *pstate, *pphase;


	GetIntData( pti, intData::state, &pstate );

	GetRealData( pti, realData::radius, &pradius );

	GetPosition( pti, &ppos );

	GetVectorData( pti, realData::velx, &pvel );
	GetVectorData( pti, realData::omegax, &pomega );

	mfix_output_manager(&np, &time, &dt, &xlen, &ylen, &zlen, &nstep,
			    pstate, pradius, ppos, pvel,  pomega, &finish);


	RestoreIntData( &pstate );

	RestoreRealData( &pradius );

	RestorePosition( pti, &ppos );

	RestoreVectorData( pti, realData::velx, &pvel );
	RestoreVectorData( pti, realData::omegax, &pomega );

    }

}
