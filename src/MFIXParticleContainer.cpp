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
    : ParticleContainer<0,0,PIdx::nattribs>
      (amr_core->GetParGDB())
{
    ReadStaticParameters();
    
    this->SetVerbose(0);
}

void
MFIXParticleContainer::AllocData () {
    reserveData();
    resizeData();
}



void
MFIXParticleContainer::InitParticlesAscii(const std::string& file) {

    // only read the file on the IO proc
    if (ParallelDescriptor::MyProc() ==  ParallelDescriptor::IOProcessorNumber()) {
	std::ifstream ifs;
	ifs.open(file.c_str(), std::ios::in);

	if (!ifs.good())
	    amrex::FileOpenFailed(file);

	// int numberOfParticles = 0;
	ifs >> numberOfParticles >> std::ws;

	// Issue an error if nparticles = 0 is specified
	if ( numberOfParticles == 0 ){
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

	std::array<Real, PIdx::nattribs> attribs;

	pstate = 1;

	
	for (int i = 0; i < numberOfParticles; i++) {

	    // Init to 0 all attributes
	    attribs.fill(0.0);

	    // Read from input file
	    ifs >> pphase;
	    ifs >> p.pos(0);
	    ifs >> p.pos(1);
	    ifs >> p.pos(2);
	    ifs >> pradius;
	    ifs >> pdensity; 
	    ifs >> attribs[PIdx::velx];
	    ifs >> attribs[PIdx::vely];
	    ifs >> attribs[PIdx::velz];

	    // Set id and cpu for this particle
	    p.id()  = ParticleType::NextID();
	    p.cpu() = ParallelDescriptor::MyProc();

	    // Compute other particle properties
	    mfix_set_particle_properties( &pstate, &pradius, &pdensity, &pvolume, &pmass, &pomoi);

	    // Set other particle properties
	    attribs[PIdx::phase]    = pphase;
	    attribs[PIdx::state]    = pstate;
	    attribs[PIdx::volume]   = pvolume;
	    attribs[PIdx::density]  = pdensity;
	    attribs[PIdx::mass]     = pmass;
	    attribs[PIdx::oneOverI] = pomoi;
	    attribs[PIdx::radius]   = pradius;

	    // Add everything to the data structure
	    auto& particle_tile = GetParticles(lev)[std::make_pair(grid,tile)];
	    particle_tile.push_back(p);
	    particle_tile.push_back(attribs);
	}
    }
    Redistribute();

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
MFIXParticleContainer:: GetParticlesAttributes (
    Array<int>& particle_state, Array<int>& particle_phase,
    Array<Real>& des_radius,  Array<Real>& ro_sol,
    Array<Real>& pvol,   Array<Real>& pmass,
    Array<Real>& omoi,   Array<Real>& des_vel_new,
    Array<Real>& omega_new,   Array<Real>& des_acc_old,
    Array<Real>& rot_acc_old,   Array<Real>& drag_fc) {


    int lev = 0;
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {

	auto& attribs = pti.GetAttribs();

	des_radius = attribs[PIdx::radius];
	ro_sol     = attribs[PIdx::density];
	pvol       = attribs[PIdx::volume];
	pmass      = attribs[PIdx::mass];
	omoi       = attribs[PIdx::oneOverI];	
	
	auto& state    = attribs[PIdx::state];
	auto& phase    = attribs[PIdx::phase];
	auto& velx     = attribs[PIdx::velx];
	auto& vely     = attribs[PIdx::vely];
	auto& velz     = attribs[PIdx::velz];	
	auto& omegax   = attribs[PIdx::omegax];
	auto& omegay   = attribs[PIdx::omegay];
	auto& omegaz   = attribs[PIdx::omegaz];
	auto& accx     = attribs[PIdx::accx];
	auto& accy     = attribs[PIdx::accy];
	auto& accz     = attribs[PIdx::accz];	
	auto& alphax   = attribs[PIdx::alphax];
	auto& alphay   = attribs[PIdx::alphay];
	auto& alphaz   = attribs[PIdx::alphaz];	
	auto& dragx    = attribs[PIdx::dragx];
	auto& dragy    = attribs[PIdx::dragy];
	auto& dragz    = attribs[PIdx::dragz];
	
	Pack3DArrays(des_vel_new, velx, vely, velz);
	Pack3DArrays(omega_new, omegax, omegay, omegaz);
	Pack3DArrays(des_acc_old, accx, accy, accz);
	Pack3DArrays(rot_acc_old, alphax, alphay, alphaz);
	Pack3DArrays(drag_fc, dragx, dragy, dragz);
	
	// State and phase need type conversion
	for ( int i = 0; i < numberOfParticles; i++ )
	{
	    particle_state[i]  = state[i];
	    particle_phase[i]  = phase[i];

	}

    }
}


void
MFIXParticleContainer:: printParticles() {
    const int lev = 0;
    const auto& plevel = GetParticles(lev);
    for (const auto& kv : plevel) {
      const auto& particle_structs = kv.second.GetArrayOfStructs();
      const auto& particle_attribs = kv.second.GetStructOfArrays();
      for (unsigned i = 0; i < particle_structs.numParticles(); ++i) {
	  std::cout << "Particle ID  = " << i << " " << std::endl;
          for (int j = 0; j < PIdx::nattribs; j++) {
	      std::cout << particle_attribs[j][i] << " " << std::endl;
          }
          std::cout << std::endl;
      }
    }
}


void
MFIXParticleContainer::ReadStaticParameters ()
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


void
MFIXParticleContainer::EvolveParticles(Array< unique_ptr<MultiFab> >& ep_g, 
				       const Array< unique_ptr<MultiFab> >& u_g,
				       const Array< unique_ptr<MultiFab> >& v_g,
				       const Array< unique_ptr<MultiFab> >& w_g,
				       const Array< unique_ptr<MultiFab> >& p_g,
				       const Array< unique_ptr<MultiFab> >& ro_g,
				       const Array< unique_ptr<MultiFab> >& mu_g,
				       int lev, int nstep, Real dt, Real time )
{
    
    int pair_count = 0;

    Box domain(Geom(lev).Domain());

    Real dx = Geom(lev).CellSize(0);
    Real dy = Geom(lev).CellSize(1);
    Real dz = Geom(lev).CellSize(2);

    Real xlen = Geom(lev).ProbHi(0) - Geom(lev).ProbLo(0);
    Real ylen = Geom(lev).ProbHi(1) - Geom(lev).ProbLo(1);
    Real zlen = Geom(lev).ProbHi(2) - Geom(lev).ProbLo(2);

    const int np = NumberOfParticlesAtLevel( lev, 1, 1 );


    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) 
    {

	const Box& sbx = (*ep_g[lev])[pti].box();
	const Box& bx  = pti.validbox();
	
	Box ubx((*u_g[lev])[pti].box());
	Box vbx((*v_g[lev])[pti].box());
	Box wbx((*w_g[lev])[pti].box());

	
	// particle struct data
	auto& structs = pti.GetArrayOfStructs();
	
	// particle array data
	auto& attribs = pti.GetStructOfArrays();
	
	//number of particles
	const int np = structs.size();

	// Temporary arrays
	Array<int>   state, phase;
	Array<Real>   vel, pos, omega, acc, alpha, drag;

	state.resize(np);
	phase.resize(np);
	vel.resize(3*np);
	pos.resize(3*np);
	omega.resize(3*np);
	acc.resize(3*np);
	alpha.resize(3*np);
	drag.resize(3*np);


	for( unsigned i=0; i < np ; ++i )
	{
	    pos[i]      = structs[i].pos(0);
	    pos[i+np]   = structs[i].pos(1);
	    pos[i+2*np] = structs[i].pos(2);
	    state[i]    = int(attribs[PIdx::state][i]);
	    phase[i]    = int(attribs[PIdx::phase][i]);
	}


	Pack3DArrays(vel,   attribs[PIdx::velx], attribs[PIdx::vely], attribs[PIdx::velz] );
	Pack3DArrays(acc,   attribs[PIdx::accx], attribs[PIdx::accy], attribs[PIdx::accz] );
	Pack3DArrays(omega, attribs[PIdx::omegax], attribs[PIdx::omegay], attribs[PIdx::omegaz] );
	Pack3DArrays(alpha, attribs[PIdx::alphax], attribs[PIdx::alphay], attribs[PIdx::alphaz] );
	Pack3DArrays(drag,  attribs[PIdx::dragx], attribs[PIdx::dragy], attribs[PIdx::dragz] );


	mfix_des_time_march( &np, 
			     sbx.loVect(), sbx.hiVect(),
			     ubx.loVect(), ubx.hiVect(),
			     vbx.loVect(), vbx.hiVect(),
			     wbx.loVect(), wbx.hiVect(),
			     bx.loVect(),  bx.hiVect(),
			     domain.loVect(), domain.hiVect(),
			     (*ep_g[lev])[pti].dataPtr(), (*p_g[lev])[pti].dataPtr(),
			     (*u_g[lev])[pti].dataPtr(),  (*v_g[lev])[pti].dataPtr(), (*w_g[lev])[pti].dataPtr(),
			     (*ro_g[lev])[pti].dataPtr(), (*mu_g[lev])[pti].dataPtr(),
			     state.dataPtr(),    phase.dataPtr(),
			     attribs[PIdx::radius].dataPtr(), attribs[PIdx::volume].dataPtr(), attribs[PIdx::mass].dataPtr(),
			     attribs[PIdx::oneOverI].dataPtr(), pos.dataPtr(),    vel.dataPtr(),   omega.dataPtr(),
			     acc.dataPtr(),    alpha.dataPtr(), drag.dataPtr(), &time, &dt, &dx, &dy, &dz, 
			     &xlen, &ylen, &zlen, &nstep);

	for( unsigned i=0; i < np ; ++i )
	{
	    structs[i].pos(0) = pos[i];
	    structs[i].pos(1) = pos[i+np] ;
	    structs[i].pos(2) = pos[i+2*np];
	    attribs[PIdx::state][i] = state[i];
	    attribs[PIdx::phase][i] = phase[i];
	}

	Unpack3DArrays(attribs[PIdx::velx], attribs[PIdx::vely], attribs[PIdx::velz], vel );
	Unpack3DArrays(attribs[PIdx::accx], attribs[PIdx::accy], attribs[PIdx::accz], acc );
	Unpack3DArrays(attribs[PIdx::omegax], attribs[PIdx::omegay], attribs[PIdx::omegaz], omega );
	Unpack3DArrays(attribs[PIdx::alphax], attribs[PIdx::alphay], attribs[PIdx::alphaz], alpha );
	Unpack3DArrays(attribs[PIdx::dragx], attribs[PIdx::dragy], attribs[PIdx::dragz], drag );

    }

}

// void
// MFIXParticleContainer::Evolve (int lev, Real dt)
// {
//     BL_PROFILE("MyPC::Evolve()");

//     const Geometry& gm  = m_gdb->Geom(lev);

// #if (BL_SPACEDIM == 3)
//     const Real* dx = gm.CellSize();
// #elif (BL_SPACEDIM == 2)
//     Real dx[3] = { gm.CellSize(0), std::numeric_limits<Real>::quiet_NaN(), gm.CellSize(1) };
// #endif
// }
