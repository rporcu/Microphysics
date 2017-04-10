#include <AMReX.H>
#include "AMReX_Particles.H"
#include <iostream>
#include <MFIXParticleContainer.H>

#include<math.h>

#include "mfix_F.H"


using namespace amrex;

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
  
    const Real four  = 4.0;
    const Real three = 3.0;
    const Real half5 = 2.5;
      
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
    	std::array<Real, PIdx::nattribs> attribs;

    	for (int i = 0; i < numberOfParticles; i++) {
     
    	    // Here we read the struct data into the particle.
    	    ifs >> p.pos(0);
    	    ifs >> p.pos(1);
    	    ifs >> p.pos(2);
    	    p.id() = ParticleType::NextID();
    	    p.cpu() = ParallelDescriptor::MyProc();

	    // Init to 0 every attributes
	    attribs.fill(0.0);
                
    	    // Initialize attributes given in input file
    	    ifs >> attribs[PIdx::radius];
    	    ifs >> attribs[PIdx::density];
    	    ifs >> attribs[PIdx::velx];
    	    ifs >> attribs[PIdx::vely];
    	    ifs >> attribs[PIdx::velz];

	    // Initialize remaning attributes
	    // State = 1 corresponds to parameter "normal_particle" defined in module discretelement
	    attribs[PIdx::phase]    = 0;
	    attribs[PIdx::state]    = 1; 
	    attribs[PIdx::volume]   = (four/three)*M_PI*pow(attribs[PIdx::radius],3);	    
	    attribs[PIdx::mass]     = attribs[PIdx::volume] * attribs[PIdx::density];
	    attribs[PIdx::oneOverI] = half5/(attribs[PIdx::mass]*pow(attribs[PIdx::radius],2));

	    // Call Fortran external to set phase of particle
	    mfix_set_phase_index( &attribs[PIdx::phase], &attribs[PIdx::radius], 
	    			  &attribs[PIdx::density], &attribs[PIdx::state], &i );
 
    	    // add them to the data structure
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

    // Fill des_pos_new
    for ( int i = 0; i < numberOfParticles; ++i )
    {
	des_pos_new[3*i]   = xp[i];
	des_pos_new[3*i+1] = yp[i];
	des_pos_new[3*i+2] = zp[i];
    }
    
}

void
MFIXParticleContainer:: GetParticlesAttributes (
    Array<int>& particle_state, Array<int>& particle_phase,
    Array<Real>& des_radius,  Array<Real>& ro_sol,
    Array<Real>& pvol,   Array<Real>& pmass,
    Array<Real>& omoi,   Array<Real>& des_vel_new,
    Array<Real>& omega_new,   Array<Real>& des_acc_old,
    Array<Real>& rot_acc_old,   Array<Real>& drag_fc,
    Array<Real>& fc,   Array<Real>& tow) {


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
	auto& cforcex  = attribs[PIdx::cforcex];
	auto& cforcey  = attribs[PIdx::cforcey];
	auto& cforcez  = attribs[PIdx::cforcez];
	auto& acforcex = attribs[PIdx::acforcex];
	auto& acforcey = attribs[PIdx::acforcey];
	auto& acforcez = attribs[PIdx::acforcez];
	
	Pack3DArrays(des_vel_new, velx, vely, velz);
	Pack3DArrays(omega_new, omegax, omegay, omegaz);
	Pack3DArrays(des_acc_old, accx, accy, accz);
	Pack3DArrays(rot_acc_old, alphax, alphay, alphaz);
	Pack3DArrays(drag_fc, dragx, dragy, dragz);
	Pack3DArrays(fc, cforcex, cforcey, cforcez);
	Pack3DArrays(tow, acforcex, acforcey, acforcez);
	
	// State and phase need type conversion
	for ( int i = 0; i < numberOfParticles; i++ )
	{
	    particle_state[i]  = int(state[i]);
	    particle_phase[i]  = int(phase[i]);

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
MFIXParticleContainer::Pack3DArrays( Array<Real>& vec, Array<Real>& comp1,
				     Array<Real>& comp2, Array<Real>& comp3 )
{
    for ( int i = 0; i < comp1.size(); i++ )
    {
	vec[3*i]   = comp1[i];
	vec[3*i+1] = comp2[i];
	vec[3*i+2] = comp3[i];
    }
    
}

    
void
MFIXParticleContainer::Evolve (int lev, Real dt)
{
    BL_PROFILE("MyPC::Evolve()");

    const Geometry& gm  = m_gdb->Geom(lev);

#if (BL_SPACEDIM == 3)
    const Real* dx = gm.CellSize();
#elif (BL_SPACEDIM == 2)
    Real dx[3] = { gm.CellSize(0), std::numeric_limits<Real>::quiet_NaN(), gm.CellSize(1) };
#endif
}




