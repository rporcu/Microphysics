#include <AMReX.H>
#include "AMReX_Particles.H"
#include <iostream>
#include <MFIXParticleContainer.H>

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
        
    // only read the file on the IO proc
    if (ParallelDescriptor::MyProc() ==  ParallelDescriptor::IOProcessorNumber()) {
    	std::ifstream ifs;
    	ifs.open(file.c_str(), std::ios::in);
        
    	if (!ifs.good())
    	    amrex::FileOpenFailed(file);
        
    	int num_particles = 0;
    	ifs >> num_particles >> std::ws;

    	// we add all the particles to grid 0 and tile 0 and let 
    	// Redistribute() put them in the right places.
    	const int lev = 0;
    	const int grid = 0;
    	const int tile = 0;

    	ParticleType p;
    	std::array<Real, PIdx::nattribs> attribs;
    	for (int i = 0; i < num_particles; i++) {
     
    	    // Here we read the struct data into the particle.
    	    ifs >> p.pos(0);
    	    ifs >> p.pos(1);
    	    ifs >> p.pos(2);
    	    p.id() = ParticleType::NextID();
    	    p.cpu() = ParallelDescriptor::MyProc();
                
    	    // here were read the additional attributes to be stored in SoA style.
    	    ifs >> attribs[PIdx::radius];
    	    ifs >> attribs[PIdx::density];
    	    ifs >> attribs[PIdx::velx];
    	    ifs >> attribs[PIdx::vely];
    	    ifs >> attribs[PIdx::velz];
	                    
    	    // add them to the data structure
    	    auto& particle_tile = GetParticles(lev)[std::make_pair(grid,tile)];
    	    particle_tile.push_back(p);
    	    particle_tile.push_back(attribs);
    	}
    }
    Redistribute();
}


void
MFIXParticleContainer:: printParticles() {
    const int lev = 0;
    const auto& plevel = GetParticles(lev);
    for (const auto& kv : plevel) {
    	const auto& particle_structs = kv.second.GetArrayOfStructs();
    	const auto& particle_attribs = kv.second.GetStructOfArrays();
    	for (unsigned i = 0; i < particle_structs.numParticles(); ++i) {
    	    std::cout << particle_structs[i] << " ";
    	    for (int j = 0; j < PIdx::nattribs; j++) {
    		std::cout << particle_attribs[j][i] << " ";
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




