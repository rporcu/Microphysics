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

	std::array<Real, PIdx::nattribs> attribs;

	pstate = 1;

	
	for (int i = 0; i < np; i++) {

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

	const Box& sbx = (*ep_g[lev])[pti].box();
	const Box& bx  = pti.validbox();
	
	Box ubx((*u_g[lev])[pti].box());
	Box vbx((*v_g[lev])[pti].box());
	Box wbx((*w_g[lev])[pti].box());
	
	//number of particles
	const int np =  NumberOfParticles(pti);

	// Temporary arrays
	Real         *pradius, *pdensity, *pvol, *pmass, *pomoi;
	Real         *ppos, *pvel, *pomega, *pacc, *palpha, *pdrag; 
    	int          *pstate, *pphase;


	GetParticlesAttributes ( pti, &pstate, &pphase, &pradius,  NULL, &pvol, &pmass,
				 &pomoi, &ppos, &pvel, &pomega, &pacc, &palpha, &pdrag ); 


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
			     pstate, pphase, pradius, pvol, pmass, pomoi, ppos,  pvel, pomega,
			     pacc,  palpha, pdrag, &time, &dt, &dx, &dy, &dz, 
			     &xlen, &ylen, &zlen, &nstep);


	RestoreParticlesAttributes  ( pti, &pstate, &pphase, &pradius,  NULL, &pvol, &pmass,
				      &pomoi, &ppos, &pvel, &pomega, &pacc, &palpha, &pdrag ); 


    }

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


	GetParticlesAttributes ( pti, &pstate, NULL, &pradius,  NULL, NULL, NULL,
				 NULL, &ppos, &pvel, &pomega, NULL, NULL, NULL ); 

	mfix_output_manager(&np, &time, &dt, &xlen, &ylen, &zlen, &nstep,
			    pstate, pradius, ppos, pvel,  pomega, &finish);


	RestoreParticlesAttributes ( pti, &pstate, NULL, &pradius,  NULL, NULL, NULL,
				     NULL, &ppos, &pvel, &pomega, NULL, NULL, NULL);  

    }

}


