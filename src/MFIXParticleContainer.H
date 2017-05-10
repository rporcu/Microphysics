#ifndef _MFIXParticleContainer_H_
#define _MFIXParticleContainer_H_


#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include "AMReX_Particles.H"


using namespace amrex;
using namespace std;


namespace {

struct GhostCommTag {

    GhostCommTag(int pid, int gid, int tid)
        : proc_id(pid), grid_id(gid), tile_id(tid)
    {}

    int proc_id;
    int grid_id;
    int tile_id;
};

bool operator<(const GhostCommTag& l, const GhostCommTag& r) {
    return (l.proc_id < r.proc_id || 
           (l.proc_id == r.proc_id && l.grid_id < r.grid_id) ||
           (l.proc_id == r.proc_id && l.grid_id == r.grid_id && l.tile_id < r.tile_id ));
}

}


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
    : public amrex::ParIter<realData::count,intData::count,0,0>
{
public:
    using amrex::ParIter<realData::count,intData::count,0,0>::ParIter;
};

class MFIXParticleContainer
    : public amrex::ParticleContainer<realData::count,intData::count,0,0>
{
public:

    using PairIndex = std::pair<int, int>;
    using GhostCommMap = std::map<GhostCommTag, Array<char> >;
    
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

    void InitLevelMask ( int lev, const Geometry &geom, const DistributionMapping &dmap, const BoxArray &ba );
    
    void InitData ();

    void InitParticlesAscii (const std::string& file);

    void printParticles ();

    void AllocData ();

    void* GetParticlesData( const int& lev, const MFIter& mfi );

    void* GetParticlesData( MFIXParIter& pti);    

    void EvolveParticles ( int lev, int nstep, Real dt, Real time );
    
    void output ( int lev, int estatus, int finish, int nstep, Real dt, Real time );

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel(int lev = 0);

private:

    static void ReadStaticParameters ();

    void fillGhosts( int lev );

    ///
    /// Pack a particle's data into the proper neighbor buffer, or put it into
    /// the structure to be sent to the other processes
    ///
    void packGhostParticle( int lev, 
			    const IntVect& neighbor_cell,
			    const BaseFab<int>& mask,
			    const ParticleType& p,
			    GhostCommMap& ghosts_to_comm);
    ///
    /// Perform the MPI communication neccesary to fill ghost buffers
    ///
    void fillGhostsMPI(GhostCommMap& ghosts_to_comm);


    void clearGhosts( int lev );
    
    // runtime parameters
    static int do_tiling;
    static IntVect tile_size;

    //
    // This members are for the handling of ghost particles
    // 
    const int ng = 2;  // number of ghost layers
    const size_t pdata_size = sizeof( ParticleType );
    // const size_t pdata_size = (BL_SPACEDIM + realData::count)*sizeof(RealType) +
    // 	(2+intData::count)*sizeof(int);
    FabArray< BaseFab<int> > mask;
    std::map<PairIndex, Array<char> > ghosts; 

};
#endif
