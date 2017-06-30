#ifndef _MFIXParticleContainer_H_
#define _MFIXParticleContainer_H_


#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include "AMReX_Particles.H"


using namespace amrex;
using namespace std;


namespace {

struct NeighborCommTag {

    NeighborCommTag(int pid, int gid, int tid)
        : proc_id(pid), grid_id(gid), tile_id(tid)
    {}

    int proc_id;
    int grid_id;
    int tile_id;
};

bool operator<(const NeighborCommTag& l, const NeighborCommTag& r) {
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
    using NeighborCommMap = std::map<NeighborCommTag, Array<char> >;

    // Get number of particles
    const int NumberOfParticles(MFIXParIter& pti) { return pti.GetArrayOfStructs().size(); }

    const int NumberOfParticles(int& lev, MFIter& mfi)
    {
      const int gridIndex = mfi.index();
      const int tileIndex = mfi.LocalTileIndex();
      auto&     particles = GetParticles(lev)[std::make_pair(gridIndex,tileIndex)];
      auto&     structs   = particles.GetArrayOfStructs();
      return structs.size();
    }

    MFIXParticleContainer (AmrCore* amr_core);

    virtual ~MFIXParticleContainer() {};

    void InitLevelMask ( int lev, const Geometry &geom, const DistributionMapping &dmap, const BoxArray &ba );

    void InitData ();

    void InitParticlesAscii (const std::string& file);

    void printParticles ();

    void AllocData ();

    void EvolveParticles ( int lev, int nstep, Real dt, Real time );

    void GetParticleAvgProp ( int lev, Real (&avg_dp)[10], Real (&avg_ro)[10]);

    ///
    /// Write out all particles at a level in the format fcompare_par can read
    ///
    void writeAllForComparison(int lev);

    void output ( int lev, int estatus, int finish, int nstep, Real dt, Real time );

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel(int lev = 0);

private:

    static void ReadStaticParameters ();

    ///
    /// Apply periodic shift to particle position, so naive distance calculation
    /// between neighbors and regular particles will be correct.
    ///
    void applyPeriodicShift(int lev, ParticleType& p, const amrex::IntVect& neighbor_cell);

    void fillNeighbors( int lev );

    ///
    /// Pack a particle's data into the proper neighbor buffer, or put it into
    /// the structure to be sent to the other processes
    ///
    void packNeighborParticle( int lev,
          const IntVect& neighbor_cell,
          const BaseFab<int>& mask,
          const ParticleType& p,
          NeighborCommMap& neighbors_to_comm);
    ///
    /// Perform the MPI communication neccesary to fill neighbor buffers
    ///
    void fillNeighborsMPI(NeighborCommMap& neighbors_to_comm);


    void clearNeighbors( int lev );

    //
    // This members are for the handling of neighbor particles
    //
    const int ng = 2;  // number of neighbor layers
    const size_t pdata_size = sizeof( ParticleType );
    // const size_t pdata_size = (BL_SPACEDIM + realData::count)*sizeof(RealType) +
    //  (2+intData::count)*sizeof(int);
    FabArray< BaseFab<int> > mask;
    std::map<PairIndex, Array<char> > neighbors;

};
#endif
