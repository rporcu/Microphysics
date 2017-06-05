#include <AMReX.H>
#include "AMReX_Particles.H"
#include "AMReX_RealVect.H"
#include <iostream>
#include <MFIXParticleContainer.H>

#include<math.h>

#include "mfix_F.H"


using namespace amrex;
using namespace std;


int     MFIXParticleContainer::do_tiling = 0;
IntVect MFIXParticleContainer::tile_size   { D_DECL(1024000,8,8) };

MFIXParticleContainer::MFIXParticleContainer (AmrCore* amr_core)
    : ParticleContainer<realData::count,intData::count,0,0>
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

void MFIXParticleContainer::InitLevelMask ( int lev,
              const Geometry &geom,
              const DistributionMapping &dmap,
              const BoxArray &ba )
{
    BL_ASSERT( lev == 0 );

    mask.define(ba, dmap, 2, ng);
    mask.setVal(-1, ng);
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        mask.setVal(grid_id, box, 0, 1);
        mask.setVal(tile_id, box, 1, 1);
    }
    mask.FillBoundary(geom.periodicity());
}

void* MFIXParticleContainer::GetParticlesData( const int& lev, const MFIter& mfi ) {

    const int gridIndex = mfi.index();
    const int tileIndex = mfi.LocalTileIndex();
    auto&     particles = GetParticles(lev)[std::make_pair(gridIndex,tileIndex)];

    void* ptr = NULL;

    if ( particles.GetArrayOfStructs().size() > 0 )
  ptr = particles.GetArrayOfStructs().data();

    return  ptr; //particles.GetArrayOfStructs().data();
}


void* MFIXParticleContainer::GetParticlesData( MFIXParIter& pti ) {

    void* ptr = NULL;

    if ( pti.GetArrayOfStructs().size() > 0 )
  ptr = pti.GetArrayOfStructs().data();

    return  ptr;
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
  Real       pradius, pdensity, pvolume, pomoi, pmass, pomega;

  pstate = 1;

  for (int i = 0; i < np; i++) {

      // Read from input file
      ifs >> pphase;
      ifs >> p.pos(0);
      ifs >> p.pos(1);
      ifs >> p.pos(2);
      ifs >> pradius;
      ifs >> pdensity;
      ifs >> p.rdata(realData::velx);
      ifs >> p.rdata(realData::vely);
      ifs >> p.rdata(realData::velz);

      // Set id and cpu for this particle
      p.id()  = ParticleType::NextID();
      p.cpu() = ParallelDescriptor::MyProc();

      // Compute other particle properties
      mfix_set_particle_properties( &pstate, &pradius, &pdensity, &pvolume, &pmass, &pomoi, &pomega);

      // Set other particle properties
      p.idata(intData::phase)     = pphase;
      p.idata(intData::state)     = pstate;
      p.rdata(realData::volume)   = pvolume;
      p.rdata(realData::density)  = pdensity;
      p.rdata(realData::mass)     = pmass;
      p.rdata(realData::oneOverI) = pomoi;
      p.rdata(realData::radius)   = pradius;
      p.rdata(realData::omegax)   = pomega;
      p.rdata(realData::omegay)   = pomega;
      p.rdata(realData::omegaz)   = pomega;

      // Add everything to the data structure
      auto& particle_tile = GetParticles(lev)[std::make_pair(grid,tile)];
      particle_tile.push_back(p);
  }
    }
    Redistribute();

}


void MFIXParticleContainer:: printParticles() {
    const int lev = 0;
    const auto& plevel = GetParticles(lev);
    for (const auto& kv : plevel) {
  const auto& particles = kv.second.GetArrayOfStructs();

  for (unsigned i = 0; i < particles.numParticles(); ++i) {

      std::cout << "Particle ID  = " << i << " " << std::endl;
      std::cout << "X            = " << particles[i].pos(0) << " " << std::endl;
      std::cout << "Y            = " << particles[i].pos(1) << " " << std::endl;
      std::cout << "Z            = " << particles[i].pos(2) << " " << std::endl;
      std::cout << "state        = " << particles[i].idata(intData::state) << " " << std::endl;
      std::cout << "phase        = " << particles[i].idata(intData::phase) << " " << std::endl;
      std::cout << "Real properties = " << std::endl;

      for (int j = 0; j < realData::count; j++) {
    std::cout << "property " << j << "  = " << particles[i].rdata(j) << " " << std::endl;
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

void MFIXParticleContainer::EvolveParticles( int lev, int nstep, Real dt, Real time ) {


    Box domain(Geom(lev).Domain());

    Real dx = Geom(lev).CellSize(0);
    Real dy = Geom(lev).CellSize(1);
    Real dz = Geom(lev).CellSize(2);

    Real xlen = Geom(lev).ProbHi(0) - Geom(lev).ProbLo(0);
    Real ylen = Geom(lev).ProbHi(1) - Geom(lev).ProbLo(1);
    Real zlen = Geom(lev).ProbHi(2) - Geom(lev).ProbLo(2);

    int   nsubsteps;
    Real  subdt;

    mfix_des_init_time_loop( &time, &dt, &nsubsteps, &subdt );

    for ( int n = 0; n < nsubsteps; ++n ) {

        fillNeighbors(lev);

  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

      // Real particles
      const int np     = NumberOfParticles(pti);
      void* particles  = GetParticlesData(pti);

      // Neighbor particles
      int nstride = pti.GetArrayOfStructs().dataShape().first;
      PairIndex index(pti.index(), pti.LocalTileIndex());
      int ng = neighbors[index].size() / pdata_size;

      mfix_des_time_loop_ops( &np, particles, &ng, neighbors[index].dataPtr(),
            &subdt, &dx, &dy, &dz,
            &xlen, &ylen, &zlen, &nstep );

      if ( mfix_des_continuum_coupled () == 0 ) {
    Real stime;
    stime = time + (n+1)*subdt;
    mfix_output_manager( &np, &stime, &subdt,  &xlen, &ylen, &zlen,
             &n, particles, 0 );

      }

      mfix_call_usr2_des( &np, particles );
  }

        clearNeighbors(lev);

        Redistribute();
    }

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

  const int np     = NumberOfParticles(pti);
  void* particles  = GetParticlesData(pti);

  mfix_call_usr3_des( &np, particles );

    }

    if ( mfix_des_continuum_coupled () != 0 ) {
      nstep = nsubsteps;
      time  = time + nsubsteps * subdt ;
    }

    // Redistribute();
}


void MFIXParticleContainer::output(int lev, int estatus, int finish, int nstep, Real dt, Real time)
{

    Real xlen = Geom(lev).ProbHi(0) - Geom(lev).ProbLo(0);
    Real ylen = Geom(lev).ProbHi(1) - Geom(lev).ProbLo(1);
    Real zlen = Geom(lev).ProbHi(2) - Geom(lev).ProbLo(2);


    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {

  //number of particles
  const int     np = NumberOfParticles(pti);
  void*  particles = GetParticlesData(pti);

  mfix_output_manager( &np, &time, &dt, &xlen, &ylen, &zlen, &nstep,
           particles, &finish);
    }

}

void MFIXParticleContainer::fillNeighbors( int lev ) {
    NeighborCommMap neighbors_to_comm;
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {
        const Box& tile_box = pti.tilebox();
        const IntVect& lo = tile_box.smallEnd();
        const IntVect& hi = tile_box.bigEnd();

        Box shrink_box = pti.tilebox();
        shrink_box.grow(-ng);

        auto& particles = pti.GetArrayOfStructs();
        for (unsigned i = 0; i < pti.numParticles(); ++i) {
            const ParticleType& p = particles[i];
            const IntVect& iv = Index(p, lev);

            // if the particle is more than one cell away from
            // the tile boundary, it's not anybody's neighbor
            if (shrink_box.contains(iv)) continue;

            // shift stores whether we are near the tile boundary in each direction.
            // -ng means lo, ng means hi, 0 means not near the boundary
            IntVect shift = IntVect::TheZeroVector();
            for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
                if (iv[idim] == lo[idim])
                    shift[idim] = -ng;
                else if (iv[idim] == hi[idim])
                    shift[idim] = ng;
            }

            // Based on the value of shift, we add the particle to a map to be sent
            // to the neighbors. A particle can be sent to up to 3 neighbors in 2D
            // and up to 7 in 3D, depending on whether is near the tile box corners,
            // edges, or just the faces. First, add the particle for the "face" neighbors
            for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
                if (shift[idim] == 0) continue;
                IntVect neighbor_cell = iv;
                neighbor_cell.shift(idim, shift[idim]);
                BL_ASSERT(mask[pti].box().contains(neighbor_cell));
                packNeighborParticle(lev,neighbor_cell, mask[pti], p, neighbors_to_comm);
            }

            // Now add the particle to the "edge" neighbors
            for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
                for (int jdim = 0; jdim < idim; ++jdim) {
                    if (shift[idim] != 0 and shift[jdim] != 0) {
                        IntVect neighbor_cell = iv;
                        neighbor_cell.shift(idim, shift[idim]);
                        neighbor_cell.shift(jdim, shift[jdim]);
                        BL_ASSERT(mask[pti].box().contains(neighbor_cell));
                        packNeighborParticle(lev,neighbor_cell, mask[pti], p, neighbors_to_comm);
                    }
                }
            }

            // Finally, add the particle for the "vertex" neighbors (only relevant in 3D)
            if (shift[0] != 0 and shift[1] != 0 and shift[2] != 0) {
                IntVect neighbor_cell = iv;
                neighbor_cell.shift(shift);
                BL_ASSERT(mask[pti].box().contains(neighbor_cell));
                packNeighborParticle(lev,neighbor_cell, mask[pti], p, neighbors_to_comm);
            }

        }
    }
    fillNeighborsMPI(neighbors_to_comm);
}

void MFIXParticleContainer::applyPeriodicShift(int lev, ParticleType& p,
                                              const IntVect& neighbor_cell) {

    const Periodicity& periodicity = Geom(lev).periodicity();
    if (not periodicity.isAnyPeriodic()) return;

    const Box& domain = Geom(lev).Domain();
    const IntVect& lo = domain.smallEnd();
    const IntVect& hi = domain.bigEnd();
    const RealBox& prob_domain = Geom(lev).ProbDomain();

    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
        if (not periodicity.isPeriodic(dir)) continue;
        if (neighbor_cell[dir] < lo[dir]) {
            p.pos(dir) += prob_domain.length(dir);
        }
        else if (neighbor_cell[dir] > hi[dir]) {
            p.pos(dir) -= prob_domain.length(dir);
        }
    }
}

void MFIXParticleContainer::packNeighborParticle(int lev,
                const IntVect& neighbor_cell,
                const BaseFab<int>& mask,
                const ParticleType& p,
                NeighborCommMap& neighbors_to_comm) {
    const int neighbor_grid = mask(neighbor_cell, 0);
    if (neighbor_grid >= 0) {
        const int who = ParticleDistributionMap(lev)[neighbor_grid];
        const int MyProc = ParallelDescriptor::MyProc();
        const int neighbor_tile = mask(neighbor_cell, 1);
        PairIndex dst_index(neighbor_grid, neighbor_tile);
        ParticleType particle = p;
        applyPeriodicShift(lev, particle, neighbor_cell);
        if (who == MyProc) {
            size_t old_size = neighbors[dst_index].size();
            size_t new_size = neighbors[dst_index].size() + pdata_size;
            neighbors[dst_index].resize(new_size);
            std::memcpy(&neighbors[dst_index][old_size], &particle, pdata_size);
        } else {
            NeighborCommTag tag(who, neighbor_grid, neighbor_tile);
            Array<char>& buffer = neighbors_to_comm[tag];
            size_t old_size = buffer.size();
            size_t new_size = buffer.size() + pdata_size;
            buffer.resize(new_size);
            std::memcpy(&buffer[old_size], &particle, pdata_size);
        }
    }
}

void MFIXParticleContainer::fillNeighborsMPI(NeighborCommMap& neighbors_to_comm) {

#ifdef BL_USE_MPI
    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();

    // count the number of tiles to be sent to each proc
    std::map<int, int> tile_counts;
    for (const auto& kv: neighbors_to_comm) {
        tile_counts[kv.first.proc_id] += 1;
    }

    // flatten all the data for each proc into a single buffer
    // once this is done, each dst proc will have an Array<char>
    // the buffer will be packed like:
    // ntiles, gid1, tid1, size1, data1....  gid2, tid2, size2, data2... etc.
    std::map<int, Array<char> > send_data;
    for (const auto& kv: neighbors_to_comm) {
        Array<char>& buffer = send_data[kv.first.proc_id];
        buffer.resize(sizeof(int));
        std::memcpy(&buffer[0], &tile_counts[kv.first.proc_id], sizeof(int));
    }

    for (auto& kv : neighbors_to_comm) {
        int data_size = kv.second.size();
        Array<char>& buffer = send_data[kv.first.proc_id];
        size_t old_size = buffer.size();
        size_t new_size = buffer.size() + 2*sizeof(int) + sizeof(int) + data_size;
        buffer.resize(new_size);
        char* dst = &buffer[old_size];
        std::memcpy(dst, &(kv.first.grid_id), sizeof(int)); dst += sizeof(int);
        std::memcpy(dst, &(kv.first.tile_id), sizeof(int)); dst += sizeof(int);
        std::memcpy(dst, &data_size,          sizeof(int)); dst += sizeof(int);
        if (data_size == 0) continue;
        std::memcpy(dst, &kv.second[0], data_size);
        Array<char>().swap(kv.second);
    }

    // each proc figures out how many bytes it will send, and how
    // many it will receive
    Array<long> snds(NProcs, 0), rcvs(NProcs, 0);
    long num_snds = 0;
    for (const auto& kv : send_data) {
        num_snds      += kv.second.size();
        snds[kv.first] = kv.second.size();
    }
    ParallelDescriptor::ReduceLongMax(num_snds);
    if (num_snds == 0) return;

    // communicate that information
    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

    BL_MPI_REQUIRE( MPI_Alltoall(snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(rcvs[MyProc] == 0);

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    Array<int> RcvProc;
    Array<std::size_t> rOffset; // Offset (in bytes) in the receive buffer

    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes);
            TotRcvBytes += rcvs[i];
        }
    }

    const int nrcvs = RcvProc.size();
    Array<MPI_Status>  stats(nrcvs);
    Array<MPI_Request> rreqs(nrcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();

    // Allocate data for rcvs as one big chunk.
    Array<char> recvdata(TotRcvBytes);

    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt    = rcvs[Who];

        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);

        rreqs[i] = ParallelDescriptor::Arecv(&recvdata[offset], Cnt, Who, SeqNum).req();
    }

    // Send.
    for (const auto& kv : send_data) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();

        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());

        ParallelDescriptor::Send(kv.second.data(), Cnt, Who, SeqNum);
    }

    // unpack the received data and put them into the proper neighbor buffers
    if (nrcvs > 0) {
        BL_MPI_REQUIRE( MPI_Waitall(nrcvs, rreqs.data(), stats.data()) );
        for (int i = 0; i < nrcvs; ++i) {
            const int offset = rOffset[i];
            char* buffer = &recvdata[offset];
            int num_tiles, gid, tid, size;
            std::memcpy(&num_tiles, buffer, sizeof(int)); buffer += sizeof(int);
            for (int j = 0; j < num_tiles; ++j) {
                std::memcpy(&gid,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&tid,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&size, buffer, sizeof(int)); buffer += sizeof(int);

                if (size == 0) continue;

                PairIndex dst_index(gid, tid);
                size_t old_size = neighbors[dst_index].size();
                size_t new_size = neighbors[dst_index].size() + size;
                neighbors[dst_index].resize(new_size);
                std::memcpy(&neighbors[dst_index][old_size], buffer, size); buffer += size;
            }
        }
    }
#endif
}

void MFIXParticleContainer::clearNeighbors( int lev )
{
    neighbors.clear();
}

void MFIXParticleContainer::writeAllAtLevel(int lev)
{
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
	auto& particles = pti.GetArrayOfStructs();
	size_t Np = pti.numParticles();
	cout << "Particles: " << Np << " << at level " << lev << endl;
	for (unsigned i = 0; i < Np; ++i)
	{
	    const ParticleType& p = particles[i];
	    const IntVect& iv = Index(p, lev);

	    RealVect xyz(p.pos(0), p.pos(1), p.pos(2));

	    cout << "[" << i << "]: id " << p.id()
		 << " index " << iv
		 << " position " << xyz << endl;
	}
    }
}

void MFIXParticleContainer::writeAllForComparison(int lev)
{
  size_t Np_tot = 0;

  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
      Np_tot += pti.numParticles();
 
  cout << Np_tot << std::endl;
 
  Real dummy = 0.;
 
  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& particles = pti.GetArrayOfStructs();
      size_t Np = pti.numParticles();
      for (unsigned i = 0; i < Np; ++i)
        {
          const ParticleType& p = particles[i];
 
          cout << p.pos(0)   << " " << p.pos(1)   << " " << p.pos(2) <<  " " <<
                  p.rdata(1) << " " << p.rdata(2) << " " << p.rdata(2) <<  " " <<
                  ParallelDescriptor::MyProc() << " " << p.id() << std::endl;
        }
    }
}

