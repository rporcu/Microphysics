#include <AMReX.H>
#include "AMReX_Particles.H"
#include "AMReX_RealVect.H"
#include <iostream>
#include <MFIXParticleContainer.H>
#include <AMReX_LoadBalanceKD.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB_F.H>

#include <math.h>

#include "mfix_F.H"
#include "mfix_des_F.H"
#include "mfix_eb_F.H"

using namespace amrex;
using namespace std;

bool MFIXParticleContainer::use_neighbor_list  {true};
bool MFIXParticleContainer::sort_neighbor_list {false};

MFIXParticleContainer::MFIXParticleContainer (AmrCore* amr_core)
    : NeighborParticleContainer<realData::count,intData::count>
      (amr_core->GetParGDB(), 1)
{
    ReadStaticParameters();

    this->SetVerbose(0);

    // turn off certain components for ghost particle communication
    setRealCommComp(4, false);
    setRealCommComp(5, false);
    setRealCommComp(6, false);
    setRealCommComp(7, false);
    setRealCommComp(14, false);
    setRealCommComp(15, false);
    setRealCommComp(16, false);

    setIntCommComp(0, false);
    setIntCommComp(1, false);
    setIntCommComp(3, false);

    nlev = amr_core->maxLevel() + 1;
}

void MFIXParticleContainer::AllocData ()
{
    reserveData();
    resizeData();
}

void MFIXParticleContainer::InitParticlesAscii(const std::string& file) {

  // only read the file on the IO proc
  if (ParallelDescriptor::IOProcessor())  {
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

    auto& particle_tile = GetParticles(lev)[std::make_pair(grid,tile)];

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
      set_particle_properties( &pstate, &pradius, &pdensity,
                               &pvolume, &pmass, &pomoi, &pomega);

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

      // Initialize these for I/O purposes
      p.rdata(realData::dragx)    = 0.0;
      p.rdata(realData::dragy)    = 0.0;
      p.rdata(realData::dragz)    = 0.0;

      // Add everything to the data structure
      particle_tile.push_back(p);

      if (!ifs.good())
          amrex::Abort("Error initializing particles from Ascii file. \n");
    }
  }
  Redistribute();
}

void MFIXParticleContainer::InitParticlesAuto()
{
  int lev = 0;

  Box domain(Geom(lev).Domain());

  Real dx = Geom(lev).CellSize(0);
  Real dy = Geom(lev).CellSize(1);
  Real dz = Geom(lev).CellSize(2);

  int total_np = 0;

  // This uses the particle tile size. Note that the default is to tile so if we
  //      remove the true and don't explicitly add false it will still tile
  for (MFIter mfi = MakeMFIter(lev,true); mfi.isValid(); ++mfi) {

      // This is particles per grid so we reset to 0
      int pcount = 0;

      // Define the real particle data for one grid-at-a-time's worth of particles
      // We don't know pcount (number of particles per grid) before this call

      const Box& tilebx = mfi.tilebox();

      mfix_particle_generator(&pcount, tilebx.loVect(), tilebx.hiVect(), &dx, &dy, &dz);

      const int grid_id = mfi.index();
      const int tile_id = mfi.LocalTileIndex();

      // Now that we know pcount, go ahead and create a particle container for this
      // grid and add the particles to it
      auto&  particles = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

      ParticleType p_new;
      for (int i = 0; i < pcount; i++) {
        // Set id and cpu for this particle
        p_new.id()  = ParticleType::NextID();
        p_new.cpu() = ParallelDescriptor::MyProc();

        // Add to the data structure
        particles.push_back(p_new);
      }

      const int np = pcount;
      total_np += np;

      // Now define the rest of the particle data and store it directly in the particles
      // std::cout << pcount << " particles " << " in grid " << grid_id << std::endl;

      if (pcount > 0)
         mfix_particle_generator_prop(&np, particles.GetArrayOfStructs().data());
  }

  ParallelDescriptor::ReduceIntSum(total_np,ParallelDescriptor::IOProcessorNumber());
  amrex::Print() << "Total number of generated particles: " << total_np << std::endl;

  // We shouldn't need this if the particles are tiled with one tile per grid, but otherwise
  // we do need this to move particles from tile 0 to the correct tile.
  Redistribute();

}

void MFIXParticleContainer::RemoveOutOfRange(int lev, const EBFArrayBoxFactory * ebfactory,
                                             const MultiFab * ls_phi, const iMultiFab * ls_valid,
                                             int ls_refinement)
{

    // Only call the routine for wall collisions if we actually have walls
    if (ebfactory != NULL) {

        Box domain(Geom(lev).Domain());
        const Real * dx = Geom(lev).CellSize();
        MultiFab dummy;

        // amrex::Print() << ParticleBoxArray(lev) << std::endl;

        dummy.define(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                     1, 0, MFInfo(), * ebfactory);

        for (MFIXParIter pti(* this, lev); pti.isValid(); ++pti) {
            // Real particles
            const int nrp = NumberOfParticles(pti);

            void * particles  = pti.GetArrayOfStructs().data();

            const auto & sfab = static_cast <EBFArrayBox const&>((dummy)[pti]);
            const auto & flag = sfab.getEBCellFlagFab();

            const Box & bx = pti.tilebox();

            // Remove particles outside of or touching the walls
            if (flag.getType(bx) != FabType::regular)
            {
                if (flag.getType(bx) == FabType::covered)
                {
                    for (auto & p: pti.GetArrayOfStructs())
                        p.id() = -1;

                }
                else
                {
                    rm_wall_collisions_eb(particles, &nrp,
                                          BL_TO_FORTRAN_3D((*ls_valid)[pti]),
                                          BL_TO_FORTRAN_3D((*ls_phi)[pti]),
                                          BL_TO_FORTRAN_3D(flag),
                                          Geom(lev).ProbLo(),
                                          dx, & ls_refinement);
                }
            }
        }

        Redistribute();

        long fin_np = 0;
        for (MFIXParIter pti(* this, lev); pti.isValid(); ++pti) {
            long np = pti.numParticles();
            fin_np += np;
        }

        ParallelDescriptor::ReduceLongSum(fin_np,ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "Final number of particles: " << fin_np << std::endl;
    }
}

void MFIXParticleContainer::PrintParticleCounts() {
  const int lev = 0;
  amrex::AllPrintToFile("load_balance") << "Particles on each box: \n";
  long local_count = 0;
  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      long np = pti.numParticles();
      local_count += np;
      amrex::AllPrintToFile("load_balance") << "Box:" << pti.index() << ", count: " << np << std::endl;
    }
  amrex::AllPrintToFile("load_balance") << "Total for this process: " << local_count << std::endl << std::endl;
}

void MFIXParticleContainer::Replicate(IntVect& Nrep, Geometry& geom, DistributionMapping& dmap, BoxArray& ba)
{
    int lev = 0;

    Vector<Real> orig_domain_size;
    orig_domain_size.resize(BL_SPACEDIM);
    for (int d = 0; d < BL_SPACEDIM; d++)
        orig_domain_size[d] = (geom.ProbHi(d) - geom.ProbLo(d)) / Nrep[d];

    ParticleType p_rep;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();

        for (const auto& p: particles)
        {
           //
           // Shift the position.
           //
           for (int k = 0; k < Nrep[2]; k++) {
               for (int j = 0; j < Nrep[1]; j++) {
                 for (int i = 0; i < Nrep[0]; i++) {

                   if ( !(i == 0 && j == 0 && k == 0) )
                   {
                    p_rep.m_rdata.pos[0] = p.m_rdata.pos[0] + i * orig_domain_size[0];
                    p_rep.m_rdata.pos[1] = p.m_rdata.pos[1] + j * orig_domain_size[1];
                    p_rep.m_rdata.pos[2] = p.m_rdata.pos[2] + k * orig_domain_size[2];

                    p_rep.rdata(realData::velx)   = p.rdata(realData::velx);
                    p_rep.rdata(realData::vely)   = p.rdata(realData::vely);
                    p_rep.rdata(realData::velz)   = p.rdata(realData::velz);

                    // Set other particle properties
                    p_rep.idata(intData::phase)     = p.idata(intData::phase);
                    p_rep.idata(intData::state)     = p.idata(intData::state);
                    p_rep.rdata(realData::volume)   = p.rdata(realData::volume);
                    p_rep.rdata(realData::density)  = p.rdata(realData::density);
                    p_rep.rdata(realData::mass)     = p.rdata(realData::mass);
                    p_rep.rdata(realData::oneOverI) = p.rdata(realData::oneOverI);
                    p_rep.rdata(realData::radius)   = p.rdata(realData::radius);
                    p_rep.rdata(realData::omegax)   = p.rdata(realData::omegax);
                    p_rep.rdata(realData::omegay)   = p.rdata(realData::omegay);
                    p_rep.rdata(realData::omegaz)   = p.rdata(realData::omegaz);
                    p_rep.rdata(realData::dragx)    = p.rdata(realData::dragx);
                    p_rep.rdata(realData::dragy)    = p.rdata(realData::dragy);
                    p_rep.rdata(realData::dragz)    = p.rdata(realData::dragz);

                    // Set id and cpu for this particle
                    p_rep.id()  = ParticleType::NextID();
                    p_rep.cpu() = ParallelDescriptor::MyProc();

                    // Add everything to the data structure
                    particles.push_back(p_rep);

                   } // not copying itself
                 } // i
              } // j
           } // k
        } // p
    } // pti

    Redistribute();
}

void MFIXParticleContainer:: printParticles()
{
    const int lev = 0;
    const auto& plevel = GetParticles(lev);

    for (const auto& kv : plevel)
    {
       const auto& particles = kv.second.GetArrayOfStructs();

       for (unsigned i = 0; i < particles.numParticles(); ++i)
       {
          std::cout << "Particle ID  = " << i << " " << std::endl;
          std::cout << "X            = " << particles[i].pos(0) << " " << std::endl;
          std::cout << "Y            = " << particles[i].pos(1) << " " << std::endl;
          std::cout << "Z            = " << particles[i].pos(2) << " " << std::endl;
          std::cout << "state        = " << particles[i].idata(intData::state) << " " << std::endl;
          std::cout << "phase        = " << particles[i].idata(intData::phase) << " " << std::endl;
          std::cout << "Real properties = " << std::endl;

          for (int j = 0; j < realData::count; j++)
            std::cout << "property " << j << "  = " << particles[i].rdata(j) << " " << std::endl;

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

        do_tiling = true;  // because the default in amrex is false

        pp.query("do_tiling",  do_tiling);

        Vector<int> ts(BL_SPACEDIM);

        if (pp.queryarr("tile_size", ts))
            tile_size = IntVect(ts);

        pp.query("use_neighbor_list", use_neighbor_list);
        pp.query("sort_neighbor_list", sort_neighbor_list);

        initialized = true;
    }
}

void
MFIXParticleContainer::InitData()
{
}

std::unique_ptr<MultiFab> MFIXParticleContainer::EBNormals(int lev,
                                                           EBFArrayBoxFactory * ebfactory, MultiFab * dummy)
{
    // Container for normal data
    std::unique_ptr<MultiFab> normal = std::unique_ptr<MultiFab>(new MultiFab);

    // Only call the routine for wall collisions if the box has a wall
    if (ebfactory != NULL) {
        dummy->define(ParticleBoxArray(lev), ParticleDistributionMap(lev), 1, 0, MFInfo(), * ebfactory);
        std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac = ebfactory->getAreaFrac();

        // We pre-compute the normals
        normal->define(ParticleBoxArray(lev), ParticleDistributionMap(lev), 3, 2);

        for(MFIter mfi(* normal, true); mfi.isValid(); ++mfi){
            Box tile_box = mfi.tilebox();
            const int* lo = tile_box.loVect();
            const int* hi = tile_box.hiVect();

            const auto& sfab = static_cast <EBFArrayBox const&>((*dummy)[mfi]);
            const auto& flag = sfab.getEBCellFlagFab();

            if (flag.getType(amrex::grow(tile_box,1)) == FabType::singlevalued)
            {
               BL_PROFILE_VAR("compute_normals()", compute_normals);
                amrex_eb_compute_normals(lo, hi,
                                         BL_TO_FORTRAN_3D(flag),
                                         BL_TO_FORTRAN_3D((* normal)[mfi]),
                                         BL_TO_FORTRAN_3D((* areafrac[0])[mfi]),
                                         BL_TO_FORTRAN_3D((* areafrac[1])[mfi]),
                                         BL_TO_FORTRAN_3D((* areafrac[2])[mfi])   );
               BL_PROFILE_VAR_STOP(compute_normals);
            }
        }
        normal->FillBoundary(Geom(0).periodicity());
    }

    return normal;
}


void MFIXParticleContainer::EvolveParticles(int lev, int nstep, Real dt, Real time,
        EBFArrayBoxFactory * ebfactory, MultiFab * eb_normals,
        const MultiFab * ls_phi, const iMultiFab * ls_valid, const int ls_refinement,
        MultiFab * dummy, MultiFab * cost, std::string & knapsack_weight_type, int subdt_io)
{
    BL_PROFILE_REGION_START("mfix_dem::EvolveParticles()");
    BL_PROFILE("mfix_dem::EvolveParticles()");

    amrex::Print() << "Evolving particles... " << std::endl;

    /****************************************************************************
     * DEBUG flag toggles:                                                      *
     *   -> Print number of collisions                                          *
     *   -> Print max (over substeps) particle velocity at each time step       *
     *   -> Print max particle-wall and particle-particle forces                *
     ****************************************************************************/

    // Debug level controls the detail of debug outut:
    //   -> debug_level = 0 : no debug output
    //   -> debug_level = 1 : debug output for every fluid step
    //   -> debug_level = 2 : debug output for every substep
    const int debug_level = 0;

    /****************************************************************************
     * Geometry                                                                 *
     ****************************************************************************/

    Box domain(Geom(lev).Domain());

    const Real* dx = Geom(lev).CellSize();

    Real xlen = Geom(lev).ProbHi(0) - Geom(lev).ProbLo(0);
    Real ylen = Geom(lev).ProbHi(1) - Geom(lev).ProbLo(1);
    Real zlen = Geom(lev).ProbHi(2) - Geom(lev).ProbLo(2);

    /****************************************************************************
     * Init substeps                                                            *
     ***************************************************************************/

    int   nsubsteps;
    Real  subdt, stime = time;
    des_init_time_loop( &time, &dt, &nsubsteps, &subdt, &subdt_io );

    /****************************************************************************
     * Init temporary storage:                                                  *
     *   -> particle-particle, and particle-wall forces                         *
     *   -> particle-particle, and particle-wall torques                        *
     ***************************************************************************/

    std::map<PairIndex, Vector<Real>> tow;
    std::map<PairIndex, Vector<Real>> fc, pfor, wfor;
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        tow[index] = Vector<Real>();
        fc[index] = Vector<Real>();
        pfor[index] = Vector<Real>();
        wfor[index] = Vector<Real>();
    }

    /****************************************************************************
     * Iterate over sub-steps                                                   *
     ***************************************************************************/

    int ncoll_total = 0;  // Counts total number of collisions
    loc_maxvel  = RealVect(0., 0., 0.);  // Tracks max (absolute) velocity
    loc_maxpfor = RealVect(0., 0., 0.);  // Tracks max particle-particle force
    loc_maxwfor = RealVect(0., 0., 0.);  // Tracks max particle-wall force
    int n = 0; // Counts sub-steps

    while (n < nsubsteps)
    {
      int ncoll = 0;  // Counts number of collisions (over sub-steps)

      // Redistribute particles ever so often BUT always update the neighbour
      // list (Note that this fills the neighbour list after every redistribute
      // operation)
      if (n % 25 == 0) {
          if (lev == 0) clearNeighbors(lev);
          Redistribute();
          if (lev == 0) fillNeighbors(lev);
          if (lev == 0) buildNeighborList(lev, MFIXCheckPair, sort_neighbor_list);
      } else {
          if (lev == 0) updateNeighbors(lev);
      }

#ifdef _OPENMP
#pragma omp parallel reduction(+:ncoll)
#endif
      for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
      {

        /***********************************************************************
         * Particle-tile's book-keeping                                        *
         ***********************************************************************/

         // Timer used for load-balancing
         Real wt = ParallelDescriptor::second();

         // Real particles
         const int nrp    = NumberOfParticles(pti);
         RealType* particles  = pti.GetArrayOfStructs().data();

         // Neighbor particles
         PairIndex index(pti.index(), pti.LocalTileIndex());
         int size_ng = neighbors[lev][index].size();
         int size_nl = neighbor_list[lev][index].size();

         // Number of particles including neighbor particles
         int ntot = nrp + size_ng;

         // Particle tile box: used to check if the tile-box contains walls =>
         // if not, then don't check for wall collisions.
         const Box& bx = pti.tilebox();

         // Particle-particle (and particle-wall) forces and torques. We need
         // these to be zero every time we start a new batch (i.e tile and
         // substep) of particles.
         tow[index].clear();
          fc[index].clear();
         tow[index].resize(ntot*3,0.0);
          fc[index].resize(ntot*3,0.0);

         // For debugging: keep track of particle-particle (pfor) and
         // particle-wall (wfor) forces
         pfor[index].clear();
         wfor[index].clear();
         pfor[index].resize(3 * ntot, 0.0);
         wfor[index].resize(3 * ntot, 0.0);

         /***********************************************************************
         * Particle-Wall collision forces (and torques)                        *
         ***********************************************************************/

         // Only call the routine for wall collisions if we actually have walls
         if (ebfactory != NULL)
         {
             // Get particle EB geometric info (the inputted `dummy` is defined
             // on the wrong grids)
             MultiFab  dummy(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                             1, 0, MFInfo(), * ebfactory);

             const auto & sfab = static_cast <EBFArrayBox const &>(dummy[pti]);
             const auto & flag = sfab.getEBCellFlagFab();

             if (flag.getType(amrex::grow(bx,1)) == FabType::singlevalued)
             {
                 const MultiCutFab * bndrycent = &(ebfactory->getBndryCent());

                 // Calculate forces and torques from particle-wall collisions
                 BL_PROFILE_VAR("calc_wall_collisions()", calc_wall_collisions);
                 if(legacy__eb_collisions) {
                     calc_wall_collisions(particles, & ntot, & nrp,
                                          tow[index].dataPtr(), fc[index].dataPtr(), & subdt,
                                          BL_TO_FORTRAN_3D(flag),
                                          BL_TO_FORTRAN_3D((* eb_normals)[pti]),
                                          BL_TO_FORTRAN_3D((* bndrycent)[pti]),
                                          dx);
                 } else {
                     calc_wall_collisions_ls(particles, & ntot, & nrp,
                                             tow[index].dataPtr(), fc[index].dataPtr(), & subdt,
                                             BL_TO_FORTRAN_3D((* ls_valid)[pti]),
                                             BL_TO_FORTRAN_3D((* ls_phi)[pti]),
                                             dx, & ls_refinement);
                 }

                 // Debugging: copy data from the fc (all forces) vector to the
                 // wfor (wall forces) vector.
                 if (debug_level > 0) {
                     for (int i = 0; i < wfor[index].size(); i++ ) {
                         wfor[index][i] = fc[index][i];
                     }
                 }

                 BL_PROFILE_VAR_STOP(calc_wall_collisions);
             }
         }

        /***********************************************************************
        * Particle-Particle collision forces (and torques)                    *
        ***********************************************************************/

#if 1
         BL_PROFILE_VAR("calc_particle_collisions()", calc_particle_collisions);

         if (lev == 0)
         calc_particle_collisions ( particles                     , &nrp,
                                    neighbors[lev][index].dataPtr()    , &size_ng,
                                    neighbor_list[lev][index].dataPtr(), &size_nl,
                                    tow[index].dataPtr(), fc[index].dataPtr(),
                                    &subdt, &ncoll);

         // Debugging: copy data from the fc (all forces) vector to the wfor
         // (wall forces) vector. Note that since fc already contains the wall
         // forces, these need to be subtracted here.
         if (debug_level > 0) {
             for (int i = 0; i < pfor[index].size(); i++ ) {
                 pfor[index][i] = fc[index][i] - wfor[index][i];
             }
         }

         BL_PROFILE_VAR_STOP(calc_particle_collisions);
#else
         Vector<Real> x(nrp);
         Vector<Real> y(nrp);
         Vector<Real> z(nrp);
         particle_get_position (particles, nrp, x, y, z);

         BL_PROFILE_VAR("calc_particle_collisions()", calc_particle_collisions);
         calc_particle_collisions_soa ( particles                     , &nrp,
                                        neighbors[lev][index].dataPtr()    , &size_ng,
                                        neighbor_list[lev][index].dataPtr(), &size_nl,
                                        tow[index].dataPtr(), fc[index].dataPtr(), &subdt, &ncoll);
         BL_PROFILE_VAR_STOP(calc_particle_collisions);
#endif

        /***********************************************************************
         * Move particles based on collision forces and torques                *
         ***********************************************************************/

        BL_PROFILE_VAR("des_time_loop()", des_time_loop);
#if 1
        des_time_loop ( &nrp     , particles,
                        &ntot, tow[index].dataPtr(), fc[index].dataPtr(), &subdt,
                        &xlen, &ylen, &zlen, &stime, &n);
#else
        des_time_loop_soa ( &nrp, particles,
                            &ntot, tow[index].dataPtr(), fc[index].dataPtr(), &subdt,
                            &xlen, &ylen, &zlen, &stime, &n);
#endif
        BL_PROFILE_VAR_STOP(des_time_loop);

        /***********************************************************************
        * Update runtime cost (used in load-balancing)                        *
        ***********************************************************************/

        //if (lev == 0)
        //if (cost) {
        //     // Runtime cost is either (weighted by tile box size):
        //     //   * time spent
        //     //   * number of particles
        //     const Box& tbx = pti.tilebox();
        //     if (knapsack_weight_type == "RunTimeCosts")
        //     {
        //        wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
        //     }
        //     else if (knapsack_weight_type == "NumParticles")
        //     {
        //        wt = nrp / tbx.d_numPts();
        //     }
        //     (*cost)[pti].plus(wt, tbx);
        // }
      }

      // Update substep count
      n += 1;

     /**************************************************************************
      * DEBUG: output the number of collisions in current substep              *
      *        output the max velocity (and forces) in current substep         *
      *        update max velocities and forces                                *
      **************************************************************************/

      if (debug_level > 1) {
         ncoll_total += ncoll;
         ParallelDescriptor::ReduceIntSum(ncoll, ParallelDescriptor::IOProcessorNumber());
         Print() << "Number of collisions: " << ncoll << " at step " << n << std::endl;
      }

      if (debug_level > 0){
          UpdateMaxVelocity();
          UpdateMaxForces(pfor, wfor);
      }

      if (debug_level > 1) {
          RealVect max_vel = GetMaxVelocity();
          Vector<RealVect> max_forces = GetMaxForces();

          const Real * dx_crse = Geom(0).CellSize();
          amrex::Print() << "Maximum distance traveled:"
                         << std::endl
                         <<  "x= " << max_vel[0] * dt
                         << " y= " << max_vel[1] * dt
                         << " z= " << max_vel[2] * dt
                         << " and note that "
                         << " dx= " << dx_crse[0] << std::endl;

          amrex::Print() << "Maximum particle-particle (pp) and particle-wall (pw) forces:"
                         << std::endl
                         <<  "ppx= " << max_forces[0][0]
                         << " ppy= " << max_forces[0][1]
                         << " ppz= " << max_forces[0][2] << std::endl
                         <<  "pwx= " << max_forces[1][0]
                         << " pwy= " << max_forces[1][1]
                         << " pwz= " << max_forces[1][2] << std::endl;

      }

    } // end of loop over substeps

    // Redistribute particles at the end of all substeps (note that the
    // particle neighbour list needs to be reset when redistributing).
    if (lev == 0) clearNeighbors(lev);
    Redistribute();

   /****************************************************************************
    * DEBUG: output the total number of collisions over all substeps           *
    *        output the maximum velocity and forces over all substeps          *
    ****************************************************************************/
    if (debug_level > 0) {
       ParallelDescriptor::ReduceIntSum(ncoll_total, ParallelDescriptor::IOProcessorNumber());
       Print() << "Number of collisions: " << ncoll_total << " in " << nsubsteps << " substeps " << std::endl;
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

      const int nrp   = NumberOfParticles(pti);
      void* particles = pti.GetArrayOfStructs().data();

      call_usr3_des( &nrp, particles );
    }

    if (debug_level > 0) {
        RealVect max_vel = GetMaxVelocity();
        Vector<RealVect> max_forces = GetMaxForces();

        const Real * dx_crse = Geom(0).CellSize();
        amrex::Print() << "Maximum possible distance traveled:" << std::endl
                       <<  "x= " << max_vel[0] * dt
                       << " y= " << max_vel[1] * dt
                       << " z= " << max_vel[2] * dt
                       << " and note that "
                       << " dx= " << dx_crse[0] << std::endl;

        amrex::Print() << "Maximum particle-particle (pp) and particle-wall (pw) forces:" << std::endl
                       <<  "ppx= " << max_forces[0][0]
                       << " ppy= " << max_forces[0][1]
                       << " ppz= " << max_forces[0][2] << std::endl
                       <<  "pwx= " << max_forces[1][0]
                       << " pwy= " << max_forces[1][1]
                       << " pwz= " << max_forces[1][2] << std::endl;
    }

    amrex::Print() << "done. \n";

    BL_PROFILE_REGION_STOP("mfix_dem::EvolveParticles()");
}

void MFIXParticleContainer::CalcVolumeFraction(amrex::MultiFab& mf_to_be_filled,
                                               const EBFArrayBoxFactory& ebfactory,
                                               IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                                               IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                                               IArrayBox& bc_klo, IArrayBox& bc_khi,
                                               int nghost )
{

    // NOTE: ebfactory HAS to be the PARTICLE EB factory!
    int fortran_volume_comp = 5;
    PICDeposition(mf_to_be_filled, ebfactory,
                  bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo,bc_khi,
                  fortran_volume_comp,nghost);

    // Now define this mf = (1 - particle_vol)
    mf_to_be_filled.mult(-1.0,mf_to_be_filled.nGrow());
    mf_to_be_filled.plus( 1.0,mf_to_be_filled.nGrow());

    // Impose a lower bound on volume fraction
    CapSolidsVolFrac(mf_to_be_filled);

}

void MFIXParticleContainer::CalcDragOnFluid(amrex::MultiFab& beta_mf,
                                            amrex::MultiFab& beta_vel_mf,
                                            const EBFArrayBoxFactory& ebfactory,
                                            IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                                            IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                                            IArrayBox& bc_klo, IArrayBox& bc_khi,
                                            int nghost )
{
    int fortran_beta_comp = 15;
    int fortran_vel_comp  =  9;
    PICMultiDeposition(beta_mf, beta_vel_mf, ebfactory,
                       bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo,bc_khi,
                       fortran_beta_comp, fortran_vel_comp, nghost);
}

void MFIXParticleContainer::PICDeposition(amrex::MultiFab& mf_to_be_filled,
                                          const EBFArrayBoxFactory& ebfactory,
                                          IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                                          IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                                          IArrayBox& bc_klo, IArrayBox& bc_khi,
                                          int fortran_particle_comp, int nghost )
{
    BL_PROFILE("MFIXParticleContainer::PICDeposition()");

    int   lev = 0;
    int ncomp = 1;

    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf_to_be_filled)) {
      // If we are already working with the internal mf defined on the
      // particle_box_array, then we just work with this.
      mf_pointer = &mf_to_be_filled;
    }
    else {
      // If mf_to_be_filled is not defined on the particle_box_array, then we need
      // to make a temporary here and copy into mf_to_be_filled at the end.
      mf_pointer = new MultiFab(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                                ncomp, mf_to_be_filled.nGrow());
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (mf_pointer->nGrow() < 1)
       amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx          = gm.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*mf_pointer, true); mfi.isValid(); ++mfi) {
        (*mf_pointer)[mfi].setVal(0);
    }

    using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;

    // Get particle EB geometric info
    MultiFab      dummy(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                        1, 0, MFInfo(), ebfactory);

    const amrex::MultiFab*                    volfrac;
    volfrac = &(ebfactory.getVolFrac());

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox local_vol;
        for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {
            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long nrp = pti.numParticles();
            FArrayBox& fab = (*mf_pointer)[pti];
            Real* data_ptr;
            const int *lo, *hi;
#ifdef _OPENMP
            Box tile_box = pti.tilebox();
            tile_box.grow(1);
            local_vol.resize(tile_box,ncomp);
            local_vol = 0.0;
            data_ptr = local_vol.dataPtr();
            lo = tile_box.loVect();
            hi = tile_box.hiVect();
#else
            data_ptr = fab.dataPtr();
            const Box& box = fab.box();
            lo = box.loVect();
            hi = box.hiVect();
#endif

            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox&  dummy_fab = static_cast<EBFArrayBox const&>(dummy[pti]);
            const EBCellFlagFab&    flags = dummy_fab.getEBCellFlagFab();

            const Box& bx  = pti.tilebox(); // I need a box without ghosts

            if (flags.getType(bx) != FabType::covered )
            {
                mfix_deposit_cic_eb(particles.data(), nstride, nrp, ncomp, data_ptr,
                                    lo, hi,
                                    BL_TO_FORTRAN_ANYD((*volfrac)[pti]),
                                    BL_TO_FORTRAN_ANYD(flags),
                                    plo, dx, &fortran_particle_comp);
            }


#ifdef _OPENMP
            amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_vol),
                                        BL_TO_FORTRAN_3D(fab), ncomp);
#endif

        }
    }

    // Move any field deposited outside the domain back into the domain
    // when BC is pressure inlet and mass inflow.
    Box domain(Geom(lev).Domain());

    for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {

      const Box& sbx = (*mf_pointer)[mfi].box();

      flip_particle_vol(sbx.loVect(), sbx.hiVect(),
                        (*mf_pointer)[mfi].dataPtr(),
                        bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                        bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                        bc_klo.dataPtr(), bc_khi.dataPtr(),
                        domain.loVect(), domain.hiVect(),
                        &nghost );
    }

    mf_pointer->SumBoundary(gm.periodicity());

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled) {
      mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
      delete mf_pointer;
    }

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::PICDeposition time: " << stoptime << '\n';
    }
}

void MFIXParticleContainer::PICMultiDeposition(amrex::MultiFab& beta_mf,
                                               amrex::MultiFab& beta_vel_mf,
                                               const EBFArrayBoxFactory& ebfactory,
                                               IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                                               IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                                               IArrayBox& bc_klo, IArrayBox& bc_khi,
                                               int fortran_beta_comp, int fortran_vel_comp, int nghost )
{
    int lev = 0;
    AMREX_ASSERT(OnSameGrids(lev,beta_mf)==OnSameGrids(lev,beta_vel_mf));

    BL_PROFILE("MFIXParticleContainer::PICMultiDeposition()");

    MultiFab*  beta_ptr;
    MultiFab*  beta_vel_ptr;

    if (OnSameGrids(lev,beta_mf))
    {
      beta_ptr     = &beta_mf;
      beta_vel_ptr = &beta_vel_mf;
    }
    else
    {
        // Make temporaries here and copy into beta_mf and beta_vel_mf at the end.
        beta_ptr     = new MultiFab(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                                    beta_mf.nComp(), beta_mf.nGrow());
        beta_vel_ptr = new MultiFab(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                                    beta_vel_mf.nComp(), beta_vel_mf.nGrow());
    }

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx          = gm.CellSize();

    beta_ptr->setVal(0.0);
    beta_vel_ptr->setVal(0.0);

    using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;
    const int* lo;
    const int* hi;
    Real* bx_dataptr;
    Real* bu_dataptr;

    // Get particle EB geometric info
    MultiFab  dummy(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                    1, 0, MFInfo(), ebfactory);

    const amrex::MultiFab*   volfrac;
    volfrac = &(ebfactory.getVolFrac());

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox local_x_vol, local_u_vol;
        for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {

            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np = pti.numParticles();

            FArrayBox& beta_fab = (*beta_ptr)[pti];
            FArrayBox& beta_vel_fab = (*beta_vel_ptr)[pti];

#ifdef _OPENMP
            // Note that we actually grow the tilebox rather than calling growntilebox
            //     because we need the overlap even in the interior.
            Box grown_tilebox = pti.tilebox();
            grown_tilebox.grow(1);

            int ncomp = BL_SPACEDIM;
            local_x_vol.resize(grown_tilebox,1);
            local_u_vol.resize(grown_tilebox,ncomp);

            local_x_vol = 0.0;
            local_u_vol = 0.0;

            bx_dataptr = local_x_vol.dataPtr();
            bu_dataptr = local_u_vol.dataPtr();

            lo = grown_tilebox.loVect();
            hi = grown_tilebox.hiVect();
#else
            bx_dataptr = beta_fab.dataPtr();
            bu_dataptr = beta_vel_fab.dataPtr();

            const Box& bx  = beta_fab.box();

            lo = bx.loVect();
            hi = bx.hiVect();
#endif

            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox&  dummy_fab = static_cast<EBFArrayBox const&>(dummy[pti]);
            const EBCellFlagFab&    flags = dummy_fab.getEBCellFlagFab();

            const Box& box = pti.tilebox(); // I need a box without ghosts

            if (flags.getType(box) != FabType::covered )
            {
                mfix_multi_deposit_cic_eb(particles.data(), nstride, np,
                                          bx_dataptr, bu_dataptr,
                                          lo, hi,
                                          BL_TO_FORTRAN_ANYD((*volfrac)[pti]),
                                          plo, dx, &fortran_beta_comp,
                                          &fortran_vel_comp);
           }

#ifdef _OPENMP
            amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_x_vol),
                                        BL_TO_FORTRAN_3D(beta_fab), beta_fab.nComp());
            amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_u_vol),
                                        BL_TO_FORTRAN_3D(beta_vel_fab), beta_vel_fab.nComp());
#endif

        }
    }

    beta_ptr->SumBoundary(gm.periodicity());
    beta_vel_ptr->SumBoundary(gm.periodicity());

    if ( &beta_mf != beta_ptr)
    {
        // Copy back from mf_pointer
        beta_mf.copy    (*beta_ptr,0,0,beta_mf.nComp());
        beta_vel_mf.copy(*beta_vel_ptr,0,0,beta_vel_mf.nComp());

        delete beta_ptr;
        delete beta_vel_ptr;
    }

    const Box domain(Geom(lev).Domain());

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::PICMultiDeposition time: " << stoptime << '\n';
    }
}

void MFIXParticleContainer::output(int estatus, int finish, int nstep, Real dt, Real time)
{
    for (int lev = 0; lev < nlev; lev++)
    {
       Real xlen = Geom(lev).ProbHi(0) - Geom(lev).ProbLo(0);
       Real ylen = Geom(lev).ProbHi(1) - Geom(lev).ProbLo(1);
       Real zlen = Geom(lev).ProbHi(2) - Geom(lev).ProbLo(2);

       // Not threaded because its writing output
       for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
       {

         //number of particles
         const int     np = NumberOfParticles(pti);
         void* particles  = pti.GetArrayOfStructs().data();

         output_manager( &np, &time, &dt, &xlen, &ylen, &zlen, &nstep,
                              particles, &finish);
       }
    }
}

void MFIXParticleContainer::writeAllAtLevel(int lev)
{
    // Not threaded because its print to terminal
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();

        for (const auto& p: particles)
        {
           const IntVect& iv = Index(p, lev);

           RealVect xyz(p.pos(0), p.pos(1), p.pos(2));
           cout << " id " << p.id()
                << " index " << iv
                << " position " << xyz << endl;
       }
    }
}

void
MFIXParticleContainer::writeAllForComparison(int lev)
{
  int Np_tot = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:Np_tot)
#endif
  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
      Np_tot += pti.numParticles();

  ParallelDescriptor::ReduceIntSum(Np_tot,ParallelDescriptor::IOProcessorNumber());

  cout << Np_tot << std::endl;

  // Not threaded because its print to terminal
  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
  {
      auto& particles = pti.GetArrayOfStructs();

      for (const auto& p: particles)
         cout << p.pos(0)   << " " << p.pos(1)   << " " << p.pos(2) <<  " " <<
                 p.rdata(1) << " " << p.rdata(2) << " " << p.rdata(2) <<  " " <<
                 ParallelDescriptor::MyProc() << " " << p.id() << std::endl;
  }
}

void
MFIXParticleContainer::WriteAsciiFileForInit (const std::string& filename)
{
    BL_ASSERT(!filename.empty());

    int lev = 0;
    long nparticles = NumberOfParticlesAtLevel (lev);

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Have I/O processor open file and write out particle metadata.
        //
        std::ofstream File;

        File.open(filename.c_str(), std::ios::out|std::ios::trunc);

        if (!File.good())
            amrex::FileOpenFailed(filename);

        File << nparticles  << '\n';

        File.flush();

        File.close();

        if (!File.good())
            amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::WriteAsciiFile(): problem writing file");
    }

    ParallelDescriptor::Barrier();

    const int MyProc = ParallelDescriptor::MyProc();

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (MyProc == i)
        {
            //
            // Each CPU opens the file for appending and adds its particles.
            //

            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            std::ofstream File;

            File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            File.open(filename.c_str(), std::ios::out|std::ios::app);

            File.precision(15);

            if (!File.good())
                amrex::FileOpenFailed(filename);

            for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

              auto& particles = pti.GetArrayOfStructs();

    int index = 0;
                for (const auto& p: particles)
                {
        if (p.id() > 0) {
                        File << p.idata(intData::phase) << ' ';
                        File << p.pos(0) << ' ';
                        File << p.pos(1) << ' ';
                        File << p.pos(2) << ' ';
                        File << p.rdata(realData::radius) << ' ';
                        File << p.rdata(realData::density) << ' ';
                        File << p.rdata(realData::velx) << ' ';
                        File << p.rdata(realData::vely) << ' ';
                        File << p.rdata(realData::velz) << ' ';

                        File << '\n';

                        index++;
                    }
                }
              }

            File.flush();

            File.close();

            if (!File.good())
                amrex::Abort("MFIXParticleContainer::WriteAsciiFileForInit(): problem writing file");

        }
        ParallelDescriptor::Barrier();
    }
}

void MFIXParticleContainer::GetParticleAvgProp(Real (&avg_dp)[10], Real (&avg_ro)[10])
{
   // The number of phases was previously hard set at 10, however lowering
   //  this number would make this code faster.
   int num_of_phases_in_use = 10; //Number of different phases being simulated

   // Cycle through the different phases, starting from 1
   for (int phse = 1; phse <= num_of_phases_in_use; ++phse)
   {
     Real p_num  = 0.0; //number of particle
     Real p_diam = 0.0; //particle diameters
     Real p_dens = 0.0; //particle density

     for (int lev = 0; lev < nlev; lev++)
     {
#ifdef _OPENMP
#pragma omp parallel reduction(+:p_num, p_diam, p_dens)
#endif
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
          auto& particles = pti.GetArrayOfStructs();

          for (const auto& p: particles){
            if ( phse==p.idata(intData::phase) ){
              p_num  += 1.0;
              p_diam += p.rdata(realData::radius) * 2.0;
              p_dens += p.rdata(realData::density);
            }
          }
        }
     }

     // A single MPI call passes all three variables
     ParallelDescriptor::ReduceRealSum({p_num,p_diam,p_dens});

     //calculate averages or set = zero if no particles of that phase
     if (p_num==0){
       avg_dp[phse-1] = 0.0;
       avg_ro[phse-1] = 0.0;
     } else {
       avg_dp[phse-1] = p_diam/p_num;
       avg_ro[phse-1] = p_dens/p_num;
     }
   }
}

void MFIXParticleContainer::UpdateMaxVelocity()
{
    Real max_vel_x = loc_maxvel[0], max_vel_y = loc_maxvel[1], max_vel_z = loc_maxvel[2];

#ifdef _OPENMP
#pragma omp parallel reduction(max:max_vel_x,max_vel_y,max_vel_z)
#endif
    for (int lev = 0; lev < nlev; lev++)
    {
       for(MFIXParIter pti(* this, lev); pti.isValid(); ++ pti)
       {
           auto & particles = pti.GetArrayOfStructs();
           for(const auto & particle : particles)
           {
              max_vel_x = std::max(Real(std::fabs(particle.rdata(realData::velx))), max_vel_x);
              max_vel_y = std::max(Real(std::fabs(particle.rdata(realData::vely))), max_vel_y);
              max_vel_z = std::max(Real(std::fabs(particle.rdata(realData::velz))), max_vel_z);
           }
       }
    }
    loc_maxvel = RealVect(max_vel_x, max_vel_y, max_vel_z);
}

void MFIXParticleContainer::UpdateMaxForces( std::map<PairIndex, Vector<Real>> pfor,
                                             std::map<PairIndex, Vector<Real>> wfor)
{
    Real max_pfor_x = loc_maxpfor[0], max_pfor_y = loc_maxpfor[1], max_pfor_z = loc_maxpfor[2];
    Real max_wfor_x = loc_maxwfor[0], max_wfor_y = loc_maxwfor[1], max_wfor_z = loc_maxwfor[2];

    for (int lev = 0; lev < nlev; lev++)
    {
#ifdef _OPENMP
#pragma omp parallel reduction(max:max_pfor_x,max_pfor_y,max_pfor_z,max_wfor_x,max_wfor_y,max_wfor_z)
#endif
       for(MFIXParIter pti(* this, lev); pti.isValid(); ++ pti)
       {
        PairIndex index(pti.index(), pti.LocalTileIndex());

        // Note the particle force data layout:
        //      p1_x, p2_x, ..., pn_x, p1_y, p2_y, ..., pn_y, p1_z, p2_z, ..., pn_z
        // Where n is the total number of particle and neighbor particles.
        const int nrp     = NumberOfParticles(pti);
        const int size_ng = neighbors[lev][index].size();
        // Number of particles including neighbor particles
        const int ntot = nrp + size_ng;

        // Find max (abs) of particle-particle forces:
        for(int i = 0; i < ntot; i++ )
            max_pfor_x = std::max(Real(std::fabs(pfor[index][i])), max_pfor_x);
        for(int i = ntot; i < 2 * ntot; i++ )
            max_pfor_y = std::max(Real(std::fabs(pfor[index][i])), max_pfor_y);
        for(int i = 2 * ntot; i < 3 * ntot; i++ )
            max_pfor_z = std::max(Real(std::fabs(pfor[index][i])), max_pfor_z);

        // Find max (abs) of particle-wall forces:
        for(int i = 0; i < ntot; i++ )
            max_wfor_x = std::max(Real(std::fabs(wfor[index][i])), max_wfor_x);
        for(int i = ntot; i < 2 * ntot; i++ )
            max_wfor_y = std::max(Real(std::fabs(wfor[index][i])), max_wfor_y);
        for(int i = 2 * ntot; i < 3 * ntot; i++ )
            max_wfor_z = std::max(Real(std::fabs(wfor[index][i])), max_wfor_z);
       }
    }

    loc_maxpfor = RealVect(max_pfor_x, max_pfor_y, max_pfor_z);
    loc_maxwfor = RealVect(max_wfor_x, max_wfor_y, max_wfor_z);
}

RealVect MFIXParticleContainer::GetMaxVelocity()
{
    Real max_vel_x = loc_maxvel[0], max_vel_y = loc_maxvel[1], max_vel_z = loc_maxvel[2];

    ParallelDescriptor::ReduceRealMax({max_vel_x, max_vel_y, max_vel_z},
                                      ParallelDescriptor::IOProcessorNumber());

    RealVect max_vel(max_vel_x, max_vel_y, max_vel_z);

    return max_vel;
};

Vector<RealVect> MFIXParticleContainer::GetMaxForces()
{
    Real max_pfor_x = loc_maxpfor[0], max_pfor_y = loc_maxpfor[1], max_pfor_z = loc_maxpfor[2];
    Real max_wfor_x = loc_maxwfor[0], max_wfor_y = loc_maxwfor[1], max_wfor_z = loc_maxwfor[2];


    ParallelDescriptor::ReduceRealMax({ max_pfor_x, max_pfor_y, max_pfor_z,
                                        max_wfor_x, max_wfor_y, max_wfor_z      },
                                      ParallelDescriptor::IOProcessorNumber());

    Vector<RealVect> max_forces(2);
    max_forces[0] = RealVect(max_pfor_x, max_pfor_y, max_pfor_z);
    max_forces[1] = RealVect(max_wfor_x, max_wfor_y, max_wfor_z);

    return max_forces;
}

void
MFIXParticleContainer::BalanceParticleLoad_KDTree()
{
  int lev = 0;
  bool verbose = true;
  BoxArray old_ba = ParticleBoxArray(lev);

  if (NumberOfParticlesAtLevel(lev) == 0)
  {
     amrex::Print() << "No particles so can't use KDTree approach " << std::endl;
     return;
  }

  if (verbose)
  {
     Vector<long> num_part;
     num_part = NumberOfParticlesInGrid(0);
     long min_number = num_part[0];
     long max_number = num_part[0];
     for (int i = 0; i < old_ba.size(); i++)
     {
        max_number = std::max(max_number, num_part[i]);
        min_number = std::min(min_number, num_part[i]);
     }
     amrex::Print() << "Before KDTree: BA had " << old_ba.size() << " GRIDS " << std::endl;
     amrex::Print() << "Before KDTree: MIN/MAX NUMBER OF PARTICLES PER GRID  " <<
                        min_number << " " << max_number << std::endl;
  }

  Vector<Real> box_costs;

  BoxArray new_ba;
  Real cell_weight = 0.;
  loadBalanceKD::balance<MFIXParticleContainer>(*this, new_ba, ParallelDescriptor::NProcs(), cell_weight, box_costs);

  // Create a new DM to go with the new BA
  DistributionMapping new_dm = DistributionMapping::makeKnapSack(box_costs);

  Regrid(new_dm, new_ba);

  if (verbose)
  {
     Vector<long> num_part;
     num_part = NumberOfParticlesInGrid(0);
     long min_number = num_part[0];
     long max_number = num_part[0];
     for (int i = 0; i < new_ba.size(); i++)
     {
        max_number = std::max(max_number, num_part[i]);
        min_number = std::min(min_number, num_part[i]);
     }
     amrex::Print() << "After  KDTree: BA had " << new_ba.size() << " GRIDS " << std::endl;
     amrex::Print() << "After  KDTree: MIN/MAX NUMBER OF PARTICLES PER GRID  " <<
                        min_number << " " << max_number << std::endl;
  }
}

void MFIXParticleContainer::ComputeAverageVelocities ( const int lev,
                   const amrex::Real time,
                   const string&  basename,
                   const vector<Real>& avg_region_x_w,
                   const vector<Real>& avg_region_x_e,
                   const vector<Real>& avg_region_y_s,
                   const vector<Real>& avg_region_y_n,
                   const vector<Real>& avg_region_z_b,
                   const vector<Real>& avg_region_z_t )
{

    // Count number of calls -- Used to determin when to create file from scratch
    static int ncalls = 0;
    ++ncalls;

    int  nregions = avg_region_x_w.size();

    //
    // Check the regions are defined correctly
    //
    if (  ( avg_region_x_e.size() != nregions ) ||
    ( avg_region_y_s.size() != nregions ) ||
    ( avg_region_y_n.size() != nregions ) ||
    ( avg_region_z_b.size() != nregions ) ||
    ( avg_region_z_t.size() != nregions )  )
    {
  amrex::Print () << "ComputeAverageVelocities: some regions are not properly defined: skipping.";
  return;
    }

    vector<long> region_np (nregions, 0);
    vector<Real> region_velx (nregions, 0.0);
    vector<Real> region_vely (nregions, 0.0);
    vector<Real> region_velz (nregions, 0.0);

    for ( int nr = 0; nr < nregions; ++nr )
    {
  // Create Real box for this region
  RealBox avg_region ( {AMREX_D_DECL(avg_region_x_w[nr],avg_region_y_s[nr],avg_region_z_b[nr])},
           {AMREX_D_DECL(avg_region_x_e[nr],avg_region_y_n[nr],avg_region_z_t[nr])} );

  // Jump to next iteration if this averaging region is not valid
  if ( !avg_region.ok () )
  {
      amrex::Print() << "ComputeAverageVelocities: region "<< nr <<" is invalid: skipping\n";
      continue;
  }

  long sum_np     = 0;    // Number of particle in avg region
  Real sum_velx   = 0.;
  Real sum_vely   = 0.;
  Real sum_velz   = 0.;

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum_np,sum_velx,sum_vely,sum_velz)
#endif
  for ( MFIXParIter pti(*this, lev); pti.isValid(); ++ pti)
  {
      Box bx       = pti.tilebox ();
      RealBox tile_region ( bx, Geom(lev).CellSize (), Geom(lev).ProbLo() );

      if ( tile_region.intersects ( avg_region ) )
      {
    const int np         = NumberOfParticles(pti);
    const AoS &particles = pti.GetArrayOfStructs();

    for (int p = 0; p < np; ++p )
    {
        if ( avg_region.contains ( &(particles[p].m_rdata.pos[0]) ) )
        {
      sum_np++;
      sum_velx += particles[p].rdata(realData::velx);
      sum_vely += particles[p].rdata(realData::vely);
      sum_velz += particles[p].rdata(realData::velz);
        }

    }

      }

  }

  region_np[nr]    = sum_np;
  region_velx[nr]  = sum_velx;
  region_vely[nr]  = sum_vely;
  region_velz[nr]  = sum_velz;
    }

    // Compute parallel reductions
    ParallelDescriptor::ReduceLongSum ( region_np.data(),   nregions );
    ParallelDescriptor::ReduceRealSum ( region_velx.data(), nregions );
    ParallelDescriptor::ReduceRealSum ( region_vely.data(), nregions );
    ParallelDescriptor::ReduceRealSum ( region_velz.data(), nregions );

    // Only the IO processor takes care of the output
    if (ParallelDescriptor::IOProcessor())
    {

  for ( int nr = 0; nr < nregions; ++nr )
  {

      //
      // Compute averages (NaN if NP=0 )
      //
      region_velx[nr] /= region_np[nr];
      region_vely[nr] /= region_np[nr];
      region_velz[nr] /= region_np[nr];

      //
      // Print to file
      //
      std::ofstream  ofs;
      std::string    fname;

      fname = basename + std::to_string(nr) + ".dat";

      // Open file
      if ( ncalls == 1 )
      {
    // Create output files only the first time this function is called
    // Use ios:trunc to delete previous contect
    ofs.open ( fname.c_str(), ios::out | ios::trunc );
      }
      else
      {
    // If this is not the first time we write to this file
    // we append to it
    ofs.open ( fname.c_str(), ios::out | ios::app );
      }

      // Check if file is good
      if ( !ofs.good() )
    amrex::FileOpenFailed ( fname );

      // Print header if first access
      if ( ncalls == 1 )
    ofs << "#  Time   NP  U  V  W" << std::endl;

      ofs << time << " "
    << region_np[nr] << " "
    << region_velx[nr] << " "
    << region_vely[nr] << " "
    << region_velz[nr] << std::endl;

      ofs.close();
  }
    }
}

void MFIXParticleContainer::CapSolidsVolFrac(amrex::MultiFab& mf_to_be_filled)
{
    int   lev = 0;

    MultiFab* eps = &mf_to_be_filled;

    Box domain(Geom(lev).Domain());

    for (MFIter mfi(*eps); mfi.isValid(); ++mfi) {
      const Box& sbx = (*eps)[mfi].box();

      mfix_cap_eps(sbx.loVect(), sbx.hiVect(), (*eps)[mfi].dataPtr());
    }

}
