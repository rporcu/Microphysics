#include <AMReX.H>
#include "AMReX_Particles.H"
#include "AMReX_RealVect.H"
#include <iostream>
#include <MFIXParticleContainer.H>
#include <AMReX_LoadBalanceKD.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB_F.H>

#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBSupport.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

#include <math.h>

#include "mfix_F.H"
#include "mfix_des_F.H"
#include "mfix_eb_F.H"
#include "mfix_util_F.H"
#include "mfix_des_K.H"
#include "MFIX_DEM_Parms.H"

using namespace amrex;
using namespace std;

bool MFIXParticleContainer::use_neighbor_list  {true};
bool MFIXParticleContainer::sort_neighbor_list {false};
Real MFIXParticleContainer::gravity[3] {0.0};

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

void MFIXParticleContainer::InitParticlesAscii(const std::string& file)
{

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
      set_particle_properties( pstate, pradius, pdensity, pvolume, pmass, pomoi, pomega);

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

  Real dx = Geom(lev).CellSize(0);
  Real dy = Geom(lev).CellSize(1);
  Real dz = Geom(lev).CellSize(2);

  int total_np = 0;

  // This uses the particle tile size. Note that the default is to tile so if we
  //      remove the true and don't explicitly add false it will still tile
  for (MFIter mfi = MakeMFIter(lev,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

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
        Gpu::HostVector<ParticleType> host_particles;
        Cuda::thrust_copy(particles.begin(), particles.end(), host_particles.begin());

        for (const auto& p: host_particles)
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

    get_gravity(gravity);

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

void MFIXParticleContainer::EvolveParticles(int lev, int nstep, Real dt, Real time,
                                            EBFArrayBoxFactory * ebfactory,
                                            const MultiFab * ls_phi, const iMultiFab * ls_valid,
                                            const int ls_refinement,
                                            MultiFab * cost, std::string & knapsack_weight_type)
{
    BL_PROFILE_REGION_START("mfix_dem::EvolveParticles()");
    BL_PROFILE("mfix_dem::EvolveParticles()");

    amrex::Print() << "Evolving particles on level: " << lev << " ... " << std::endl;

    /****************************************************************************
     * DEBUG flag toggles:                                                      *
     *   -> Print number of collisions                                          *
     *   -> Print max (over substeps) particle velocity at each time step       *
     *   -> Print max particle-wall and particle-particle forces                *
     ***************************************************************************/

    // Debug level controls the detail of debug outut:
    //   -> debug_level = 0 : no debug output
    //   -> debug_level = 1 : debug output for every fluid step
    //   -> debug_level = 2 : debug output for every substep
    const int debug_level = 0;

    /****************************************************************************
     * Geometry                                                                 *
     ***************************************************************************/

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
    des_init_time_loop( &time, &dt, &nsubsteps, &subdt );

    /****************************************************************************
     * Init temporary storage:                                                  *
     *   -> particle-particle, and particle-wall forces                         *
     *   -> particle-particle, and particle-wall torques                        *
     ***************************************************************************/

    std::map<PairIndex, Gpu::ManagedDeviceVector<Real>> tow;
    std::map<PairIndex, Gpu::ManagedDeviceVector<Real>> fc, pfor, wfor;
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        tow[index]  = Gpu::ManagedDeviceVector<Real>();
        fc[index]   = Gpu::ManagedDeviceVector<Real>();
        pfor[index] = Gpu::ManagedDeviceVector<Real>();
        wfor[index] = Gpu::ManagedDeviceVector<Real>();
    }

    /****************************************************************************
     * Get particle EB geometric info
     ***************************************************************************/
    const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());

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
        // list (Note that this fills the neighbour list after every
        // redistribute operation)
        if (n % 25 == 0) {
            clearNeighbors();
            Redistribute();
            fillNeighbors();
            buildNeighborList(MFIXCheckPair(), sort_neighbor_list);
        } else {
            updateNeighbors();
        }

#ifdef _OPENMP
#pragma omp parallel reduction(+:ncoll) if (Gpu::notInLaunchRegion())
#endif
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            // Timer used for load-balancing
            Real wt = ParallelDescriptor::second();

            const Box& bx = pti.tilebox();
            PairIndex index(pti.index(), pti.LocalTileIndex());
            
            const int nrp = GetParticles(lev)[index].numRealParticles();            
            RealType* particles  = pti.GetArrayOfStructs().data();

            auto& plev = GetParticles(lev);
            auto& ptile = plev[index];
            auto& aos   = ptile.GetArrayOfStructs();
            ParticleType* pstruct = aos().dataPtr();

            // Neighbor particles
#ifdef AMREX_USE_CUDA
            int size_ng = aos.numNeighborParticles();            
#else
            int size_ng = neighbors[lev][index].size();
            int size_nl = neighbor_list[lev][index].size();
#endif

            // Number of particles including neighbor particles
            int ntot = nrp + size_ng;

            // Particle-particle (and particle-wall) forces and torques. We need
            // these to be zero every time we start a new batch (i.e tile and
            // substep) of particles.
            tow[index].clear();
            fc[index].clear();
            tow[index].resize(ntot*3,0.0);
            fc[index].resize(ntot*3,0.0);

            Real* fc_ptr = fc[index].dataPtr();
            Real* tow_ptr = tow[index].dataPtr();

            // For debugging: keep track of particle-particle (pfor) and
            // particle-wall (wfor) forces
            pfor[index].clear();
            wfor[index].clear();
            pfor[index].resize(3 * ntot, 0.0);
            wfor[index].resize(3 * ntot, 0.0);

            /********************************************************************
             * Particle-Wall collision forces (and torques)                     *
             *******************************************************************/

            // Only call the routine for wall collisions if we actually have walls
            bool has_wall = false;
            if ((ebfactory != NULL)
                && ((*flags)[pti].getType(amrex::grow(bx,1)) == FabType::singlevalued))
            {
                has_wall = true;
            }
            else
            {
                int int_has_wall = 0;
                Real tol = std::min(dx[0], std::min(dx[1], dx[2])) / 2;
                ls_has_walls(& int_has_wall, BL_TO_FORTRAN_3D((* ls_phi)[pti]), & tol);
                has_wall = (int_has_wall > 0);
            }

            if (has_wall)
            {
                // Calculate forces and torques from particle-wall collisions
                BL_PROFILE_VAR("calc_wall_collisions()", calc_wall_collisions);

                auto& geom = this->Geom(lev);
                const auto dxi = geom.InvCellSizeArray();
                const auto plo = geom.ProbLoArray();
                const auto phiarr = ls_phi->array(pti);

                AMREX_FOR_1D ( nrp, i,
                {
                    ParticleType& p = pstruct[i];
                    Real rp = p.rdata(realData::radius);
                    
                    Real ls_value = interp_level_set(p, ls_refinement, phiarr, plo, dxi);

                    Real overlap_n = rp - ls_value;

                    if (ls_value < rp)
                    {
                        Real normal[3];
                        level_set_normal(p, ls_refinement, &normal[0], phiarr, plo, dxi);

                        normal[0] *= -1;
                        normal[1] *= -1;
                        normal[2] *= -1;

                        Real v_rot[3];                       
                        v_rot[0] = ls_value * p.rdata(realData::omegax);
                        v_rot[1] = ls_value * p.rdata(realData::omegay);
                        v_rot[2] = ls_value * p.rdata(realData::omegaz);
                        
                        Real vreltrans[3];
                        Real cprod[3];
                        
                        cross_product(v_rot, normal, cprod);
                        vreltrans[0] = p.rdata(realData::velx) + cprod[0];
                        vreltrans[1] = p.rdata(realData::vely) + cprod[1];
                        vreltrans[2] = p.rdata(realData::velz) + cprod[2];

                        Real vreltrans_norm = dot_product(vreltrans, normal);
                        
                        Real vrel_t[3];
                        vrel_t[0] = vreltrans[0] - vreltrans_norm*normal[0];
                        vrel_t[1] = vreltrans[1] - vreltrans_norm*normal[1];
                        vrel_t[2] = vreltrans[2] - vreltrans_norm*normal[2];

                        int phase = p.idata(intData::phase);

                        Real kn_des_w;
                        Real kt_des_w; 
                        Real etan_des_w;
                        Real etat_des_w;

                        if (DEMParams::CollisionModel == DEMParams::HERTZIAN)
                        {                            
                            amrex::Abort("Not implemented");
                        }
                        else
                        {
                            kn_des_w   = DEMParams::kn_w;
                            kt_des_w   = DEMParams::kt_w;
                            etan_des_w = DEMParams::etan_w[phase-1];
                            etat_des_w = DEMParams::etat_w[phase-1];
                        }

                        Real fn[3];
                        Real ft[3];
                        Real overlap_t[3];
                        Real mag_overlap_t;

                        // calculate the normal contact force
                        fn[0] = -(kn_des_w*overlap_n*normal[0] + etan_des_w*vreltrans_norm*normal[0]);
                        fn[1] = -(kn_des_w*overlap_n*normal[1] + etan_des_w*vreltrans_norm*normal[1]);
                        fn[2] = -(kn_des_w*overlap_n*normal[2] + etan_des_w*vreltrans_norm*normal[2]);

                        // calculate the tangential displacement
                        overlap_t[0] = subdt*vrel_t[0];
                        overlap_t[1] = subdt*vrel_t[1];
                        overlap_t[2] = subdt*vrel_t[2];

                        mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));
                        
                        if (mag_overlap_t > 0.0) {                            
                            Real fnmd = DEMParams::mew * sqrt(dot_product(fn, fn));
                            Real tangent[3];
                            tangent[0] = overlap_t[0]/mag_overlap_t;
                            tangent[1] = overlap_t[1]/mag_overlap_t;
                            tangent[2] = overlap_t[2]/mag_overlap_t;
                            ft[0] = -fnmd * tangent[0];
                            ft[1] = -fnmd * tangent[1];
                            ft[2] = -fnmd * tangent[2];
                        } else {
                            ft[0] = 0.0;
                            ft[1] = 0.0;
                            ft[2] = 0.0;
                        }

                        // each particle updates its force (no need for atomics)
                        fc_ptr[i         ] += fn[0] + ft[0];
                        fc_ptr[i + ntot  ] += fn[1] + ft[1];
                        fc_ptr[i + 2*ntot] += fn[2] + ft[2];
                        
                        Real tow_force[3];

                        cross_product(normal, ft, tow_force);

                        tow_ptr[i         ] += ls_value*tow_force[0];
                        tow_ptr[i + ntot  ] += ls_value*tow_force[1];
                        tow_ptr[i + 2*ntot] += ls_value*tow_force[2];                        
                    }
                });                

                // Debugging: copy data from the fc (all forces) vector to
                // the wfor (wall forces) vector.
                if (debug_level > 0) {
                    for (int i = 0; i < wfor[index].size(); i++ ) {
                        wfor[index][i] = fc[index][i];
                    }
                }
                BL_PROFILE_VAR_STOP(calc_wall_collisions);
            }

            /********************************************************************
             * Particle-Particle collision forces (and torques)                 *
             *******************************************************************/

            BL_PROFILE_VAR("calc_particle_collisions()", calc_particle_collisions);

#ifdef AMREX_USE_CUDA
            auto nbor_data = m_neighbor_list[index].data();
            
            constexpr Real small_number = 1.0e-15;
            long ncoll = 0;
            long* pncoll = &ncoll;

            Real eps = std::numeric_limits<Real>::epsilon();
            
            // now we loop over the neighbor list and compute the forces
            AMREX_FOR_1D ( nrp, i,
            {
                ParticleType& p1 = pstruct[i];
                for (const auto& p2 : nbor_data.getNeighbors(i))
                {
                    Real dx = p2.pos(0) - p1.pos(0);
                    Real dy = p2.pos(1) - p1.pos(1);
                    Real dz = p2.pos(2) - p1.pos(2);
                    
                    Real r2 = dx*dx + dy*dy + dz*dz;
                    Real r_lm = p1.rdata(realData::radius) + p2.rdata(realData::radius);

                    if ( r2 <= (r_lm - small_number)*(r_lm - small_number) )
                    {
                        Cuda::Atomic::Add(pncoll, 1);
                        Real dist_mag = sqrt(r2);
                        AMREX_ASSERT(dist_mag >= eps);
                        
                        Real normal[3];
                        normal[0] = dx / dist_mag;
                        normal[1] = dy / dist_mag;
                        normal[2] = dz / dist_mag;
                        
                        Real overlap_n = r_lm - dist_mag;
                        Real vrel_trans_norm;
                        Real vrel_t[3];

                        cfrelvel(p1, p2, vrel_trans_norm, vrel_t, normal, dist_mag);

                        Real kn_des;
                        Real kt_des; 
                        Real etan_des;
                        Real etat_des;

                        int phase1 = p1.idata(intData::phase);
                        int phase2 = p2.idata(intData::phase);

                        if (DEMParams::CollisionModel == DEMParams::HERTZIAN)
                        {                            
                            amrex::Abort("Not implemented");
                        }
                        else
                            {
                            kn_des = DEMParams::kn;
                            kt_des = DEMParams::kt;
                            etan_des = DEMParams::etan[phase1-1][phase2-1];
                            etat_des = DEMParams::etat[phase1-1][phase2-1];
                        }

                        Real fn[3];
                        Real ft[3];
                        Real overlap_t[3];
                        Real mag_overlap_t;

                        // calculate the normal contact force
                        fn[0] = -(kn_des*overlap_n*normal[0] + etan_des*vrel_trans_norm*normal[0]);
                        fn[1] = -(kn_des*overlap_n*normal[1] + etan_des*vrel_trans_norm*normal[1]);
                        fn[2] = -(kn_des*overlap_n*normal[2] + etan_des*vrel_trans_norm*normal[2]);

                        // calculate the tangential overlap
                        overlap_t[0] = subdt*vrel_t[0];
                        overlap_t[1] = subdt*vrel_t[1];
                        overlap_t[2] = subdt*vrel_t[2];
                        mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));
                        
                        if (mag_overlap_t > 0.0) {                            
                            Real fnmd = DEMParams::mew * sqrt(dot_product(fn, fn));
                            Real tangent[3];
                            tangent[0] = overlap_t[0]/mag_overlap_t;
                            tangent[1] = overlap_t[1]/mag_overlap_t;
                            tangent[2] = overlap_t[2]/mag_overlap_t;
                            ft[0] = -fnmd * tangent[0];
                            ft[1] = -fnmd * tangent[1];
                            ft[2] = -fnmd * tangent[2];
                        } else {
                            ft[0] = 0.0;
                            ft[1] = 0.0;
                            ft[2] = 0.0;
                        }

                        // each particle updates its force (no need for atomics)
                        fc_ptr[i         ] += fn[0] + ft[0];
                        fc_ptr[i + ntot  ] += fn[1] + ft[1];
                        fc_ptr[i + 2*ntot] += fn[2] + ft[2];
                        
                        Real r1 = p1.rdata(realData::radius);
                        Real r2 = p2.rdata(realData::radius);
                                               
                        Real dist_cl = 0.5 * (dist_mag + (r1*r1 - r2*r2) / dist_mag);
                        dist_cl = dist_mag - dist_cl;
                        
                        Real tow_force[3];

                        cross_product(normal, ft, tow_force);

                        tow_ptr[i         ] += dist_cl*tow_force[0];
                        tow_ptr[i + ntot  ] += dist_cl*tow_force[1];
                        tow_ptr[i + 2*ntot] += dist_cl*tow_force[2];
                        } 
                    }
            });
#else            
            calc_particle_collisions ( particles                          , &nrp,
                                       neighbors[lev][index].dataPtr()    , &size_ng,
                                       neighbor_list[lev][index].dataPtr(), &size_nl,
                                       tow[index].dataPtr(), fc[index].dataPtr(),
                                       &subdt, &ncoll);
            
            // Debugging: copy data from the fc (all forces) vector to the wfor
            // (wall forces) vector. Note that since fc already contains the
            // wall forces, these need to be subtracted here.
            if (debug_level > 0)
            {
                for (int i = 0; i < pfor[index].size(); i++ ) {
                    pfor[index][i] = fc[index][i] - wfor[index][i];
                }
            }
#endif            
            BL_PROFILE_VAR_STOP(calc_particle_collisions);

            /********************************************************************
             * Move particles based on collision forces and torques             *
             *******************************************************************/

            Real grav[3];
            grav[0] = gravity[0];
            grav[1] = gravity[1];
            grav[2] = gravity[2];

            Real* grav_ptr = &grav[0];

            AMREX_FOR_1D ( nrp, i,
            {
                ParticleType& p = pstruct[i];
                
                p.rdata(realData::velx) += subdt * (
                    (p.rdata(realData::dragx) + fc_ptr[i       ]) /  p.rdata(realData::mass) + grav_ptr[0]);
                p.rdata(realData::vely) += subdt * (
                    (p.rdata(realData::dragy) + fc_ptr[i + ntot]) /  p.rdata(realData::mass) + grav_ptr[1]);
                p.rdata(realData::velz) += subdt * (
                    (p.rdata(realData::dragz) + fc_ptr[i+2*ntot]) /  p.rdata(realData::mass) + grav_ptr[2]);
                
                p.rdata(realData::omegax) += subdt * p.rdata(realData::oneOverI) * tow_ptr[i       ];
                p.rdata(realData::omegay) += subdt * p.rdata(realData::oneOverI) * tow_ptr[i+  ntot];
                p.rdata(realData::omegaz) += subdt * p.rdata(realData::oneOverI) * tow_ptr[i+2*ntot];
                
                p.pos(0) += subdt * p.rdata(realData::velx);
                p.pos(1) += subdt * p.rdata(realData::vely);
                p.pos(2) += subdt * p.rdata(realData::velz);
            });

            /********************************************************************
             * Update runtime cost (used in load-balancing)                     *
             *******************************************************************/
            
            if (cost)
            {
                // Runtime cost is either (weighted by tile box size):
                //   * time spent
                //   * number of particles
                const Box& tbx = pti.tilebox();
                if (knapsack_weight_type == "RunTimeCosts")
                {
                    wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
                }
                else if (knapsack_weight_type == "NumParticles")
                {
                    wt = nrp / tbx.d_numPts();
                }
                (*cost)[pti].plus(wt, tbx);
            }
        }

        // Update substep count
        n += 1;
        
        /************************************************************************
         * DEBUG: output the number of collisions in current substep            *
         *        output the max velocity (and forces) in current substep       *
         *        update max velocities and forces                              *
         ***********************************************************************/

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

    // Redistribute particles at the end of all substeps (note that the particle
    // neighbour list needs to be reset when redistributing).
    clearNeighbors();
    Redistribute();
    
    /****************************************************************************
     * DEBUG: output the total number of collisions over all substeps           *
     *        output the maximum velocity and forces over all substeps          *
     ***************************************************************************/
    if (debug_level > 0) {
        ParallelDescriptor::ReduceIntSum(ncoll_total, ParallelDescriptor::IOProcessorNumber());
        Print() << "Number of collisions: " << ncoll_total << " in " << nsubsteps << " substeps " << std::endl;
    }
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {        
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

void MFIXParticleContainer::CalcVolumeFraction(const Vector<std::unique_ptr<MultiFab>> & mf_to_be_filled,
                                               const Vector<std::unique_ptr<EBFArrayBoxFactory>> & ebfactory,
                                               const Vector<std::unique_ptr<IArrayBox>> & bc_ilo,
                                               const Vector<std::unique_ptr<IArrayBox>> & bc_ihi,
                                               const Vector<std::unique_ptr<IArrayBox>> & bc_jlo,
                                               const Vector<std::unique_ptr<IArrayBox>> & bc_jhi,
                                               const Vector<std::unique_ptr<IArrayBox>> & bc_klo,
                                               const Vector<std::unique_ptr<IArrayBox>> & bc_khi,
                                               int nghost )
{

    // NOTE: ebfactory HAS to be the PARTICLE EB factory!
    int fortran_volume_comp = 5;
    PICDeposition(mf_to_be_filled, ebfactory,
                  bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo,bc_khi,
                  fortran_volume_comp,nghost);

    for (int lev = 0; lev < nlev; lev++)
    {
        // Now define this mf = (1 - particle_vol)
        mf_to_be_filled[lev]->mult(-1.0,mf_to_be_filled[lev]->nGrow());
        mf_to_be_filled[lev]->plus( 1.0,mf_to_be_filled[lev]->nGrow());

        // We set ep_g to 1 rather than 0 in covered cells so that when we divide by ep_g
        //    following the projection we don't have to protect against divide by 0.
        EB_set_covered(*mf_to_be_filled[lev],1.0);

        // Impose a lower bound on volume fraction
        CapSolidsVolFrac(*mf_to_be_filled[lev]);
    }

    // HACK -- we really should average down (ep_g * volfrac) not ep_g.
    for (int lev = nlev - 1; lev > 0; lev --)
    {
        amrex::EB_average_down(* mf_to_be_filled[lev], * mf_to_be_filled[lev - 1],
                               0, 1, m_gdb->refRatio(lev - 1));
    }
}

void MFIXParticleContainer::CalcDragOnFluid(const Vector<std::unique_ptr<MultiFab>> & beta_mf,
                                            const Vector<std::unique_ptr<MultiFab>> & beta_vel_mf,
                                            const Vector<std::unique_ptr<EBFArrayBoxFactory>> & ebfactory,
                                            const Vector<std::unique_ptr<IArrayBox>> & bc_ilo,
                                            const Vector<std::unique_ptr<IArrayBox>> & bc_ihi,
                                            const Vector<std::unique_ptr<IArrayBox>> & bc_jlo,
                                            const Vector<std::unique_ptr<IArrayBox>> & bc_jhi,
                                            const Vector<std::unique_ptr<IArrayBox>> & bc_klo,
                                            const Vector<std::unique_ptr<IArrayBox>> & bc_khi,
                                            int nghost )
{
    int fortran_beta_comp = 15;
    int fortran_vel_comp  =  9;
    PICMultiDeposition(beta_mf, beta_vel_mf, ebfactory,
                       bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo, bc_khi,
                       fortran_beta_comp, fortran_vel_comp, nghost);
}

void MFIXParticleContainer::PICDeposition(const amrex::Vector< std::unique_ptr<MultiFab> >& mf_to_be_filled,
                                          const amrex::Vector< std::unique_ptr<EBFArrayBoxFactory>  >& ebfactory,
                                          const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
                                          const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
                                          const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
                                          const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
                                          const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
                                          const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
                                          int fortran_particle_comp, int nghost )
{
    BL_PROFILE("MFIXParticleContainer::PICDeposition()");

    int ncomp = 1;

    MultiFab* mf_pointer[nlev];

    // Start the timers ...
    const Real      strttime    = ParallelDescriptor::second();

    if (nlev > 2)
      amrex::Abort("For right now MFIXParticleContainer::PICDeposition can only handle up to 2 levels");

    for (int lev = 0; lev < nlev; lev++)
    {
       if (lev == 0 && OnSameGrids(lev, *mf_to_be_filled[lev])) {
          // If we are already working with the internal mf defined on the
          // particle_box_array, then we just work with this.
          mf_pointer[lev] = mf_to_be_filled[lev].get();

       } else if (lev == 0 && !OnSameGrids(lev, *mf_to_be_filled[lev]))  {
          // If mf_to_be_filled is not defined on the particle_box_array, then we need
          // to make a temporary here and copy into mf_to_be_filled at the end.
          mf_pointer[lev] = new MultiFab(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                                         ncomp, mf_to_be_filled[lev]->nGrow());

       } else {
          // If lev > 0 we make a temporary at the coarse resolution
          BoxArray ba_crse(amrex::coarsen(ParticleBoxArray(lev),this->m_gdb->refRatio(0)));
          mf_pointer[lev] = new MultiFab(ba_crse, ParticleDistributionMap(lev),
                                         ncomp, 1);
       }

       // We must have ghost cells for each FAB so that a particle in one grid can spread
       // its effect to an adjacent grid by first putting the value into ghost cells of its
       // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
       // to another grid's valid region.
       if (mf_pointer[lev]->nGrow() < 1)
          amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

       mf_pointer[lev]->setVal(0.0,0,1,mf_pointer[lev]->nGrow());
    }

    // We always use the coarse dx
    const Geometry& gm          = Geom(0);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx          = gm.CellSize();

    using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;

    const FabArray<EBCellFlagFab>* flags;

    for (int lev = 0; lev < nlev; lev++)
    {

       Vector<int> ngrow = {1,1,1};
       std::unique_ptr<EBFArrayBoxFactory> crse_factory;

       if (lev == 0) {
          flags   = &(ebfactory[lev]->getMultiEBCellFlagFab());
       } else {
          // std::unique_ptr<EBFArrayBoxFactory>
          crse_factory = makeEBFabFactory(
                 gm, mf_pointer[lev]->boxArray(), mf_pointer[lev]->DistributionMap(),
                 ngrow, EBSupport::volume);
          flags   = &(crse_factory->getMultiEBCellFlagFab());
       }

       const MultiFab* volfrac = (lev == 0) ? &(ebfactory[lev]->getVolFrac()) : &(crse_factory->getVolFrac());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       {
        FArrayBox local_vol;
        for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {
            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long nrp = pti.numParticles();
            FArrayBox& fab = (*mf_pointer[lev])[pti];
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

            const Box& bx  = pti.tilebox(); // I need a box without ghosts

            if ((*flags)[pti].getType(bx) != FabType::covered ) {
                mfix_deposit_cic_eb(particles.data(), nstride, nrp,
                                    data_ptr, lo, hi,
                                    BL_TO_FORTRAN_ANYD((*volfrac)[pti]),
                                    BL_TO_FORTRAN_ANYD((*flags)[pti]),
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

       for (MFIter mfi(*mf_pointer[lev]); mfi.isValid(); ++mfi) {

         const Box& sbx = (*mf_pointer[lev])[mfi].box();

         flip_particle_vol(sbx.loVect(), sbx.hiVect(),
                           (*mf_pointer[lev])[mfi].dataPtr(),
                           (*bc_ilo[lev]).dataPtr(), (*bc_ihi[lev]).dataPtr(),
                           (*bc_jlo[lev]).dataPtr(), (*bc_jhi[lev]).dataPtr(),
                           (*bc_klo[lev]).dataPtr(), (*bc_khi[lev]).dataPtr(),
                           domain.loVect(), domain.hiVect(),
                           &nghost );
       }
    }

    int  src_nghost = 1;
    int dest_nghost = 0;
    for (int lev = 1; lev < nlev; lev++)
        mf_pointer[0]->copy(*mf_pointer[lev],0,0,ncomp,src_nghost,dest_nghost,gm.periodicity(),FabArrayBase::ADD);

    mf_pointer[0]->SumBoundary(gm.periodicity());

    if (nlev > 1)
    {
        IntVect ref_ratio(this->m_gdb->refRatio(0));

        // Now interpolate from the coarse grid to define the fine grid ep-g
        Interpolater* mapper = &cell_cons_interp;
        int lo_bc[] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
        int hi_bc[] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
        Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

        BndryFuncArray bfunc(phifill);

        Real time = 0.0;
        for (int lev = 1; lev < nlev; lev++)
        {
            PhysBCFunct<BndryFuncArray> cphysbc(Geom(lev-1), bcs, bfunc);
            PhysBCFunct<BndryFuncArray> fphysbc(Geom(lev  ), bcs, bfunc);
            mf_to_be_filled[lev]->setVal(0.0);
            amrex::InterpFromCoarseLevel(*mf_to_be_filled[lev], time, *mf_pointer[lev-1],
                                         0, 0, 1, Geom(lev-1), Geom(lev),
                                         cphysbc, 0, fphysbc, 0,
                                         ref_ratio, mapper,
                                         bcs, 0);
        }
    }

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.

    if (mf_pointer[0] != mf_to_be_filled[0].get())
       mf_to_be_filled[0]->copy(*mf_pointer[0],0,0,ncomp);

    for (int lev = 0; lev < nlev; lev++)
       if (mf_pointer[lev] != mf_to_be_filled[lev].get())
          delete mf_pointer[lev];

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::PICDeposition time: " << stoptime << '\n';
    }
}

void MFIXParticleContainer::PICMultiDeposition(const amrex::Vector< std::unique_ptr<MultiFab> >& beta_mf,
                                               const amrex::Vector< std::unique_ptr<MultiFab> >& beta_vel_mf,
                                               const amrex::Vector< std::unique_ptr<EBFArrayBoxFactory>  >& ebfactory,
                                               const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
                                               const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
                                               const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
                                               const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
                                               const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
                                               const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
                                               int fortran_beta_comp, int fortran_vel_comp, int nghost )
{
    BL_PROFILE("MFIXParticleContainer::PICMultiDeposition()");

    const Real      strttime    = ParallelDescriptor::second();

    if (nlev > 2)
      amrex::Abort("For right now MFIXParticleContainer::PICMultiDeposition can only handle up to 2 levels");

    for (int lev = 0; lev < nlev; lev++)
       AMREX_ASSERT(OnSameGrids(lev,*beta_mf[lev])==OnSameGrids(lev,*beta_vel_mf[lev]));

    MultiFab*  beta_ptr[nlev];
    MultiFab*  beta_vel_ptr[nlev];

    for (int lev = 0; lev < nlev; lev++)
    {
       if (lev == 0 && OnSameGrids(lev, *beta_mf[lev])) {
          // If we are already working with the internal mf defined on the
          // particle_box_array, then we just work with this.
              beta_ptr[lev] = beta_mf[lev].get();
          beta_vel_ptr[lev] = beta_vel_mf[lev].get();

       } else if (lev == 0 && !OnSameGrids(lev, *beta_mf[lev]))  {
          // If beta_mf is not defined on the particle_box_array, then we need
          // to make a temporary here and copy into beta_mf at the end.
          beta_ptr[lev]     = new MultiFab(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                                           beta_mf[lev]->nComp(), beta_mf[lev]->nGrow());
          beta_vel_ptr[lev] = new MultiFab(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                                           beta_vel_mf[lev]->nComp(), beta_vel_mf[lev]->nGrow());

       } else {
          // If lev > 0 we make a temporary at the coarse resolution
          BoxArray ba_crse(amrex::coarsen(ParticleBoxArray(lev),this->m_gdb->refRatio(0)));
              beta_ptr[lev] = new MultiFab(ba_crse, ParticleDistributionMap(lev),beta_mf[lev]->nComp()    ,1);
          beta_vel_ptr[lev] = new MultiFab(ba_crse, ParticleDistributionMap(lev),beta_vel_mf[lev]->nComp(),1);
       }

       // We must have ghost cells for each FAB so that a particle in one grid can spread
       // its effect to an adjacent grid by first putting the value into ghost cells of its
       // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
       // to another grid's valid region.
       if (beta_ptr[lev]->nGrow() < 1)
          amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

           beta_ptr[lev]->setVal(0.0,0,1,    beta_ptr[lev]->nGrow());
       beta_vel_ptr[lev]->setVal(0.0,0,3,beta_vel_ptr[lev]->nGrow());
    }


    // We always use the coarse dx
    const Geometry& gm          = Geom(0);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx          = gm.CellSize();

    using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;
    const FabArray<EBCellFlagFab>* flags;

    for (int lev = 0; lev < nlev; lev++)
    {
        Vector<int> ngrow = {1,1,1};
        std::unique_ptr<EBFArrayBoxFactory> crse_factory;

        if (lev == 0) {
           flags   = &(ebfactory[lev]->getMultiEBCellFlagFab());
        } else {
           // std::unique_ptr<EBFArrayBoxFactory>
           crse_factory = makeEBFabFactory(
                  gm, beta_ptr[lev]->boxArray(), beta_ptr[lev]->DistributionMap(),
                  ngrow, EBSupport::volume);
           flags   = &(crse_factory->getMultiEBCellFlagFab());
        }

        const MultiFab* volfrac = (lev == 0) ? &(ebfactory[lev]->getVolFrac()) : &(crse_factory->getVolFrac());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {

        const int* lo;
        const int* hi;
        Real* bx_dataptr;
        Real* bu_dataptr;

        FArrayBox local_x_vol, local_u_vol;
         for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {

            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np = pti.numParticles();

            FArrayBox& beta_fab = (*beta_ptr[lev])[pti];
            FArrayBox& beta_vel_fab = (*beta_vel_ptr[lev])[pti];

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
            const Box& box = pti.tilebox(); // I need a box without ghosts

            if ((*flags)[pti].getType(box) != FabType::covered )
            {
                mfix_multi_deposit_cic_eb(particles.data(), nstride, np,
                                          bx_dataptr, bu_dataptr, lo, hi,
                                          BL_TO_FORTRAN_ANYD((*volfrac)[pti]),
                                          BL_TO_FORTRAN_ANYD((*flags)[pti]),
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
    }

    int  src_nghost = 1;
    int dest_nghost = 0;
    for (int lev = 1; lev < nlev; lev++)
    {
            beta_ptr[0]->copy(    *beta_ptr[lev],0,0,beta_ptr[0]->nComp()    ,src_nghost,dest_nghost,gm.periodicity(),FabArrayBase::ADD);
        beta_vel_ptr[0]->copy(*beta_vel_ptr[lev],0,0,beta_vel_ptr[0]->nComp(),src_nghost,dest_nghost,gm.periodicity(),FabArrayBase::ADD);
    }

    beta_ptr[0]->SumBoundary(gm.periodicity());
    beta_vel_ptr[0]->SumBoundary(gm.periodicity());

    if (nlev > 1)
    {
        IntVect ref_ratio(this->m_gdb->refRatio(0));

        // Now interpolate from the coarse grid to define the fine grid ep-g
        Interpolater* mapper = &cell_cons_interp;
        int lo_bc[] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
        int hi_bc[] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
        Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

        BndryFuncArray bfunc(phifill);

        Real time = 0.0;
        for (int lev = 1; lev < nlev; lev++)
        {
            PhysBCFunct<BndryFuncArray> cphysbc(Geom(lev-1), bcs, bfunc);
            PhysBCFunct<BndryFuncArray> fphysbc(Geom(lev  ), bcs, bfunc);
                beta_mf[lev]->setVal(0.0);
            beta_vel_mf[lev]->setVal(0.0);
            amrex::InterpFromCoarseLevel(*beta_mf[lev], time, *beta_ptr[lev-1],
                                         0, 0, 1, Geom(lev-1), Geom(lev),
                                         cphysbc, 0, fphysbc, 0,
                                         ref_ratio, mapper,
                                         bcs, 0);
            amrex::InterpFromCoarseLevel(*beta_vel_mf[lev], time, *beta_vel_ptr[lev-1],
                                         0, 0, 1, Geom(lev-1), Geom(lev),
                                         cphysbc, 0, fphysbc, 0,
                                         ref_ratio, mapper,
                                         bcs, 0);
        }
    }

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.

    if (beta_ptr[0] != beta_mf[0].get())
    {
           beta_mf[0]->copy(    *beta_ptr[0],0,0,beta_mf[0]->nComp());
       beta_vel_mf[0]->copy(*beta_vel_ptr[0],0,0,beta_vel_mf[0]->nComp());
    }

    for (int lev = 0; lev < nlev; lev++)
    {
       if (beta_ptr[lev] != beta_mf[lev].get())
          delete beta_ptr[lev];
       if (beta_vel_ptr[lev] != beta_vel_mf[lev].get())
          delete beta_vel_ptr[lev];
    }

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::PICMultiDeposition time: " << stoptime << '\n';
    }
}


void MFIXParticleContainer::writeAllAtLevel(int lev)
{
    // Not threaded because its print to terminal
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        Gpu::HostVector<ParticleType> host_particles;
        Cuda::thrust_copy(particles.begin(), particles.end(), host_particles.begin());

        for (const auto& p: host_particles)
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
#pragma omp parallel reduction(+:Np_tot) if (Gpu::notInLaunchRegion())
#endif
  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
      Np_tot += pti.numParticles();

  ParallelDescriptor::ReduceIntSum(Np_tot,ParallelDescriptor::IOProcessorNumber());

  cout << Np_tot << std::endl;

  // Not threaded because its print to terminal
  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
  {
      auto& particles = pti.GetArrayOfStructs();
      Gpu::HostVector<ParticleType> host_particles;
      Cuda::thrust_copy(particles.begin(), particles.end(), host_particles.begin());

      for (const auto& p: host_particles)
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

              Gpu::HostVector<ParticleType> host_particles;
              Cuda::thrust_copy(particles.begin(), particles.end(), host_particles.begin());

              int index = 0;
              for (const auto& p: host_particles)
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
#pragma omp parallel reduction(+:p_num, p_diam, p_dens) if (Gpu::notInLaunchRegion())
#endif
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& particles = pti.GetArrayOfStructs();
            
            Gpu::HostVector<ParticleType> host_particles(pti.numParticles());
            Cuda::thrust_copy(particles.begin(), particles.end(), host_particles.begin());
            
            for (const auto& p: host_particles){
                if ( phse==p.idata(intData::phase) )
                {
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
#pragma omp parallel reduction(max:max_vel_x,max_vel_y,max_vel_z) if (Gpu::notInLaunchRegion())
#endif
    for (int lev = 0; lev < nlev; lev++)
    {
       for(MFIXParIter pti(* this, lev); pti.isValid(); ++ pti)
       {
           auto & particles = pti.GetArrayOfStructs();
           Gpu::HostVector<ParticleType> host_particles;
           Cuda::thrust_copy(particles.begin(), particles.end(), host_particles.begin());

           for(const auto & particle : host_particles)
           {
              max_vel_x = std::max(Real(std::fabs(particle.rdata(realData::velx))), max_vel_x);
              max_vel_y = std::max(Real(std::fabs(particle.rdata(realData::vely))), max_vel_y);
              max_vel_z = std::max(Real(std::fabs(particle.rdata(realData::velz))), max_vel_z);
           }
       }
    }
    loc_maxvel = RealVect(max_vel_x, max_vel_y, max_vel_z);
}

void MFIXParticleContainer::UpdateMaxForces( std::map<PairIndex, Gpu::ManagedDeviceVector<Real>> pfor,
                                             std::map<PairIndex, Gpu::ManagedDeviceVector<Real>> wfor)
{
    Real max_pfor_x = loc_maxpfor[0], max_pfor_y = loc_maxpfor[1], max_pfor_z = loc_maxpfor[2];
    Real max_wfor_x = loc_maxwfor[0], max_wfor_y = loc_maxwfor[1], max_wfor_z = loc_maxwfor[2];

    for (int lev = 0; lev < nlev; lev++)
    {
#ifdef _OPENMP
#pragma omp parallel reduction(max:max_pfor_x,max_pfor_y,max_pfor_z,max_wfor_x,max_wfor_y,max_wfor_z) if (Gpu::notInLaunchRegion())
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

void MFIXParticleContainer::
ComputeAverageVelocities ( const int lev,
                           const amrex::Real time,
                           const string&  basename,
                           const Vector<int>& avg_vel_p,
                           const Gpu::ManagedDeviceVector<Real>& avg_region_x_w,
                           const Gpu::ManagedDeviceVector<Real>& avg_region_x_e,
                           const Gpu::ManagedDeviceVector<Real>& avg_region_y_s,
                           const Gpu::ManagedDeviceVector<Real>& avg_region_y_n,
                           const Gpu::ManagedDeviceVector<Real>& avg_region_z_b,
                           const Gpu::ManagedDeviceVector<Real>& avg_region_z_t )
{

  // Count number of calls -- Used to determin when to create file from scratch
  static int ncalls = 0;
  ++ncalls;

  int  nregions = avg_region_x_w.size();

  if(avg_vel_p.size() > 0)
    {

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

          amrex::Print() << "size of avg_vel_p " << avg_vel_p[nr] << "\n";

          // This region isn't needed for particle data.
          if( avg_vel_p[nr] == 0) continue;

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
#pragma omp parallel reduction(+:sum_np,sum_velx,sum_vely,sum_velz) if (Gpu::notInLaunchRegion())
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

              // Skip this region.
              if( avg_vel_p[nr] == 0 ) continue;

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

              fname = basename + "_vel_p_" + std::to_string(nr) + ".dat";

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

}

void MFIXParticleContainer::CapSolidsVolFrac(amrex::MultiFab& mf_to_be_filled)
{
    for (MFIter mfi(mf_to_be_filled); mfi.isValid(); ++mfi) {
       const Box& sbx = mf_to_be_filled[mfi].box();
       mfix_cap_eps(sbx.loVect(), sbx.hiVect(), (mf_to_be_filled)[mfi].dataPtr());
    }
}

void MFIXParticleContainer::set_particle_properties(int pstate, Real pradius, Real pdensity,
                                                    Real& pvol, Real& pmass, Real& omoi, Real& omega)
{
    pvol  = (4.0/3.0)*M_PI*(pradius*pradius*pradius);
    pmass = pvol * pdensity;
    omoi  = 2.5/(pmass * (pradius*pradius));
    omega = 0.0;
}
