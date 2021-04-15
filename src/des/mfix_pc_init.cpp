#include <mfix_solids_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_species_parms.H>
#include <mfix_ic_parms.H>

#include <mfix_particle_generator.H>

using namespace amrex;

void MFIXParticleContainer::InitParticlesAscii (const std::string& file)
{
  // only read the file on the IO proc
  if (ParallelDescriptor::IOProcessor())
  {
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

    auto& particles = DefineAndReturnParticleTile(lev,grid,tile);
    particles.resize(np);

    Gpu::HostVector<ParticleType> host_particles(np);

    std::array<Gpu::HostVector<Real>, SoArealData::count> host_realarrays;
    std::array<Gpu::HostVector<int>, SoAintData::count> host_intarrays;

    for (int comp(0); comp < SoArealData::count; ++comp)
      host_realarrays[comp].resize(np);

    for (int comp(0); comp < SoAintData::count; ++comp)
      host_intarrays[comp].resize(np);

    int  pstate, pphase;
    Real velx, vely, velz;
    Real pradius, pdensity, pvolume, pomoi, pmass, pomega;

    pstate = 1;

    int max_particle_phase(-1);

    for (int i = 0; i < np; i++)
    {
      // Read from input file
      ifs >> pphase;

      max_particle_phase = amrex::max(max_particle_phase, pphase);

      ifs >> host_particles[i].pos(0);
      ifs >> host_particles[i].pos(1);
      ifs >> host_particles[i].pos(2);
      ifs >> pradius;
      ifs >> pdensity;
      ifs >> velx;
      ifs >> vely;
      ifs >> velz;

      host_realarrays[SoArealData::velx][i]   = velx;
      host_realarrays[SoArealData::vely][i]   = vely;
      host_realarrays[SoArealData::velz][i]   = velz;

      // Compute other particle properties
      set_particle_properties(pstate, pradius, pdensity, pvolume, pmass, pomoi, pomega);

      // Set id and cpu for this particle
      host_particles[i].id()  = ParticleType::NextID();
      host_particles[i].cpu() = ParallelDescriptor::MyProc();

      // Set other particle properties
      host_intarrays[SoAintData::phase][i]        = pphase;
      host_intarrays[SoAintData::state][i]        = pstate;
      host_realarrays[SoArealData::volume][i]     = pvolume;
      host_realarrays[SoArealData::density][i]    = pdensity;
      host_realarrays[SoArealData::mass][i]       = pmass;
      host_realarrays[SoArealData::oneOverI][i]   = pomoi;
      host_realarrays[SoArealData::radius][i]     = pradius;
      host_realarrays[SoArealData::omegax][i]     = pomega;
      host_realarrays[SoArealData::omegay][i]     = pomega;
      host_realarrays[SoArealData::omegaz][i]     = pomega;
      host_realarrays[SoArealData::statwt][i]   = 1.0;

      // Initialize these for I/O purposes
      host_realarrays[SoArealData::dragcoeff][i]     = 0.0;
      host_realarrays[SoArealData::dragx][i]         = 0.0;
      host_realarrays[SoArealData::dragy][i]         = 0.0;
      host_realarrays[SoArealData::dragz][i]         = 0.0;
      host_realarrays[SoArealData::c_ps][i]          = 0.0;
      host_realarrays[SoArealData::temperature][i]   = 0.0;
      host_realarrays[SoArealData::convection][i]    = 0.0;

      if (!ifs.good())
          amrex::Abort("Error initializing particles from Ascii file. \n");
    }

    // NOTE : No need to do a ParallelDescriptor::ReduceIntMax on
    // max_particle_phase because we're reading the particle_input.dat only on
    // the IO proc

    if (max_particle_phase > DEM::NPHASE)
      amrex::Abort("One or more particle in the particle_input.dat has a phase number that is not present in the inputs file");

    auto& aos = particles.GetArrayOfStructs();
    Gpu::DeviceVector<ParticleType>& gpu_particles = aos();

    // Copy particles from host to device
    Gpu::copyAsync(Gpu::hostToDevice, host_particles.begin(), host_particles.end(), gpu_particles.begin());

    auto& soa = particles.GetStructOfArrays();
    auto p_realarray = soa.realarray();
    auto p_intarray = soa.intarray();

    // Copy particles from host to device
    for (int comp(0); comp < SoArealData::count; ++comp) {
      Gpu::copyAsync(Gpu::hostToDevice, host_realarrays[comp].begin(),
          host_realarrays[comp].end(), &(p_realarray[comp][0]));
    }

    // Copy particles from host to device
    for (int comp(0); comp < SoAintData::count; ++comp) {
      Gpu::copyAsync(Gpu::hostToDevice, host_intarrays[comp].begin(),
          host_intarrays[comp].end(), &(p_intarray[comp][0]));
    }

    // Add real components for solid species
    if (SOLIDS::solve_species)
    {
      // Add SOLIDS::nspecies components for each of the new species vars
      for (int n_s(0); n_s < SOLIDS::nspecies; ++n_s)
        particles.push_back_real(n_s, np, 0.);
    }

    // Add real components for chemical reactions rates
    if (SOLIDS::solve_species && REACTIONS::solve)
    {
      // Add components after solid species
      const int gap = SOLIDS::nspecies;

      // Add SOLIDS::nspecies components for each of the reactions
      for (int n_s(0); n_s < SOLIDS::nspecies; ++n_s) {
        for(int q(0); q < REACTIONS::nreactions; ++q) {
          const int comp = gap + n_s*REACTIONS::nreactions + q;

          particles.push_back_real(comp, np, 0.);
        }
      }
    }

  }

  Redistribute();
}

void MFIXParticleContainer::InitParticlesAuto ()
{
  int lev = 0;

  Real dx = Geom(lev).CellSize(0);
  Real dy = Geom(lev).CellSize(1);
  Real dz = Geom(lev).CellSize(2);

  int total_np = 0;

  ParticlesGenerator particles_generator;

  // This uses the particle tile size. Note that the default is to tile so if we
  //      remove the true and don't explicitly add false it will still tile
  for (MFIter mfi = MakeMFIter(lev,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      // This is particles per grid so we reset to 0
      int pcount = 0;

      // Define the real particle data for one grid-at-a-time's worth of particles
      // We don't know pcount (number of particles per grid) before this call

      const Box& tilebx = mfi.tilebox();

      const IntVect lo(tilebx.loVect());
      const IntVect hi(tilebx.hiVect());

      const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();

      particles_generator.generate(pcount, lo, hi, dx, dy, dz, plo);

      // Now that we know pcount, go ahead and create a particle container for this
      // grid and add the particles to it
      auto& particles = DefineAndReturnParticleTile(lev,mfi);

      SuperParticleType p_new;

      // Set p_new.pos() to avoid warnings
      p_new.m_pos[0] = 0;
      p_new.m_pos[1] = 0;
      p_new.m_pos[2] = 0;
      // Set p_new.id() to avoid warnings
      p_new.m_idata[0] = 0;
      // Set p_new.cpu() to avoid warnings
      p_new.m_idata[1] = 0;

      // If possible Parallelize this
      for (int i = 0; i < pcount; i++) {
        // Set id and cpu for this particle
        p_new.id()  = ParticleType::NextID();
        p_new.cpu() = ParallelDescriptor::MyProc();

        // Add to the data structure
        particles.push_back(p_new);

        // Add real components for solid species
        if (SOLIDS::solve_species)
        {
          // Add SOLIDS::nspecies components for each of the new species vars
          for (int n_s(0); n_s < SOLIDS::nspecies; ++n_s)
            particles.push_back_real(n_s, 0.);
        }

        // Add real components for solid species
        if (SOLIDS::solve_species && REACTIONS::solve)
        {
          const int gap = SOLIDS::nspecies;

          // Add SOLIDS::nspecies components for each of the reactions
          for (int n_s(0); n_s < SOLIDS::nspecies; ++n_s) {
            for(int q(0); q < REACTIONS::nreactions; ++q) {
              const int comp = gap + n_s*REACTIONS::nreactions + q;
              particles.push_back_real(comp, 0.);
            }
          }
        }
      }

      const int np = pcount;
      total_np += np;

      // Now define the rest of the particle data and store it directly in the particles
      if (pcount > 0)
        particles_generator.generate_prop(np, particles);
  }

  ParallelDescriptor::ReduceIntSum(total_np,ParallelDescriptor::IOProcessorNumber());
  amrex::Print() << "Total number of generated particles: " << total_np << std::endl;

  // We shouldn't need this if the particles are tiled with one tile per grid, but otherwise
  // we do need this to move particles from tile 0 to the correct tile.
  Redistribute();

}

void MFIXParticleContainer::InitParticlesEnthalpy ()
{
  int lev = 0;
  const auto dx  = Geom(lev).CellSizeArray();
  const auto idx = Geom(lev).InvCellSizeArray();

  const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();

  // Create a temporary copy of IC particle temperatures mapped
  // to the particle type.
  Gpu::HostVector<Real> h_temperature_loc(SOLIDS::NMAX);
  Gpu::HostVector<Real> h_cp0_loc(SOLIDS::NMAX);
  
  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
  {
    for(int phase(0); phase<SOLIDS::names.size(); phase++) {
      h_cp0_loc[phase] = SOLIDS::cp_p0[phase];
    }

    Gpu::AsyncArray<Real> d_dcp0_loc(h_cp0_loc.data(), h_cp0_loc.size());
    Real* p_cp0_loc = d_dcp0_loc.data();

    // Set the initial conditions.
    for(int icv(0); icv < IC::ic.size(); ++icv)
    {
      IC::IC_t ic = IC::ic[icv];

      // This is a round about way to address what is likely overly complex
      // logic in the particle generator. We take the region specified by
      // a user and convert it index space (i,j,k). That in turn is turned
      // back into a physical region that may be a little larger than what
      // was actually defined to account for spatial discretization.
      const IntVect bx_lo(amrex::Math::floor((ic.region->lo(0)-plo[0])*idx[0] + 0.5),
                          amrex::Math::floor((ic.region->lo(1)-plo[1])*idx[1] + 0.5),
                          amrex::Math::floor((ic.region->lo(2)-plo[2])*idx[0] + 0.5));

      const IntVect bx_hi(amrex::Math::floor((ic.region->hi(0)-plo[0])*idx[0] + 0.5),
                          amrex::Math::floor((ic.region->hi(1)-plo[1])*idx[1] + 0.5),
                          amrex::Math::floor((ic.region->hi(2)-plo[2])*idx[0] + 0.5));

      // Start/end of IC domain bounds
      const RealVect ic_lo = {plo[0]+bx_lo[0]*dx[0],
                                     plo[1]+bx_lo[1]*dx[1],
                                     plo[2]+bx_lo[2]*dx[2]};

      const RealVect ic_hi = {plo[0]+bx_hi[0]*dx[0],
                                     plo[1]+bx_hi[1]*dx[1],
                                     plo[2]+bx_hi[2]*dx[2]};

      const Box ic_box(bx_lo, bx_hi);

      if ( pti.tilebox().intersects ( ic_box ) )
      {
        for(int solid_type(0); solid_type<SOLIDS::names.size(); solid_type++) {
          // Initialize to zero
          h_temperature_loc[solid_type] = 0.0;

          // Loop through IC solids looking for match.
          for(int ics(0); ics < IC::ic[icv].solids.size(); ics++) {
            SOLIDS::SOLIDS_t ic_solid = IC::ic[icv].solids[ics];
            if(SOLIDS::names[solid_type] == ic_solid.name) {
              h_temperature_loc[solid_type] = ic_solid.temperature;
            }
          }
        }

        auto& particles = pti.GetArrayOfStructs();
        int np = pti.numParticles();

        Gpu::AsyncArray<Real> d_temperature_loc(h_temperature_loc.data(), h_temperature_loc.size());
        Real* p_temperature_loc = d_temperature_loc.data();

        auto& soa = pti.GetStructOfArrays();
        auto p_realarray = soa.realarray();
        auto p_intarray = soa.intarray();

        auto particles_ptr = particles().dataPtr();

        amrex::ParallelFor(np,
          [particles_ptr,p_realarray,p_intarray,p_temperature_loc,p_cp0_loc,ic_lo,ic_hi]
          AMREX_GPU_DEVICE (int ip) noexcept
        {
          MFIXParticleContainer::ParticleType& p = particles_ptr[ip];

          if(ic_lo[0] <= p.pos(0) && p.pos(0) <= ic_hi[0] &&
             ic_lo[1] <= p.pos(1) && p.pos(1) <= ic_hi[1] &&
             ic_lo[2] <= p.pos(2) && p.pos(2) <= ic_hi[2])
          {
            const int phase = p_intarray[SoAintData::phase][ip];
            p_realarray[SoArealData::temperature][ip] = p_temperature_loc[phase-1];
            p_realarray[SoArealData::c_ps][ip] = p_cp0_loc[phase-1];
          }
        });

      } // Intersecting Boxes
    } // IC regions
  } // MFIXParIter

}


void MFIXParticleContainer::InitParticlesSpecies ()
{
  int lev = 0;
  const auto dx  = Geom(lev).CellSizeArray();
  const auto idx = Geom(lev).InvCellSizeArray();

  const Real * plo = Geom(lev).ProbLo();

  const int nspecies_s = SOLIDS::nspecies;

  // Create a temporary copy of IC particle mass fractions mapped
  // to the particle type.
  Gpu::HostVector<Real> h_mass_fractions(nspecies_s);

  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
  {
    auto& particles = pti.GetArrayOfStructs();
    int np = pti.numParticles();

    auto particles_ptr = particles().dataPtr();

    PairIndex index(pti.index(), pti.LocalTileIndex());
    auto& plev = GetParticles(lev);
    auto& ptile = plev[index];
    auto ptile_data = ptile.getParticleTileData();

    // Set the initial conditions.
    for(int icv(0); icv < IC::ic.size(); ++icv)
    {
      IC::IC_t ic = IC::ic[icv];

      // This is a round about way to address what is likely overly complex
      // logic in the particle generator. We take the region specified by
      // a user and convert it index space (i,j,k). That in turn is turned
      // back into a physical region that may be a little larger than what
      // was actually defined to account for spatial discretization.
      const IntVect bx_lo(amrex::Math::floor((ic.region->lo(0)-plo[0])*idx[0] + 0.5),
                          amrex::Math::floor((ic.region->lo(1)-plo[1])*idx[1] + 0.5),
                          amrex::Math::floor((ic.region->lo(2)-plo[2])*idx[0] + 0.5));

      const IntVect bx_hi(amrex::Math::floor((ic.region->hi(0)-plo[0])*idx[0] + 0.5),
                          amrex::Math::floor((ic.region->hi(1)-plo[1])*idx[1] + 0.5),
                          amrex::Math::floor((ic.region->hi(2)-plo[2])*idx[0] + 0.5));

      // Start/end of IC domain bounds
      const RealVect ic_lo = {plo[0]+bx_lo[0]*dx[0],
                                     plo[1]+bx_lo[1]*dx[1],
                                     plo[2]+bx_lo[2]*dx[2]};

      const RealVect ic_hi = {plo[0]+bx_hi[0]*dx[0],
                                     plo[1]+bx_hi[1]*dx[1],
                                     plo[2]+bx_hi[2]*dx[2]};

      const Box ic_box(bx_lo, bx_hi);

      if (pti.tilebox().intersects(ic_box))
      {
        for(int solid_type(0); solid_type<SOLIDS::names.size(); solid_type++)
        {
          // Loop through IC solids looking for match.
          for(int ics(0); ics < IC::ic[icv].solids.size(); ics++)
          {
            const SOLIDS::SOLIDS_t& ic_solid = IC::ic[icv].solids[ics];

            if(SOLIDS::names[solid_type] == ic_solid.name) {
              for (int n_s(0); n_s < nspecies_s; n_s++) {
                h_mass_fractions[n_s] = ic_solid.species[n_s].mass_fraction;
              }
            }
          }

          Gpu::AsyncArray<Real> d_mass_fractions(h_mass_fractions.data(), h_mass_fractions.size());
          Real* p_mass_fractions = d_mass_fractions.data();

          amrex::ParallelFor(np, [particles_ptr,ptile_data,p_mass_fractions,
              ic_lo,ic_hi,nspecies_s]
            AMREX_GPU_DEVICE (int ip) noexcept
          {
            MFIXParticleContainer::ParticleType& p = particles_ptr[ip];

            if(ic_lo[0] <= p.pos(0) && p.pos(0) <= ic_hi[0] &&
               ic_lo[1] <= p.pos(1) && p.pos(1) <= ic_hi[1] &&
               ic_lo[2] <= p.pos(2) && p.pos(2) <= ic_hi[2])
            {
              for (int n_s(0); n_s < nspecies_s; n_s++)
                ptile_data.m_runtime_rdata[n_s][ip] = p_mass_fractions[n_s];
            }
          });
        }

      } // Intersecting Boxes
    } // IC regions
  } // MFIXParIter

}
