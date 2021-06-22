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
      host_realarrays[SoArealData::statwt][i]     = 1.0;

      // Initialize these for I/O purposes
      host_realarrays[SoArealData::dragcoeff][i]     = 0.0;
      host_realarrays[SoArealData::dragx][i]         = 0.0;
      host_realarrays[SoArealData::dragy][i]         = 0.0;
      host_realarrays[SoArealData::dragz][i]         = 0.0;
      host_realarrays[SoArealData::cp_s][i]          = 0.0;
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

    // Add components for each of the runtime variables
    const int start = AoSrealData::count + SoArealData::count;
    for (int comp(0); comp < m_runtimeRealData.count; ++comp)
      particles.push_back_real(start+comp, np, 0.);
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

      particles.resize(pcount);

      auto& aos = particles.GetArrayOfStructs();
      ParticleType* pstruct = aos().dataPtr();

      const int id = ParticleType::NextID();
      const int cpu = ParallelDescriptor::MyProc();

      amrex::ParallelFor(pcount, [pstruct,id,cpu]
        AMREX_GPU_DEVICE (int i) noexcept
      {
        ParticleType& part = pstruct[i];

        part.id() = id+i;
        part.cpu() = cpu;
      });

      // Update the particles NextID
      ParticleType::NextID(id+pcount);

      Gpu::synchronize();

      // Add components for each of the runtime variables
      const int start = AoSrealData::count + SoArealData::count;
      for (int comp(0); comp < m_runtimeRealData.count; ++comp)
        particles.push_back_real(start+comp, pcount, 0.);

      const int np = pcount;
      total_np += np;

      // Now define the rest of the particle data and store it directly in the particles
      if (pcount > 0)
        particles_generator.generate_prop(np, particles);
  }

  ParallelDescriptor::ReduceIntSum(total_np);
  amrex::Print() << "Total number of generated particles: " << total_np << std::endl;
  m_total_numparticle = total_np;

  // We shouldn't need this if the particles are tiled with one tile per grid, but otherwise
  // we do need this to move particles from tile 0 to the correct tile.
  Redistribute();

}

void MFIXParticleContainer::InitParticlesRuntimeVariables (const int adv_enthalpy,
                                                           const int solve_species)
{
  int lev = 0;
  const auto dx  = Geom(lev).CellSizeArray();
  const auto idx = Geom(lev).InvCellSizeArray();

  const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();

  const int nspecies_s = solids.nspecies;
  const int solid_is_a_mixture = solids.is_a_mixture;

  auto& solids_parms = *solids.parameters;

  // Create a temporary copy of IC particle temperatures mapped
  // to the particle type.
  Gpu::HostVector<Real> h_temperature_loc(SolidsPhase::NMAX);
  Gpu::HostVector<Real> h_mass_fractions(nspecies_s);

  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
  {
    auto& particles = pti.GetArrayOfStructs();
    int np = pti.numParticles();

    auto& soa = pti.GetStructOfArrays();
    auto p_realarray = soa.realarray();
    auto p_intarray = soa.intarray();

    auto particles_ptr = particles().dataPtr();

    PairIndex index(pti.index(), pti.LocalTileIndex());
    auto& plev = GetParticles(lev);
    auto& ptile = plev[index];
    auto ptile_data = ptile.getParticleTileData();

    // Set the initial conditions.
    for(int icv(0); icv < IC::ic.size(); ++icv)
    {
      IC::IC_t& loc_ic = IC::ic[icv];

      // This is a round about way to address what is likely overly complex
      // logic in the particle generator. We take the region specified by
      // a user and convert it index space (i,j,k). That in turn is turned
      // back into a physical region that may be a little larger than what
      // was actually defined to account for spatial discretization.
      const IntVect bx_lo(static_cast<int>(amrex::Math::floor((loc_ic.region->lo(0)-plo[0])*idx[0] + 0.5)),
                          static_cast<int>(amrex::Math::floor((loc_ic.region->lo(1)-plo[1])*idx[1] + 0.5)),
                          static_cast<int>(amrex::Math::floor((loc_ic.region->lo(2)-plo[2])*idx[0] + 0.5)));

      const IntVect bx_hi(static_cast<int>(amrex::Math::floor((loc_ic.region->hi(0)-plo[0])*idx[0] + 0.5)),
                          static_cast<int>(amrex::Math::floor((loc_ic.region->hi(1)-plo[1])*idx[1] + 0.5)),
                          static_cast<int>(amrex::Math::floor((loc_ic.region->hi(2)-plo[2])*idx[0] + 0.5)));

      // Start/end of IC domain bounds
      const RealVect ic_lo = {plo[0]+bx_lo[0]*dx[0],
                              plo[1]+bx_lo[1]*dx[1],
                              plo[2]+bx_lo[2]*dx[2]};

      const RealVect ic_hi = {plo[0]+bx_hi[0]*dx[0],
                              plo[1]+bx_hi[1]*dx[1],
                              plo[2]+bx_hi[2]*dx[2]};

      const Box ic_box(bx_lo, bx_hi);

      if (pti.tilebox().intersects(ic_box)) {

        for (int solid_type(0); solid_type < solids.names.size(); solid_type++) {
          // Initialize to zero
          h_temperature_loc[solid_type] = 0.0;

          // Loop through IC solids looking for match.
          for (int ics(0); ics < loc_ic.solids.size(); ics++) {

            auto& ic_solid = loc_ic.solids[ics];

            if (solids.names[solid_type] == ic_solid.name) {

              if (adv_enthalpy)
                h_temperature_loc[solid_type] = ic_solid.temperature;

              if (solve_species)
                for (int n_s(0); n_s < nspecies_s; n_s++) {
                  h_mass_fractions[n_s] = ic_solid.species[n_s].mass_fraction;
              }
            }
          }
        }

        Gpu::AsyncArray<Real> d_temperature_loc(h_temperature_loc.data(), h_temperature_loc.size());
        Real* p_temperature_loc = adv_enthalpy ? d_temperature_loc.data() : nullptr;

        Gpu::AsyncArray<Real> d_mass_fractions(h_mass_fractions.data(), h_mass_fractions.size());
        Real* p_mass_fractions = solve_species ? d_mass_fractions.data() : nullptr;

        amrex::ParallelFor(np,
          [particles_ptr,p_realarray,p_intarray,ptile_data,p_temperature_loc,
           p_mass_fractions,ic_lo,ic_hi,nspecies_s,solid_is_a_mixture,adv_enthalpy,
           solids_parms,solve_species]
          AMREX_GPU_DEVICE (int ip) noexcept
        {
          MFIXParticleContainer::ParticleType& p = particles_ptr[ip];

          if(ic_lo[0] <= p.pos(0) && p.pos(0) <= ic_hi[0] &&
             ic_lo[1] <= p.pos(1) && p.pos(1) <= ic_hi[1] &&
             ic_lo[2] <= p.pos(2) && p.pos(2) <= ic_hi[2])
          {
            const int phase = p_intarray[SoAintData::phase][ip];

            if(adv_enthalpy) {
              p_realarray[SoArealData::temperature][ip] = p_temperature_loc[phase-1];
            }

            if(adv_enthalpy) {
              if(!solid_is_a_mixture) {
                p_realarray[SoArealData::cp_s][ip] = solids_parms.calc_cp_s<RunOn::Gpu>(phase-1,
                                                                            p_temperature_loc[phase-1]);
              }
              else {
                Real cp_s_sum(0);
                for (int n_s(0); n_s < nspecies_s; n_s++)
                  cp_s_sum += p_mass_fractions[n_s]* solids_parms.calc_cp_sn<RunOn::Gpu>(p_temperature_loc[phase-1],n_s);

                p_realarray[SoArealData::cp_s][ip] = cp_s_sum;
              }
            }

            if(solve_species)
              for (int n_s(0); n_s < nspecies_s; n_s++)
                ptile_data.m_runtime_rdata[n_s][ip] = p_mass_fractions[n_s];
          }
        });

      } // Intersecting Boxes
    } // IC regions
  } // MFIXParIter

}
