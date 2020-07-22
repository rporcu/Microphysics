#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_RealVect.H>
#include <iostream>
#include <mfix_pc.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB_F.H>

#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBSupport.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

#include <cmath>

#include "mfix_des_K.H"
#include "mfix_dem_parms.H"
#include <mfix_pic_parms.H>
#include <mfix_ic_parms.H>

#include <mfix_particle_generator.H>

using namespace amrex;
using namespace std;

void MFIXParticleContainer::InitParticlesAscii (const std::string& file)
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

    auto& particle_tile = DefineAndReturnParticleTile(lev,grid,tile);
    //auto& particle_tile = GetParticles(lev)[std::make_pair(grid,tile)];

    ParticleType p;
    int  pstate, pphase;
    Real pradius, pdensity, pvolume, pomoi, pmass, pomega;

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
      set_particle_properties(pstate, pradius, pdensity, pvolume, pmass, pomoi, pomega);

      // Set other particle properties
      p.idata(intData::phase)       = pphase;
      p.idata(intData::state)       = pstate;
      p.rdata(realData::volume)     = pvolume;
      p.rdata(realData::density)    = pdensity;
      p.rdata(realData::mass)       = pmass;
      p.rdata(realData::oneOverI)   = pomoi;
      p.rdata(realData::radius)     = pradius;
      p.rdata(realData::omegax)     = pomega;
      p.rdata(realData::omegay)     = pomega;
      p.rdata(realData::omegaz)     = pomega;

      p.rdata(realData::statwt) = 1.0;

      // Initialize these for I/O purposes
      p.rdata(realData::dragcoeff) = 0.0;

      p.rdata(realData::dragx) = 0.0;
      p.rdata(realData::dragy) = 0.0;
      p.rdata(realData::dragz) = 0.0;

      p.rdata(realData::c_ps) = 0.0;
      p.rdata(realData::temperature) = 0.0;
      p.rdata(realData::convection) = 0.0;

      // Add everything to the data structure
      particle_tile.push_back(p);

      if (!ifs.good())
          amrex::Abort("Error initializing particles from Ascii file. \n");
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

      ParticleType p_new;
      for (int i = 0; i < pcount; i++) {
        // Set id and cpu for this particle
        p_new.id()  = ParticleType::NextID();
        p_new.cpu() = ParallelDescriptor::MyProc();

        // Add to the data structure
        particles.push_back(p_new);
        if (DEM::nspecies_dem > 0){
           for(int ii=0; ii < DEM::nspecies_dem; ++ii){
               particles.push_back_real(ii, -1.0);
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

  const Real * plo = Geom(lev).ProbLo();

  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
  {

    amrex::GpuArray<amrex::Real, DEM::NMAX> cp0_loc;
    for(int phase(0); phase<DEM::names.size(); phase++) {
      cp0_loc[phase] = DEM::c_p0[phase];
    }

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
      const amrex::RealVect ic_lo = {bx_lo[0]*dx[0], bx_lo[1]*dx[1], bx_lo[2]*dx[2]};
      const amrex::RealVect ic_hi = {bx_hi[0]*dx[0], bx_hi[1]*dx[1], bx_hi[2]*dx[2]};

      const Box ic_box(bx_lo, bx_hi);

      if ( pti.tilebox().intersects ( ic_box ) ){

        // Create a temporary copy of IC particle temperatures mapped
        // to the particle type.
        amrex::GpuArray<amrex::Real, DEM::NMAX> temperature_loc;
        for(int solid_type(0); solid_type<DEM::names.size(); solid_type++) {
          // Initialize to zero
          temperature_loc[solid_type] = 0.0;

          // Loop through IC solids looking for match.
          for(int ics(0); ics < IC::ic[icv].solids.size(); ics++) {
            DEM::DEM_t ic_solid = IC::ic[icv].solids[ics];
            if(DEM::names[solid_type] == ic_solid.name) {
              temperature_loc[solid_type] = ic_solid.temperature;
            }
          }
        }

        auto& particles = pti.GetArrayOfStructs();
        int np = pti.numParticles();

        auto particles_ptr = particles().dataPtr();

        amrex::ParallelFor(np,
          [particles_ptr, temperature_loc, cp0_loc, ic_lo, ic_hi]
          AMREX_GPU_DEVICE (int ip) noexcept
        {
          MFIXParticleContainer::ParticleType& p = particles_ptr[ip];

          if(ic_lo[0] <= p.pos(0) and p.pos(0) <= ic_hi[0] and
             ic_lo[1] <= p.pos(1) and p.pos(1) <= ic_hi[1] and
             ic_lo[2] <= p.pos(2) and p.pos(2) <= ic_hi[2])
          {
            const int phase = p.idata(intData::phase);
            p.rdata(realData::temperature) = temperature_loc[phase-1];
            p.rdata(realData::c_ps) = cp0_loc[phase-1];
          }
        });



      } // Intersecting Boxes
    } // IC regions
  } // MFIXParIter

}
