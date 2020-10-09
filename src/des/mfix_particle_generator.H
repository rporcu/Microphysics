#ifndef _PARTICLES_GENERATOR_HPP_
#define _PARTICLES_GENERATOR_HPP_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_GpuContainers.H>
#include <mfix_pc.H>

class ParticlesGenerator {

  //< Position............... 1,2,3
  //< Radius................. 4
  //< Density................ 5
  //< Linear velocity........ 6,7,8
  //< Stat weight ........... 9
  static constexpr int nr = 9;

  //< Type................... 1
  static constexpr int ni = 1;

  public:
    using ParticleType = MFIXParticleContainer::ParticleType;
    using ParticleTileType = MFIXParticleContainer::ParticleTileType;

    // Constructor
    explicit ParticlesGenerator ();

    // Destructor
    ~ParticlesGenerator ();

    // Particle generator
    void generate (int& pc,
                   const amrex::IntVect& lo,
                   const amrex::IntVect& hi,
                   const amrex::Real dx,
                   const amrex::Real dy,
                   const amrex::Real dz,
                   const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo);

    void generate_prop (const int nrp, ParticleTileType& particles);

    // Write particles data
    void write (const int nrp,
                const ParticleTileType& particles,
                const int nstep = 0) const;

    void hex_close_pack (const int icv,
                         const int type,
                         const amrex::IntVect& lo,
                         const amrex::IntVect& hi,
                         int& np,
                         int& pc,
                         const amrex::Real dx,
                         const amrex::Real dy,
                         const amrex::Real dz,
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo);

    void one_per_fill (const int icv,
                       const int type,
                       const amrex::IntVect& lo,
                       const amrex::IntVect& hi,
                       int& np,
                       int& pc,
                       const amrex::Real dx,
                       const amrex::Real dy,
                       const amrex::Real dz,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo);

    void eight_per_fill (const int icv,
                         const int type,
                         const amrex::IntVect& lo,
                         const amrex::IntVect& hi,
                         int& np,
                         int& pc,
                         const amrex::Real dx,
                         const amrex::Real dy,
                         const amrex::Real dz,
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo);

    void random_fill_dem (const int icv,
                          const int type,
                          const amrex::IntVect& lo,
                          const amrex::IntVect& hi,
                          int& np,
                          int& pc,
                          const amrex::Real dx,
                          const amrex::Real dy,
                          const amrex::Real dz,
                          const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
                          const bool fix_seed);

    void random_fill_pic (const int icv,
                          const int type,
                          const amrex::IntVect& lo,
                          const amrex::IntVect& hi,
                          int& np,
                          int& pc,
                          const amrex::Real dx,
                          const amrex::Real dy,
                          const amrex::Real dz,
                          const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
                          const bool fix_seed);

    void nor_rno (amrex::Gpu::DeviceVector<amrex::Real>& dp,
                  const amrex::Real mean,
                  const amrex::Real sigma,
                  const amrex::Real dp_min,
                  const amrex::Real dp_max);

    void uni_rno (amrex::Gpu::DeviceVector<amrex::Real>& dp,
                  const amrex::Real dp_min,
                  const amrex::Real dp_max);

    void grow_pdata (const int gsize);

  private:
    amrex::Gpu::DeviceVector<amrex::Real> m_rdata;
    amrex::Gpu::DeviceVector<int> m_idata;
};

#endif
