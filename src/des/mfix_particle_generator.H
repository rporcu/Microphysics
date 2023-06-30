#ifndef MFIX_PARTICLES_GENERATOR_H_
#define MFIX_PARTICLES_GENERATOR_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_GpuContainers.H>

#include <mfix_pc.H>
#include <mfix_ic.H>
#include <mfix_dem.H>


class ParticlesGenerator {

  public:
    using ParticleType = MFIXParticleContainer::ParticleType;
    using ParticleTileType = MFIXParticleContainer::ParticleTileType;

    // Constructor
    ParticlesGenerator (const amrex::IntVect& lo,
                        const amrex::IntVect& hi,
                        const amrex::RealVect& plo,
                        const amrex::RealVect& dx,
                        const int id,
                        const int cpu,
                        const int icv,
                        const int phase,
                        MFIXInitialConditions& initial_conditions,
                        MFIXDEM& dem,
                        MFIXPIC& pic);

    // Particle generator
    void generate (int& particles_count,
                   ParticleTileType& particles,
                   const amrex::RealBox* regions,
                   const int regions_nb,
                   const int allow_overlap);

    // Particle generator
    template <typename F1>
    void generate (int& particles_count,
                   ParticleTileType& particles,
                   const amrex::RealBox* regions,
                   const int regions_nb,
                   const int allow_overlap,
                   F1 positions_generator);

    // Write particles data
    void write (const int nrp,
                ParticleTileType& particles,
                const int nstep = 0) const;

    amrex::Gpu::HostVector<amrex::Real> m_h_data;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_data;

  private:
    const amrex::IntVect& m_bx_lo;
    const amrex::IntVect& m_bx_hi;
    const amrex::RealVect& m_plo;
    const amrex::RealVect& m_dx;
    const int m_id;
    const int m_cpu;
    const int m_icv;
    const int m_phase;

    MFIXInitialConditions& m_initial_conditions;
    MFIXDEM& m_dem;
    MFIXPIC& m_pic;
};

#endif
