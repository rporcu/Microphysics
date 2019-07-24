#ifndef _PARTICLES_GENERATOR_HPP_
#define _PARTICLES_GENERATOR_HPP_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_CudaContainers.H>

class ParticlesGenerator {
  //< Position............... 1,2,3
  //< Radius................. 4
  //< Density................ 5
  //< Linear velocity........ 6,7,8
  static constexpr int nr = 8;

  //< Type................... 1
  static constexpr int ni = 1;

  public:
    // Constructor
    explicit ParticlesGenerator();

    // Destructor
    ~ParticlesGenerator();

    // Particle generator
    void generate(int& pc,
                  int iter,
                  const amrex::IntVect& lo,
                  const amrex::IntVect& hi,
                  const amrex::Real dx,
                  const amrex::Real dy,
                  const amrex::Real dz);

    void generate_prop(const int nrp,
                       void* particles);

    // Write particles data
    void write(const int nrp,
               void* particles) const;

    void hex_close_pack(const int icv,
                        const int type,
                        const amrex::IntVect& lo,
                        const amrex::IntVect& hi,
                        int& np,
                        int& pc,
                        const amrex::Real dx,
                        const amrex::Real dy,
                        const amrex::Real dz);

    void one_per_fill(const int icv,
                      const int type,
                      const amrex::IntVect& lo,
                      const amrex::IntVect& hi,
                      int& np,
                      int& pc,
                      const amrex::Real dx,
                      const amrex::Real dy,
                      const amrex::Real dz);

    void eight_per_fill(const int icv,
                        const int type,
                        const amrex::IntVect& lo,
                        const amrex::IntVect& hi,
                        int& np,
                        int& pc,
                        const amrex::Real dx,
                        const amrex::Real dy,
                        const amrex::Real dz);

    void random_fill(const int icv,
                     const int type,
                     const amrex::IntVect& lo,
                     const amrex::IntVect& hi,
                     int& np,
                     int& pc,
                     const amrex::Real dx,
                     const amrex::Real dy,
                     const amrex::Real dz,
                     const bool fix_seed);

    void nor_rno(amrex::Cuda::ManagedVector<amrex::Real>& dp,
                 const amrex::Real mean,
                 const amrex::Real sigma,
                 const amrex::Real dp_min,
                 const amrex::Real dp_max);

    void uni_rno(amrex::Cuda::ManagedVector<amrex::Real>& dp,
                 const amrex::Real dp_min,
                 const amrex::Real dp_max);

    void grow_pdata(const int gsize);

  private:
    amrex::Cuda::ManagedVector<amrex::Real> m_rdata;
    amrex::Cuda::ManagedVector<int> m_idata;
};

#endif
