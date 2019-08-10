#ifndef _MFIXParticleContainer_H_
#define _MFIXParticleContainer_H_

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_IArrayBox.H>
#include "AMReX_Particles.H"
#include "AMReX_NeighborParticles.H"

#include "mfix_des_F.H"
#include "MFIX_DEM_Parms.H"
#include <particles_generator.hpp>

// Boundary conditions list auxiliary structure
//
typedef class BcList {
  public:
    BcList(const int pinf, const int pout, const int minf)
      : pinf_(pinf)
      , pout_(pout) 
      , minf_(minf)
    {}

    int get_pinf() const {return pinf_;};
    int get_pout() const {return pout_;};
    int get_minf() const {return minf_;};
  private:
    int pinf_, pout_, minf_;
} BcList;

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

///
/// We add a particle to the neighbor list if its center is within 3 particle radii
///
struct MFIXCheckPair
{
  template <class P>
  AMREX_GPU_DEVICE AMREX_INLINE
  bool operator()(const P& p1, const P& p2) const
  {
    return AMREX_D_TERM( (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0)),
                        +(p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)),
                        +(p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)))
             <= DEMParams::neighborhood;
  }
};

class MFIXParIter
    : public amrex::ParIter<realData::count,intData::count,0,0>
{
public:
    using amrex::ParIter<realData::count,intData::count,0,0>::ParIter;
};

class MFIXParticleContainer final
    : public amrex::NeighborParticleContainer<realData::count,intData::count>
{
public:
    using MyConstParIter = amrex::ParConstIter<realData::count,intData::count,0,0>;

    // Get number of particles
    int NumberOfParticles(MFIXParIter& pti){
        return pti.GetArrayOfStructs().size();
    }

    MFIXParticleContainer(amrex::AmrCore* amr_core);

    virtual ~MFIXParticleContainer() {};

    void InitData();

    void InitParticlesAscii(const std::string& file);
    void InitParticlesAuto();
    void RemoveOutOfRange(int lev,
                          const amrex::EBFArrayBoxFactory * ebfactory,
                          const amrex::MultiFab * ls_phi,
                          int ls_refinement);

    void Replicate(amrex::IntVect& Nrep,
                   amrex::Geometry& geom,
                   amrex::DistributionMapping& dmap,
                   amrex::BoxArray& ba);

    void printParticles();

    void AllocData();

    void PrintParticleCounts();

    void EvolveParticles(int lev,
                         int nstep,
                         amrex::Real dt,
                         amrex::Real time,
                         amrex::EBFArrayBoxFactory * ebfactory,
                         const amrex::MultiFab * ls_phi,
                         const amrex::iMultiFab * ls_valid,
                         const int ls_refinement,
                         amrex::MultiFab * cost,
                         std::string & knapsack_weight_type);

    void GetParticleAvgProp(amrex::Real (&min_dp)[10], amrex::Real (&min_ro)[10],
                            amrex::Real (&max_dp)[10], amrex::Real (&max_ro)[10],
                            amrex::Real (&avg_dp)[10], amrex::Real (&avg_ro)[10]);

    //amrex::Real GetMaxRadius(int lev);

    void UpdateMaxVelocity();
    void UpdateMaxForces(std::map<PairIndex, amrex::Gpu::ManagedDeviceVector<amrex::Real>> pfor,
                         std::map<PairIndex, amrex::Gpu::ManagedDeviceVector<amrex::Real>> wfor);

    amrex::RealVect GetMaxVelocity();
    amrex::Vector<amrex::RealVect> GetMaxForces();

    void BalanceParticleLoad_KDTree();

    ///
    /// Write out all particles in the format of the particle_input.dat files so we
    ///       can initialize a run with them
    ///
    void WriteAsciiFileForInit(const std::string& filename);

    void output(int estatus,
                int finish,
                int nstep,
                amrex::Real dt,
                amrex::Real time);

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel(int lev = 0);

    ///
    /// Calculate gas volume fraction using a PIC deposition of particle volume
    ///
    void CalcVolumeFraction(const amrex::Vector< std::unique_ptr<amrex::MultiFab> > & mf,
                            const amrex::Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > & ebfactory,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
                            const BcList & bc_list,
                            int nghost);

    void CapSolidsVolFrac(amrex::MultiFab& eps);

    ///
    /// Do the actual PIC deposition for CalcVolumeFraction
    ///
    void PICDeposition(const amrex::Vector< std::unique_ptr<amrex::MultiFab> > & mf,
                       const amrex::Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > & ebfactory,
                       const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
                       const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
                       const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
                       const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
                       const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
                       const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
                       const BcList & bc_list,
                       int fortran_particle_comp,
                       int nghost);

    ///
    /// Calculate particle-induced drag on the fluid
    ///
    void CalcDragOnFluid(const amrex::Vector< std::unique_ptr<amrex::MultiFab> > & drag_mf,
                         const amrex::Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > & ebfactory,
                         const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
                         const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
                         const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
                         const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
                         const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
                         const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
                         int nghost);

    ///
    /// Do the actual PIC deposition for CalcDragOnFluid
    ///
    void PICMultiDeposition(const amrex::Vector< std::unique_ptr<amrex::MultiFab> > & drag_mf,
                            const amrex::Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > & ebfactory,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
                            int fortran_beta_comp,
                            int fortran_vel_comp,
                            int nghost);

    //
    // Postprocessing utilities to compute spatially-averaged particle velocities
    //
    void ComputeAverageVelocities(const int lev,
                                  const amrex::Real time,
                                  const std::string& basename,
                                  const amrex::Vector<int>& avg_vel_p,
                                  const amrex::Gpu::ManagedDeviceVector<amrex::Real>& avg_region_x_w,
                                  const amrex::Gpu::ManagedDeviceVector<amrex::Real>& avg_region_x_e,
                                  const amrex::Gpu::ManagedDeviceVector<amrex::Real>& avg_region_y_s,
                                  const amrex::Gpu::ManagedDeviceVector<amrex::Real>& avg_region_y_n,
                                  const amrex::Gpu::ManagedDeviceVector<amrex::Real>& avg_region_z_b,
                                  const amrex::Gpu::ManagedDeviceVector<amrex::Real>& avg_region_z_t );

private:
    int nlev;

    static amrex::Real gravity[3];
    static int domain_bc[6];

    // Particles generator
    ParticlesGenerator m_particles_generator;

    ///
    /// We add a particle to the neighbor list if it's center is within 3 particle radii
    ///
    inline virtual bool check_pair(const ParticleType& p1, const ParticleType& p2) final {
        return AMREX_D_TERM( (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0)),
                            +(p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)),
                            +(p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)))
                 <= 9.0*p1.rdata(0)*p1.rdata(0);
    }

    static void ReadStaticParameters ();

    // Local containers for: maxvel - the maximum (absolute) velocity
    //                       maxpfor - the max (abs) particle-particle force
    //                       maxwfor - the max (abs) particle-wall force
    amrex::RealVect loc_maxvel;
    amrex::RealVect loc_maxpfor;
    amrex::RealVect loc_maxwfor;

    void time_advance(MFIXParIter& pti,
                      int ntot,
                      amrex::Real subdt,
                      amrex::Gpu::ManagedDeviceVector<amrex::Real>& tow,
                      amrex::Gpu::ManagedDeviceVector<amrex::Real>& fc);

    void set_particle_properties(int pstate,
                                 amrex::Real pradius,
                                 amrex::Real pdensity,
                                 amrex::Real& pvol,
                                 amrex::Real& pmass,
                                 amrex::Real& omoi,
                                 amrex::Real& omega);
};
#endif
