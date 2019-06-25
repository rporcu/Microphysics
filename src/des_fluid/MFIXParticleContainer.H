#ifndef _MFIXParticleContainer_H_
#define _MFIXParticleContainer_H_


#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_IArrayBox.H>
#include "AMReX_Particles.H"
#include "AMReX_NeighborParticles.H"

#include "mfix_des_F.H"
#include "MFIX_DEM_Parms.H"

using namespace amrex;
using namespace std;

// Boundary conditions list auxiliary structure
//
typedef struct BcList {
  public:
    BcList(): minf(0), pinf(0), pout(0) {};
    int minf;
    int pinf;
    int pout;
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
        return AMREX_D_TERM(   (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0)) ,
                             + (p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)) ,
                             + (p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)) )
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

    using MyConstParIter = ParConstIter<realData::count,intData::count,0,0>;

    // Get number of particles
    int NumberOfParticles(MFIXParIter& pti){
        return pti.GetArrayOfStructs().size();
    }

    MFIXParticleContainer(AmrCore* amr_core);

    virtual ~MFIXParticleContainer() {};

    void InitData();

    void InitParticlesAscii(const std::string& file);
    void InitParticlesAuto();
    void RemoveOutOfRange(int lev, const EBFArrayBoxFactory * ebfactory,
                          const MultiFab * ls_phi, int ls_refinement);

    void Replicate(IntVect& Nrep, Geometry& geom, DistributionMapping& dmap, BoxArray& ba);

    void printParticles();

    void AllocData();

    void PrintParticleCounts();

    void EvolveParticles(int lev, int nstep, Real dt, Real time, EBFArrayBoxFactory * ebfactory,
            const MultiFab * ls_phi, const iMultiFab * ls_valid, const int ls_refinement,
            MultiFab * cost, std::string & knapsack_weight_type);

    void GetParticleAvgProp(Real (&min_dp)[10], Real (&min_ro)[10],
                            Real (&max_dp)[10], Real (&max_ro)[10],
                            Real (&avg_dp)[10], Real (&avg_ro)[10]);

    //Real GetMaxRadius(int lev);

    void UpdateMaxVelocity();
    void UpdateMaxForces(std::map<PairIndex, Gpu::ManagedDeviceVector<Real>> pfor,
                         std::map<PairIndex, Gpu::ManagedDeviceVector<Real>> wfor);
    RealVect GetMaxVelocity();
    Vector<RealVect> GetMaxForces();

    void BalanceParticleLoad_KDTree();

    ///
    /// Write out all particles in the format of the particle_input.dat files so we
    ///       can initialize a run with them
    ///
    void WriteAsciiFileForInit(const std::string& filename);

    void output(int estatus, int finish, int nstep, Real dt, Real time);

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel(int lev = 0);

    ///
    /// Calculate gas volume fraction using a PIC deposition of particle volume
    ///
    void CalcVolumeFraction(const amrex::Vector< std::unique_ptr<amrex::MultiFab> > & mf,
                            const amrex::Vector< std::unique_ptr<EBFArrayBoxFactory> > & ebfactory,
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
                       const amrex::Vector< std::unique_ptr<EBFArrayBoxFactory> > & ebfactory,
                       const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
                       const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
                       const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
                       const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
                       const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
                       const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
                       const BcList & bc_list,
                       int fortran_particle_comp, int nghost);

    ///
    /// Calculate particle-induced drag on the fluid
    ///
    void CalcDragOnFluid(const amrex::Vector< std::unique_ptr<amrex::MultiFab> > & drag_mf,
                         const amrex::Vector< std::unique_ptr<EBFArrayBoxFactory> > & ebfactory,
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
                            const amrex::Vector< std::unique_ptr<EBFArrayBoxFactory> > & ebfactory,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
                            int fortran_beta_comp, int fortran_vel_comp, int nghost);

    //
    // Postprocessing utilities to compute spatially-averaged particle velocities
    //
    void ComputeAverageVelocities ( const int lev,
            const amrex::Real time,
            const string& basename,
            const Vector<int>& avg_vel_p,
            const Gpu::ManagedDeviceVector<Real>& avg_region_x_w,
            const Gpu::ManagedDeviceVector<Real>& avg_region_x_e,
            const Gpu::ManagedDeviceVector<Real>& avg_region_y_s,
            const Gpu::ManagedDeviceVector<Real>& avg_region_y_n,
            const Gpu::ManagedDeviceVector<Real>& avg_region_z_b,
            const Gpu::ManagedDeviceVector<Real>& avg_region_z_t );

private:

    int nlev;

    static Real gravity[3];
    static int domain_bc[6];

    ///
    /// We add a particle to the neighbor list if it's center is within 3 particle radii
    ///
    inline virtual bool check_pair(const ParticleType& p1, const ParticleType& p2) final {
        return AMREX_D_TERM(   (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0)) ,
                             + (p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)) ,
                             + (p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)) ) <= 9.0*p1.rdata(0)*p1.rdata(0);
    }

    static void ReadStaticParameters ();

    // Local containers for: maxvel - the maximum (absolute) velocity
    //                       maxpfor - the max (abs) particle-particle force
    //                       maxwfor - the max (abs) particle-wall force
    RealVect loc_maxvel;
    RealVect loc_maxpfor;
    RealVect loc_maxwfor;

    void time_advance(MFIXParIter& pti, int ntot, Real subdt,
                      Gpu::ManagedDeviceVector<Real>& tow, Gpu::ManagedDeviceVector<Real>& fc);

    void set_particle_properties(int pstate, Real pradius, Real pdensity,
                                 Real& pvol, Real& pmass, Real& omoi, Real& omega);


};
#endif
