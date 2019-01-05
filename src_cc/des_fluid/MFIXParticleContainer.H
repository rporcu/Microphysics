#ifndef _MFIXParticleContainer_H_
#define _MFIXParticleContainer_H_


#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_IArrayBox.H>
#include "AMReX_Particles.H"
#include "AMReX_NeighborParticles.H"

using namespace amrex;
using namespace std;

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
    // Data type returned by MFIXParIter::ArrayOfStructs::data (accesses
    // particle data represented by realData)
    //using RealType       = typename ParticleType::RealType;
    //using AoS            = ArrayOfStructs<realData::count, intData::count>;

    ///
    /// We add a particle to the neighbor list if its center is within 3 particle radii
    ///
    struct {
        bool operator()(const ParticleType& p1, const ParticleType& p2) const
        {
            return AMREX_D_TERM(   (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0))   ,
                                   + (p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)) ,
                                   + (p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)) )
                <= 9.0*p1.rdata(0)*p1.rdata(0);
        }
    } MFIXCheckPair;

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
                          const MultiFab * ls_phi, const iMultiFab * ls_valid, int ls_refinement);

    void Replicate(IntVect& Nrep, Geometry& geom, DistributionMapping& dmap, BoxArray& ba);

    void printParticles();

    void AllocData();

    void PrintParticleCounts();

    void EvolveParticles(int lev, int nstep, Real dt, Real time, EBFArrayBoxFactory * ebfactory,
            const MultiFab * ls_phi, const iMultiFab * ls_valid, const int ls_refinement,
            MultiFab * cost, std::string & knapsack_weight_type, int subdt_io);

    void GetParticleAvgProp(Real (&avg_dp)[10], Real (&avg_ro)[10]);

    //Real GetMaxRadius(int lev);

    void UpdateMaxVelocity();
    void UpdateMaxForces(std::map<PairIndex, Vector<Real>> pfor,
                         std::map<PairIndex, Vector<Real>> wfor);
    RealVect GetMaxVelocity();
    Vector<RealVect> GetMaxForces();

    void BalanceParticleLoad_KDTree();
    ///
    /// Write out all particles at a level in the format fcompare_par can read
    ///
    void writeAllForComparison(int lev);

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
                       int fortran_particle_comp, int nghost);

    ///
    /// Calculate particle-induced drag on the fluid
    ///
    void CalcDragOnFluid(const amrex::Vector< std::unique_ptr<amrex::MultiFab> > & beta_mf,
                         const amrex::Vector< std::unique_ptr<amrex::MultiFab> > & beta_vel_mf,
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
    void PICMultiDeposition(const amrex::Vector< std::unique_ptr<amrex::MultiFab> > & beta_mf,
                            const amrex::Vector< std::unique_ptr<amrex::MultiFab> > & beta_vel_mf,
                            const amrex::Vector< std::unique_ptr<EBFArrayBoxFactory> > & ebfactory,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
                            const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
                            int fortran_beta_comp, int fortran_vel_comp, int nghost);
    ///
    /// Flag switching from calc_wall_collisions to calc_wall_collisions_ls if true
    ///
    bool legacy__eb_collisions = false;

    //
    // Postprocessing utilities to compute spatially-averaged particle velocities
    // 
    void ComputeAverageVelocities ( const int lev,
				    const amrex::Real time,
				    const string& basename,
				    const vector<Real>& avg_region_x_w,
				    const vector<Real>& avg_region_x_e,
				    const vector<Real>& avg_region_y_s,
				    const vector<Real>& avg_region_y_n,
				    const vector<Real>& avg_region_z_b,
				    const vector<Real>& avg_region_z_t );

private:

    int nlev;

    ///
    /// We add a particle to the neighbor list if it's center is within 3 particle radii
    ///
    inline virtual bool check_pair(const ParticleType& p1, const ParticleType& p2) final {
        return AMREX_D_TERM(   (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0)) ,
                             + (p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)) ,
                             + (p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)) ) <= 9.0*p1.rdata(0)*p1.rdata(0);
    }

    static void ReadStaticParameters ();

    static bool use_neighbor_list;
    static bool sort_neighbor_list;

    // Local containers for: maxvel - the maximum (absolute) velocity
    //                       maxpfor - the max (abs) particle-particle force
    //                       maxwfor - the max (abs) particle-wall force
    RealVect loc_maxvel;
    RealVect loc_maxpfor;
    RealVect loc_maxwfor;

};
#endif
