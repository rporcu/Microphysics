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
    : public amrex::ParIter<realData::count, intData::count, 0, 0>
{
public:
    using amrex::ParIter<realData::count, intData::count, 0, 0>::ParIter;
};

class MFIXParticleContainer final
    : public amrex::NeighborParticleContainer<realData::count, intData::count>
{
public:
    using MyConstParIter = ParConstIter<realData::count, intData::count, 0, 0>;
    // Data type returned by MFIXParIter::ArrayOfStructs::data (accesses
    // particle data represented by realData)
    //using RealType       = typename ParticleType::RealType;
    //using AoS            = ArrayOfStructs<realData::count, intData::count>;

    // Get number of particles
    int NumberOfParticles(MFIXParIter& pti){
        return pti.GetArrayOfStructs().size();
    }

    MFIXParticleContainer(AmrCore* amr_core);

    virtual ~MFIXParticleContainer() {};

    void InitData();

    void InitParticlesAscii(const std::string& file);
    void InitParticlesAuto(int lev);
    void RemoveOutOfRange(int lev, EBFArrayBoxFactory * ebfactory);

    void Replicate(IntVect& Nrep, Geometry& geom, DistributionMapping& dmap, BoxArray& ba);

    void printParticles();

    void AllocData();

    void PrintParticleCounts();

    std::unique_ptr<MultiFab> EBNormals(int lev, EBFArrayBoxFactory * ebfactory, MultiFab * dummy);

    void EvolveParticles(int lev, int nstep, Real dt, Real time, EBFArrayBoxFactory * ebfactory,
            MultiFab * eb_normals, const MultiFab * ls_phi, const iMultiFab * ls_valid, const int ls_refinement,
            MultiFab * dummy, MultiFab * cost, std::string & knapsack_weight_type, int subdt_io);

    void GetParticleAvgProp(int lev, Real (&avg_dp)[10], Real (&avg_ro)[10]);

    //Real GetMaxRadius(int lev);

    void UpdateMaxVelocity(int lev);
    RealVect GetMaxVelocity(int lev);

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

    void output(int lev, int estatus, int finish, int nstep, Real dt, Real time);

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel(int lev = 0);

    ///
    /// Calculate gas volume fraction using a PIC deposition of particle volume
    ///
    void CalcVolumeFraction(amrex::MultiFab& mf,
                            IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                            IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                            IArrayBox& bc_klo, IArrayBox& bc_khi,
			    int nghost);

    ///
    /// Do the actual PIC deposition for CalcVolumeFraction
    ///
    void PICDeposition(amrex::MultiFab& mf,
                       IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                       IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                       IArrayBox& bc_klo, IArrayBox& bc_khi,
                       int fortran_particle_comp, int nghost );

    ///
    /// Calculate particle-induced drag on the fluid
    ///
    void CalcDragOnFluid(amrex::MultiFab& beta_mf_u,
                         amrex::MultiFab& beta_mf_v,
                         amrex::MultiFab& beta_mf_w,
                         amrex::MultiFab& beta_u_mf,
                         amrex::MultiFab& beta_v_mf,
                         amrex::MultiFab& beta_w_mf,
                         IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                         IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                         IArrayBox& bc_klo, IArrayBox& bc_khi,
			 int nghost );

    ///
    /// Do the actual PIC deposition for CalcDragOnFluid
    ///
    void PICMultiDeposition(amrex::MultiFab& beta_x_mf, amrex::MultiFab& beta_y_mf, amrex::MultiFab& beta_z_mf,
                            amrex::MultiFab& beta_u_mf, amrex::MultiFab& beta_v_mf, amrex::MultiFab& beta_w_mf,
                            IArrayBox& bc_ilo, IArrayBox& bc_ihi,
                            IArrayBox& bc_jlo, IArrayBox& bc_jhi,
                            IArrayBox& bc_klo, IArrayBox& bc_khi,
                            int fortran_beta_comp, int fortran_vel_comp, int nghost);

    ///
    /// Flag switching from calc_wall_collisions to calc_wall_collisions_ls if true
    ///
    bool legacy__eb_collisions = false;

private:

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

    // Local containers for: the maximum (absolute) velocity 
    //                       the max (abs) particle-particle force
    //                       the max (abs) particle-wall force
    RealVect loc_maxvel;
    RealVect loc_maxpfor;
    RealVect loc_maxwfor;
};
#endif
