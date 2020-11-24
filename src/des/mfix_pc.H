#ifndef _MFIX_PC_H_
#define _MFIX_PC_H_

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_FabArray.H>
#include <AMReX_EBSupport.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FabArray.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <mfix_pic_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_bc_list.H>

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
    statwt,
    dragcoeff,
    dragx,
    dragy,
    dragz,
    c_ps,
    temperature,
    convection,
    count
  };
};

struct speciesData
{
  enum { // Added Particle data stored in p.m_runtime_data
    X_sn = 0,
    count
  };
};

struct reactionsData
{
  enum { // Added Particle data stored in p.m_runtime_data
    G_sn_pg_q = 0,
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
  amrex::Real neighborhood;
  MFIXCheckPair (amrex::Real neighborhood_) : neighborhood(neighborhood_) {}

  template <class P>
  AMREX_GPU_DEVICE AMREX_INLINE
  bool operator() (const P& p1, const P& p2) const
  { return (p1.pos() - p2.pos()).radSquared() <= neighborhood; }
};

class MFIXParticleContainer final
    : public amrex::NeighborParticleContainer<realData::count,intData::count>
{
public:
    using MFIXParIter = amrex::ParIter<realData::count,intData::count,0,0>;

    // Get number of particles
    int NumberOfParticles (MFIXParIter& pti)
    { return pti.GetArrayOfStructs().size(); }

    MFIXParticleContainer (amrex::AmrCore* amr_core);

    virtual ~MFIXParticleContainer () {};

    void InitData ();

    void InitParticlesAscii (const std::string& file);

    void InitParticlesAuto ();

    void InitParticlesEnthalpy ();

    void InitParticlesSpecies ();

    void RemoveOutOfRange (int lev,
                           const amrex::EBFArrayBoxFactory * ebfactory,
                           const amrex::MultiFab * ls_phi,
                           int ls_refinement);

    void Replicate (amrex::IntVect& Nrep,
                    amrex::Geometry& geom,
                    amrex::DistributionMapping& dmap,
                    amrex::BoxArray& ba);

    void printParticles ();

    void AllocData ();

    void PrintParticleCounts ();

    void EvolveParticles (int lev,
                          int nstep,
                          amrex::Real dt,
                          amrex::Real time,
                          amrex::RealVect& gravity,
                          amrex::EBFArrayBoxFactory * ebfactory,
                          const amrex::MultiFab * ls_phi,
                          const int ls_refinement,
                          amrex::MultiFab * cost,
                          std::string & knapsack_weight_type,
                          int & nsubsteps,
                          const int advect_enthalpy,
                          const amrex::Real enthalpy_source);

    void MFIX_PC_InitCollisionParams ();

    void UpdateMaxVelocity ();

    void UpdateMaxForces (std::map<PairIndex, amrex::Gpu::DeviceVector<amrex::Real>>& pfor,
                          std::map<PairIndex, amrex::Gpu::DeviceVector<amrex::Real>>& wfor);

    amrex::RealVect GetMaxVelocity ();

    amrex::Vector<amrex::RealVect> GetMaxForces ();

    ///
    /// Write out all particles in the format of the particle_input.dat files so we
    ///       can initialize a run with them
    ///
    void WriteAsciiFileForInit (const std::string& filename);

    void output (int estatus,
                 int finish,
                 int nstep,
                 amrex::Real dt,
                 amrex::Real time);

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel (int lev = 0);

    void SolidsVolumeDeposition (int lev,
                                 amrex::MultiFab& mf,
                                 const amrex::MultiFab* volfrac,
                                 const amrex::FabArray<amrex::EBCellFlagFab>* flags);

    ///
    /// Do the actual PIC deposition for CalcVolumeFraction
    ///
    template <typename FB>
    void SolidsVolumeDeposition (FB WeightFunc,
                                 int lev,
                                 amrex::MultiFab& mf,
                                 const amrex::MultiFab* volfrac,
                                 const amrex::FabArray<amrex::EBCellFlagFab>* flags);

    ///
    /// Do the actual PIC deposition for CalcDragOnFluid
    ///
    void InterphaseTxfrDeposition (int lev,
                                   amrex::MultiFab& mf_tmp_eps,
                                   amrex::MultiFab& txfr_mf,
                                   const amrex::MultiFab* volfrac,
                                   const amrex::FabArray<amrex::EBCellFlagFab>* flags,
                                   const int advect_enthalpy);

    template <typename F>
    void InterphaseTxfrDeposition (F WeightFunc,
                                   int lev,
                                   amrex::MultiFab& mf_tmp_eps,
                                   amrex::MultiFab& txfr_mf,
                                   const amrex::MultiFab* volfrac,
                                   const amrex::FabArray<amrex::EBCellFlagFab>* flags,
                                   const int advect_enthalpy);

    ///
    /// Do the actual PIC deposition for CalcDragOnFluid
    ///
    void InterphaseChemDeposition (int lev,
                                   amrex::MultiFab& mf_tmp_eps,
                                   amrex::MultiFab& Rrates,
                                   const amrex::MultiFab* volfrac,
                                   const amrex::FabArray<amrex::EBCellFlagFab>* flags,
                                   const amrex::Vector< REACTIONS::ChemicalReaction* >&);

    template <typename F>
    void InterphaseChemDeposition (F WeightFunc,
                                   int lev,
                                   amrex::MultiFab& mf_tmp_eps,
                                   amrex::MultiFab& chem_txfr_mf,
                                   const amrex::MultiFab* volfrac,
                                   const amrex::FabArray<amrex::EBCellFlagFab>* flags,
                                   const amrex::Vector< REACTIONS::ChemicalReaction* >&);

    //
    // Postprocessing utilities to compute spatially-averaged particle velocities
    //
    void ComputeAverageVelocities (const int lev,
                                   const amrex::Real time,
                                   const std::string& basename,
                                   const amrex::Vector<amrex::Real>& avg_vel_p,
                                   const amrex::Vector<amrex::Real>& avg_region_x_w,
                                   const amrex::Vector<amrex::Real>& avg_region_x_e,
                                   const amrex::Vector<amrex::Real>& avg_region_y_s,
                                   const amrex::Vector<amrex::Real>& avg_region_y_n,
                                   const amrex::Vector<amrex::Real>& avg_region_z_b,
                                   const amrex::Vector<amrex::Real>& avg_region_z_t);


    //
    // Postprocessing utilities to compute spatially-averaged particle velocities
    //
    void ComputeAverageTemperatures (const int lev,
                                     const amrex::Real time,
                                     const std::string& basename,
                                     const amrex::Vector<amrex::Real>& avg_T_p,
                                     const amrex::Vector<amrex::Real>& avg_region_x_w,
                                     const amrex::Vector<amrex::Real>& avg_region_x_e,
                                     const amrex::Vector<amrex::Real>& avg_region_y_s,
                                     const amrex::Vector<amrex::Real>& avg_region_y_n,
                                     const amrex::Vector<amrex::Real>& avg_region_z_b,
                                     const amrex::Vector<amrex::Real>& avg_region_z_t);



    void MFIX_PC_CalcSolidsStressGrad (int lev,
                                       amrex::EBFArrayBoxFactory* ebfactory,
                                       amrex::MultiFab & ep_g_in,
                                       amrex::MultiFab & grad_tau_in,
                                       amrex::Vector< amrex::MultiFab* >& cost,
                                       std::string& knapsack_weight_type);

    ///
    /// Deposit parcel velocity to the grid
    ///
    void MFIX_PC_SolidsVelocityDeposition (int lev,
                                           amrex::MultiFab& vel_s_mf,
                                           const amrex::FabArray<amrex::EBCellFlagFab>* flags);


    void MFIX_PC_AdvanceParcels (amrex::Real dt, amrex::RealVect& gravity,
                                 amrex::Vector< amrex::MultiFab*> & avg_prop_in,
                                 amrex::Vector< amrex::MultiFab*> & cost,
                                 std::string& knapsack_weight_type,
                                 const int advect_enthalpy);

    void MFIX_PC_ImposeWalls (int lev,
                              amrex::EBFArrayBoxFactory * ebfactory,
                              const int ls_refinement,
                              const amrex::MultiFab * ls_phi,
                              amrex::MultiFab* cost,
                              std::string& knapsack_weight_type);

  private:
    int nlev;

    static int domain_bc[6];

    static void ReadStaticParameters ();

    // Local containers for: maxvel - the maximum (absolute) velocity
    //                       maxpfor - the max (abs) particle-particle force
    //                       maxwfor - the max (abs) particle-wall force
    amrex::RealVect loc_maxvel;
    amrex::RealVect loc_maxpfor;
    amrex::RealVect loc_maxwfor;

    void set_particle_properties (int pstate,
                                  amrex::Real pradius,
                                  amrex::Real pdensity,
                                  amrex::Real& pvol,
                                  amrex::Real& pmass,
                                  amrex::Real& omoi,
                                  amrex::Real& omega);

    void usr0_des ();
    void usr1_des ();

    void usr2_des (int np,
                   NeighborParticleContainer<realData::count, intData::count>::ParticleType*& particles);

    void usr3_des (int np,void*&);
};

#endif
