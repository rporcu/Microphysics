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
#include <mfix_fluid_parms.H>
#include <mfix_solids_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_bc_list.H>
#include <mfix_leveldata.H>


struct AoSrealData
{
  enum {
    count = 0
  };
};

struct AoSintData
{
  enum {
    count = 0
  };
};

struct SoArealData
{
  enum {
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
    cp_s,
    temperature,
    convection,
    count
  };
};

// Particle real data stored at runtime in a SoA style
struct runtimeRealData
{
  runtimeRealData (const int nspecies_s,
                   const int nreactions)
    : X_sn(0)                                              // nb of components = nspecies
    , mass_sn_txfr(X_sn + nspecies_s)                      // nb of components = nspecies
    , vel_s_txfr(mass_sn_txfr + nspecies_s*(nreactions>0)) // nb of components = 3
    , h_s_txfr(vel_s_txfr + 3*(nreactions>0))              // nb of components = 1
    , count(h_s_txfr + 1*(nreactions>0))
  {}

  const int X_sn;
  const int mass_sn_txfr;
  const int vel_s_txfr;
  const int h_s_txfr;
  const int count;
};


struct SoAintData
{
  enum {
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
  bool operator() (const P* particles, const int i, const int j) const
  {
    const P& p1 = particles[i];
    const P& p2 = particles[j];

    if (i < j)
      if ((p1.pos() - p2.pos()).radSquared() <= neighborhood)
        return true;

    return false;
  }
};

struct MFIXCheckFullPair
{
  amrex::Real neighborhood;
  MFIXCheckFullPair (amrex::Real neighborhood_) : neighborhood(neighborhood_) {}

  template <class P>
  AMREX_GPU_DEVICE AMREX_INLINE
  bool operator() (const P* particles, const int i, const int j) const
  {
    const P& p1 = particles[i];
    const P& p2 = particles[j];

    if ((p1.pos() - p2.pos()).radSquared() <= neighborhood)
      return true;

    return false;
  }
};

class MFIXParticleContainer final
    : public amrex::NeighborParticleContainer<AoSrealData::count,AoSintData::count,
                                              SoArealData::count,SoAintData::count>
{
public:
    using MFIXParIter = amrex::ParIter<AoSrealData::count,AoSintData::count,
                                       SoArealData::count,SoAintData::count>;
    using MFIXParConstIter = amrex::ParConstIter<AoSrealData::count,AoSintData::count,
                                                 SoArealData::count,SoAintData::count>;

    // Get number of particles
    int NumberOfParticles (MFIXParIter& pti)
    { return pti.numParticles(); }

    void setTotalNumParticles(const int num_particles)
    {
      this->m_total_numparticle = num_particles;
    }

    MFIXParticleContainer (amrex::AmrCore* amr_core,
                           SolidsPhase& solids,
                           Reactions& reactions);

    virtual ~MFIXParticleContainer () {};

    void InitData ();

    void InitParticlesAscii (const std::string& file);

    void InitParticlesAuto ();

    void InitParticlesRuntimeVariables (const int adv_enthalpy);

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

    amrex::Vector<std::map<std::pair<int, int>, amrex::NeighborList<ParticleType> > > get_neighbor_list()
    { return m_neighbor_list; }

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
                          const amrex::Real enthalpy_source,
                          const int update_mass,
                          const int update_momentum,
                          const int update_enthalpy);

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


    //
    // Postprocessing utilities to compute spatially-averaged particle velocities
    //
    void ComputeAverageDensities (const int lev,
                                  const amrex::Real time,
                                  const std::string& basename,
                                  const amrex::Vector<amrex::Real>& avg_ro_p,
                                  const amrex::Vector<amrex::Real>& avg_region_x_w,
                                  const amrex::Vector<amrex::Real>& avg_region_x_e,
                                  const amrex::Vector<amrex::Real>& avg_region_y_s,
                                  const amrex::Vector<amrex::Real>& avg_region_y_n,
                                  const amrex::Vector<amrex::Real>& avg_region_z_b,
                                  const amrex::Vector<amrex::Real>& avg_region_z_t);

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

    ///
    /// Deposit parcel velocity to the grid
    ///
    void MFIX_PC_SolidsVelocityDeposition (int lev,
                                           amrex::Array<amrex::MultiFab*,3>& vel_s_mf,
                                           const amrex::FabArray<amrex::EBCellFlagFab>* flags);

    void PICHydroStep (int lev,
                       const bool apply_forces,
                       const bool update_parcels,
                       const bool use_taylor_approx,
                       const amrex::Real advance_vel_p,
                       amrex::Real dt,
                       amrex::RealVect& gravity,
                       amrex::Vector< amrex::Array<amrex::MultiFab*,3> >& vel_s_in,
                       amrex::MultiFab & ep_s_out,
                       amrex::Array<amrex::MultiFab*,3>& vel_s_out,
                       const amrex::MultiFab * volfrac,
                       const amrex::FabArray<amrex::EBCellFlagFab>* flags,
                       amrex::EBFArrayBoxFactory* ebfactory,
                       const int ls_refinement,
                       const amrex::MultiFab* ls_phi);

    void MFIX_PC_AdvanceParcels (amrex::Real dt,
                                 const int advect_enthalpy,
                                 const amrex::Real enthalpy_source,
                                 amrex::Vector< amrex::MultiFab*> & cost,
                                 std::string& knapsack_weight_type);

    const runtimeRealData m_runtimeRealData;

    void setSortingBinSizes (amrex::IntVect sorting_bin) { m_sorting_bin = sorting_bin;  }

    void partitionParticleGrids(int lev,
                                const amrex::BoxArray& fba,
                                const amrex::DistributionMapping& fdmap,
                                const int greedy_dir,
                                amrex::Real imbalance_toler,
                                amrex::Real partition_factor);

    void setMaxGridSize(const amrex::IntVect& rhs)
    { max_grid_size = rhs;}

    void setMaxGridSize(int size_x, int size_y, int size_z)
    { max_grid_size = amrex::IntVect(AMREX_D_DECL(size_x, size_y, size_z));}

    void setSortInt(int sort_particle_int)
    { sort_int = sort_particle_int; }

    void setReduceGhostParticles(bool reduce)
    { reduceGhostParticles = reduce; };

    const amrex::IntVect& MaxGridSize() const noexcept
    { return max_grid_size; }

    // return the average/max particle count on IO process
    amrex::Real particleImbalance();

    /* count particles for each box in the box list
       The boxes are broken down into bins of the input size
       output: pcounts  - # particles in each bin
               poffsets - offsets for bins in each box
    */
    void countParticle(int lev,
                       const amrex::BoxList& bl,
                       const amrex::IntVect& binsize,
                       amrex::Vector<int>& pcounts,
                       amrex::Vector<int>& poffsets);

    // verify total # particles hasn't been changed, for debugging
    void verifyParticleCount();

    void printGhostParticleCount();

    void ReadParameters();

  private:
    int nlev;

    int sort_int = -1; // # steps to sort particles

    bool reduceGhostParticles = false;

    amrex::Real newton_abstol = 1.e-6;
    amrex::Real newton_reltol = 1.e-6;
    int newton_maxiter = 100;

    amrex::IntVect max_grid_size{-1};

    amrex::IntVect m_sorting_bin = amrex::IntVect::TheZeroVector();

    static int domain_bc[6];

    int m_total_numparticle = 0;
    amrex::Vector<int> m_pboxid_to_fboxid;

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

    void usr2_des (int np, ParticleTileType& particles);

    void usr3_des (int np,void*&);

    SolidsPhase& solids;

    Reactions& reactions;
};

#endif
