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

#include <mfix_ic.H>
#include <mfix_bc.H>
#include <mfix_pic.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_leveldata.H>
#include <mfix_run_on.H>


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
class runtimeRealData
{
  public:
    AMREX_GPU_HOST_DEVICE
    runtimeRealData (const int nspecies_s,
                     const int nreactions)
      : X_sn(0)                                       // nb of components = nspecies_s
      , vel_txfr(X_sn + nspecies_s)                   // nb of components = 3
      , h_txfr(vel_txfr + 3*(nreactions>0))           // nb of components = 1
      , mass_txfr(h_txfr + 1*(nreactions>0))          // nb of components = nspecies_s
      , count(mass_txfr + nspecies_s*(nreactions>0))
    {}

    const int X_sn;
    const int vel_txfr;
    const int h_txfr;
    const int mass_txfr;
    const int count;
};


struct SoAintData
{
  enum {
    phase = 0,
    state,
#if MFIX_POLYDISPERSE
    ptype,
#endif
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

  template <int NSR, int NSI, int NAR, int NAI>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  bool operator() (const amrex::ConstParticleTileData<amrex::Particle<NSR, NSI>, NAR, NAI>& ptile, int i, int j) const
  {
    const auto p1 = ptile.getSuperParticle(i);
    const auto p2 = ptile.getSuperParticle(j);

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

  template <int NSR, int NSI, int NAR, int NAI>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  bool operator() (const amrex::ConstParticleTileData<amrex::Particle<NSR, NSI>, NAR, NAI>& ptile, int i, int j) const
  {
    const auto p1 = ptile.getSuperParticle(i);
    const auto p2 = ptile.getSuperParticle(j);

    if ((p1.pos() - p2.pos()).radSquared() <= neighborhood)
      return true;

    return false;
  }
};


///
/// Neighbor search with multi-grid for polydisperse flows
///
struct MFIXCheckPolyPair
{
    int typeInd;
    int num_bin_types;
    amrex::Real* neighborhood;
    MFIXCheckPolyPair (int typeInd_, int num_bin_types_, amrex::Real* neighborhood_)
        : typeInd(typeInd_), num_bin_types(num_bin_types_), neighborhood(neighborhood_) {}

  template <int NSR, int NSI, int NAR, int NAI>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  bool operator() (const amrex::ConstParticleTileData<amrex::Particle<NSR, NSI>, NAR, NAI>& ptile,
                   int i, int j, int type, bool ghost_i, bool ghost_j) const
  {
    const auto p1 = ptile.getSuperParticle(i);
    const auto p2 = ptile.getSuperParticle(j);

    int type_i = p1.idata(typeInd);
    int type_j = p2.idata(typeInd);

    if (!(ghost_i && ghost_j)) {                    // One particle must be a real particle
        if ((type_i == type) || (type_j == type)) { // One of the particles matches the type being searched
            if ((i < j) || (type_i < type_j)) {     // Smaller index or type particle finds the neighbor
                // Get index for NEIGHBORHOOD from analytical partial sums
                int imin = amrex::min(type_i,type_j);
                int imax = amrex::max(type_i,type_j);
                int k    = imin - 1;
                // \Sum_{i=N-k}^{N} (i)
                int indx = num_bin_types*(k+1) - (k*(k+1))/2;
                // \Sum_{i=1}^{l} (i)
                indx += imax - imin;
                if ( (p1.pos() - p2.pos()).radSquared() <= neighborhood[indx] ) return true;
            }
        }
    }

    return false;
  }
};


class MFIXParticleContainer final
    : public amrex::NeighborParticleContainer<0,0,SoArealData::count,SoAintData::count>
{
  public:
    using MFIXParIter = amrex::ParIter<0,0,SoArealData::count,SoAintData::count>;
    using MFIXParConstIter = amrex::ParConstIter<0,0,SoArealData::count,SoAintData::count>;

    // Get number of particles
    int NumberOfParticles (MFIXParIter& pti)
    { return pti.numParticles(); }

    void setTotalNumParticles (const long num_particles)
    { this->m_total_numparticle = num_particles; }

    long getTotalNumParticles () { return m_total_numparticle; }

    MFIXParticleContainer (amrex::AmrCore* amr_core,
                           MFIXInitialConditions& initial_conditions,
                           MFIXBoundaryConditions& boundary_conditions,
                           MFIXSolidsPhase& solids,
                           MFIXDEM& dem,
                           MFIXPIC& pic,
                           MFIXFluidPhase& fluid,
                           MFIXReactions& reactions);

    MFIXParticleContainer (const amrex::Geometry& geom,
                           const amrex::DistributionMapping& dmap,
                           const amrex::BoxArray& ba,
                           const int nlev,
                           MFIXInitialConditions& initial_conditions,
                           MFIXBoundaryConditions& boundary_conditions,
                           MFIXSolidsPhase& solids,
                           MFIXDEM& dem,
                           MFIXPIC& pic,
                           MFIXFluidPhase& fluid,
                           MFIXReactions& reactions);

    void define();

    virtual ~MFIXParticleContainer () {}

    void InitData ();

    void InitParticlesAscii (const std::string& file);

    void InitParticlesAuto (amrex::EBFArrayBoxFactory* particle_ebfactory);

    void InitParticlesRuntimeVariables (const int adv_enthalpy);

    void RemoveOutOfRange (int lev,
                           const amrex::EBFArrayBoxFactory * ebfactory,
                           const amrex::MultiFab * ls_phi,
                           int ls_refinement);

    void ReportParticleGenerationStats (int lev);

    void Replicate (amrex::IntVect& Nrep,
                    amrex::Geometry& geom,
                    amrex::DistributionMapping& dmap,
                    amrex::BoxArray& ba);

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
                          amrex::EBFArrayBoxFactory * particle_ebfactory,
                          const amrex::MultiFab * ls_phi,
                          const int ls_refinement,
                          amrex::MultiFab * cost,
                          std::string & knapsack_weight_type,
                          int & nsubsteps,
                          int compute_mass_balance);

    void mfix_pc_inflow (int const lev,
                         int const is_dem,
                         int const is_pic,
                         amrex::Real dt,
                         const int advect_enthalpy,
                         amrex::EBFArrayBoxFactory* factory);

    void MFIX_PC_InitCollisionParams ();

    void UpdateCost (amrex::MultiFab* cost_mf,
                     const MFIXParIter& particle_tile_iterator,
                     const std::string& cost_type,
                     const amrex::Real& wall_time) const;

    void output (int estatus,
                 int finish,
                 int nstep,
                 amrex::Real dt,
                 amrex::Real time);

    const MFIXFluidPhase& get_fluid () const { return fluid; }

    const MFIXSolidsPhase& get_solids () const { return solids; }

    const MFIXDEM& get_dem () const { return m_dem; }

    const MFIXPIC& get_pic () const {return m_pic; }

    const MFIXReactions& get_reactions () const { return reactions; }

    amrex::Real get_mass_prod ( int const idx ) { return m_p_mass_prod[idx]; }
    amrex::Real get_mass_accum ( int const idx ) { return m_p_mass_accum[idx]; }
    amrex::Real get_mass_inflow ( int const idx ) { return m_p_mass_inflow[idx]; }
    amrex::Real get_mass_outflow ( int const idx ) { return m_p_mass_outflow[idx]; }


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
                                 amrex::Vector< amrex::MultiFab*> & cost,
                                 std::string& knapsack_weight_type);

    const runtimeRealData m_runtimeRealData;

    void setSortingBinSizes (amrex::IntVect sorting_bin) { m_sorting_bin = sorting_bin; }

    amrex::IntVect const& getSortingBinSizes () const noexcept
    { return m_sorting_bin; }

    bool sortNow (int const a_nstep) const noexcept
    { return ((sort_int > 0) && (a_nstep % sort_int == 0)) ? true : false; }

    void partitionParticleGrids(int lev,
                                const amrex::BoxArray& fba,
                                const amrex::DistributionMapping& fdmap,
                                amrex::Real imbalance_toler,
                                amrex::Real partition_factor);

    void setMaxGridSize(const amrex::IntVect& rhs)
    { max_grid_size = rhs;}

    void setMaxGridSize(int size_x, int size_y, int size_z)
    { max_grid_size = amrex::IntVect(AMREX_D_DECL(size_x, size_y, size_z));}

    void setSortInt(int sort_particle_int)
    { sort_int = sort_particle_int; }

    void setReduceGhostParticles(bool reduce)
    { reduceGhostParticles = reduce; }

    void setGreedyRegrid(int dir, bool partition_3d, int min_size=2)
    { greedy_dir = dir; greedy_3d = partition_3d; greedy_min_grid_size=min_size;}

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

    void countParticle(int lev, amrex::iMultiFab& np_mf);

    // verify total # particles hasn't been changed, for debugging
    void verifyParticleCount();

    void printGhostParticleCount();

    void ReadParameters();

    void setParticleFluidGridMap(const amrex::Vector<int>& boxmap);

    void ComputeMassAccum ( int const a_offset );

    void ComputeMassProduction (int const a_lev, amrex::Real const a_dt);

    void ComputeMassOutflow (int const a_lev);

    void ResetMassBalance (int const a_n);

    int UseConstraint ()  noexcept { return m_use_constraint.sum(); }

    void ImposeMean (int const a_lev);

  private:
    int nlev;

    int sort_int = -1; // # steps to sort particles

#ifdef AMREX_USE_GPU
    bool reduceGhostParticles = true;
#else
    bool reduceGhostParticles = false;
#endif

    amrex::Real newton_abstol = 1.e-6;
    amrex::Real newton_reltol = 1.e-6;
    int newton_maxiter = 100;

    amrex::IntVect m_use_constraint = amrex::IntVect::TheZeroVector();
    amrex::RealVect m_constraint = amrex::RealVect::TheZeroVector();

    amrex::IntVect max_grid_size{-1};

    amrex::IntVect m_sorting_bin = amrex::IntVect::TheZeroVector();

    int  greedy_dir           = 0;
    bool greedy_3d            = false;
    int  greedy_min_grid_size = 2;

    long m_total_numparticle = 0;

    amrex::Vector<amrex::Real> m_p_mass_accum;
    amrex::Vector<amrex::Real> m_p_mass_inflow;
    amrex::Vector<amrex::Real> m_p_mass_outflow;
    amrex::Vector<amrex::Real> m_p_mass_prod;

    amrex::Vector<int> m_pboxid_to_fboxid;

    static void ReadStaticParameters ();

    // Local containers for: maxvel - the maximum (absolute) velocity
    //                       maxpfor - the max (abs) particle-particle force
    //                       maxwfor - the max (abs) particle-wall force
    amrex::RealVect loc_maxvel;
    amrex::RealVect loc_maxpfor;
    amrex::RealVect loc_maxwfor;

    void usr0_des ();
    void usr1_des ();

    void usr2_des (int np, ParticleTileType& particles);

    void usr3_des (int np,void*&);

    MFIXInitialConditions& m_initial_conditions;

    MFIXBoundaryConditions& m_boundary_conditions;

    MFIXFluidPhase& fluid;

    MFIXSolidsPhase& solids;

    MFIXDEM& m_dem;

    MFIXPIC& m_pic;

    MFIXReactions& reactions;
};

#endif
