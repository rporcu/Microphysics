#ifndef MFIX_LEVEL_H_
#define MFIX_LEVEL_H_

#include <iostream>
#include <memory>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_NodalProjector.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MacProjector.H>

#include <mfix_diffusion_op.H>
#include <mfix_pc.H>
#include <mfix_bc_list.H>
#include <mfix_eb_if.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_leveldata.H>
#include <mfix_des_heterogeneous_rates_K.H>


// Drag type enumeration
struct DragType { enum {Invalid=-1, WenYu, Gidaspow, BVK2, UserDrag}; };

// Convection type enumeration
struct ConvectionType { enum {Invalid=-1, RanzMarshall, Gunn}; };

// Deposition scheme enumeration
enum DepositionScheme {trilinear, square_dpvm, true_dpvm, centroid};

// Reaction rates type enumeration
struct ReactionRatesType { enum {RRatesUser}; };

class mfix : public AmrCore
{
  public:
    // Constructor
    mfix ();

    // Destructor
    ~mfix ();

    // Face-based coefficients b in MAC projection and implicit diffusion solve
    Vector< Array<MultiFab*,3> > bcoeff;

    void InitParams ();

    bool IsSteadyState ()
    { return (m_steady_state > 0); }

    void Init (Real time);

    void InitLevelData (Real time);

    void init_advection ();

    void mfix_init_solvers ();

    void mfix_setup_solvers ();

    void PostInit (Real& dt, Real time, int restart_flag, Real stop_time);

    void ResizeArrays ();

    void InitialRedistribution (amrex::Real l_time);

    void WriteCheckPointFile (std::string & check_file_name,
                              int nstep = 0, Real dt = 0.0, Real time = 0.0);

    void WritePlotFile (std::string & plot_file_name,
                        int nstep = 0, Real time = 0.0);

    //! Save variables that don't change to plot file. The idea is that they can
    //! be saved _once_ per simulations and not repeatedly every Nth time step.
    void WriteStaticPlotFile (const std::string & plot_file_name) const;

    void Restart (std::string& restart_chkfile, int* nstep, Real* dt,
                  Real* time, IntVect& Nrep);

    void WriteParticleAscii (std::string& par_ascii_file_name, int nstep = 0) const;

    void WriteAverageRegions (std::string& avg_file, int nstep, Real time = 0.) const;

    void ComputeAverageFluidVars (const int lev,
                                  const amrex::Real time,
                                  const std::string& basename) const;

    void Regrid ();

    void Evolve (int nstep, Real & dt, Real & prev_dt, Real time, Real stop_time);

    void mfix_usr0 () const;
    void mfix_usr1 (Real time) const;
    void mfix_usr2 () const;
    void mfix_usr3 () const;

    void output (int estatus, int finish, int nstep, Real dt, Real time)
    { pc->output(estatus, finish, nstep, dt, time); };

    void make_mg_bc (int mg_bc[]);

    void make_eb_geometry ();

    void make_eb_factories ();

    void fill_eb_levelsets ();

    void intersect_ls_walls ();

    //! Construct EB levels from Geometry shop. This builds each EB level and
    //! saves a pointer to each level into `eb_levels`.
    template<class F> void build_eb_levels (EB2::GeometryShop<F> gshop)
    {
      /************************************************************************
       *                                                                      *
       * Build EB levels                                                      *
       *                                                                      *
       ***********************************************************************/
      EB2::Build(gshop, geom[nlev-1], nlev-1, 100);
      const EB2::IndexSpace& ebis = EB2::IndexSpace::top();
      for (int lev = 0; lev < nlev; lev ++)
      {
          eb_levels[lev] = &(ebis.getLevel(geom[lev]));
          particle_eb_levels[lev] = eb_levels[lev];
      }

      if (nlev == 1)
      {
          // We could skip the rebuilding for ref ratio of 1, but that makes
          // BENCH05-Size0008_wide to fail on 4 processes.
//          if (levelset_refinement == 1) {
//              eb_levels[1] = eb_levels[0];
//              particle_eb_levels[1] = eb_levels[1];
//          } else {
              Geometry geom_ls = amrex::refine(geom[0],levelset_refinement);
              EB2::Build(gshop, geom_ls, 0, 100);
              eb_levels[1] = &(EB2::IndexSpace::top().getLevel(geom_ls));
              particle_eb_levels[1] = eb_levels[1];
//          }
      }
    }

    const Geometry get_geom_ref (int lev)
    { return geom[lev]; }

    std::shared_ptr<UnionListIF<EB2::PlaneIF>> get_walls (bool & has_walls);

    std::shared_ptr<UnionListIF<EB2::PlaneIF>> get_real_walls (bool & has_real_walls);

    std::shared_ptr<EB2::TranslationIF<EB2::PolynomialIF>>
      get_poly (int max_order, std::string field_prefix);

    std::shared_ptr<EB2::IntersectionIF<EB2::PlaneIF,EB2::PlaneIF,EB2::PlaneIF>>
      make_wall (int dir, // direction (long edge) of wall
                 Real position, Real height, Real width);

    void WriteMyEBSurface ();

    void InitIOChkData ();

    void ResetIOChkData ();

    void InitIOPltData ();

    // If > 0, dump plotfiles every "plot_int" steps
    static int plot_int;

    // If > 0, dump plotfiles at as close as possible to the designated period
    // *without* changing dt
    static Real plot_per_approx;

    // If > 0, dump plotfiles at exactcly the designated period by changing dt
    static Real plot_per_exact;

    static inline std::string get_load_balance_type ()
    { return load_balance_type; }

    // flag enabling level-set restart (i.e. prevent make_eb_* from rebuilding
    // the level-set data).
    bool levelset_restart = false;

    void mfix_compute_vort ();

    void mfix_compute_dt (int nstep, Real time, Real stop_time, Real & dt, Real & prev_dt);

    ///
    /// Parameters describing the level-set grid:
    ///
    int levelset_refinement = 1;    // refinement (wrt particle grid) of the level-set's grid
    int levelset_eb_refinement = 1; // refinement of the EB facets used to generate the level-set
    int levelset_pad = 2;           // padding of the level-set grid
    int levelset_eb_pad = 2;        // padding of the EBIS used to generate the level-set (useful for
                                    // ensuring that neighbor grids are considered when filling LS)

    bool contains_ebs = false;

    void mfix_add_txfr_explicit (Real dt);

    void mfix_calc_txfr_particle (Real time,
                                  Vector< MultiFab* > const& vel_g_in,
                                  Vector< MultiFab* > const& gp_in,
                                  Vector< MultiFab* > const& T_g_in);

    void mfix_add_txfr_implicit (amrex::Real dt,
                                 amrex::Vector<amrex::MultiFab*      > const& vel_in,
                                 amrex::Vector<amrex::MultiFab*      > const& h_g_in,
                                 amrex::Vector<amrex::MultiFab*      > const& T_g_in,
                                 amrex::Vector<amrex::MultiFab const*> const& txfr_in,
                                 amrex::Vector<amrex::MultiFab const*> const& rho_in,
                                 amrex::Vector<amrex::MultiFab const*> const& ep_g_in,
                                 amrex::Vector<amrex::MultiFab const*> const& cp_g_in);

    void mfix_correct_small_cells (Vector< MultiFab*     > const& vel_in,
                                   Vector<MultiFab const*> const& ep_u_mac,
                                   Vector<MultiFab const*> const& ep_v_mac,
                                   Vector<MultiFab const*> const& ep_w_mac);

    void mfix_set_velocity_bcs (Real time,
                                Vector< MultiFab* > const& vel,
                                int extrap_dir_bcs) const;

    void mfix_set_epg_bcs (const Vector< MultiFab* >& ep_g, const int dir_bc) const;

    void mfix_set_enthalpy_bcs (Real time,
                                Vector< MultiFab* > const& h_g);

    void mfix_set_temperature_bcs (Real time,
                                   Vector< MultiFab* > const& T_g);

    void mfix_set_eb_temperature_bcs (Vector< MultiFab* > const& eb_T_g_in,
                                      Vector< MultiFab* > const& eb_k_g_in);

    void mfix_set_density_bcs (Real time,
                               Vector< MultiFab* > const& ro_g);

    void mfix_set_tracer_bcs (Real time,
                              Vector< MultiFab* > const& trac);

    void mfix_set_scalar_bcs (Real time,
                              Vector< MultiFab* > const& mu_g,
                              Vector< MultiFab* > const& cp_g,
                              Vector< MultiFab* > const& k_g,
                              Vector< MultiFab* > const& MW_g);

    void mfix_set_species_bcs (Real time,
                               Vector< MultiFab* > const& X_gk,
                               Vector< MultiFab* > const& D_gk,
                               Vector< MultiFab* > const& cp_gk,
                               Vector< MultiFab* > const& h_gk);

    void set_velocity_bcs (Real time,
                           const int lev,
                           FArrayBox& vel_fab,
                           const Box& domain,
                           const int* extrap_dir_bcs) const;

    void set_epg_bcs (const int lev,
                      FArrayBox& eps_fab,
                      const Box& domain,
                      const int* extrap_dir_bcs) const;

    void set_enthalpy_bcs (Real time,
                           const int lev,
                           FArrayBox& h_g_fab,
                           const Box& domain);

    void set_temperature_bcs (Real time,
                              const int lev,
                              FArrayBox& T_g_fab,
                              const Box& domain);

    void set_eb_temperature_bcs (const Box& sbx,
                                 const Box& bx,
                                 const Box& domain,
                                 const Real dx,
                                 const Real dy,
                                 const Real dz,
                                 const GpuArray<Real, AMREX_SPACEDIM>& plo,
                                 FArrayBox& eb_T_g_fab,
                                 FArrayBox& eb_k_g_fab,
                                 FArrayBox& k_g_fab,
                                 const EBCellFlagFab& flags_fab);

    void set_density_bcs (Real time,
                          const int lev,
                          FArrayBox& scal_fab,
                          const Box& domain);

    void set_tracer_bcs (Real time,
                         const int lev,
                         FArrayBox& trac_fab,
                         const Box& domain);

    void set_viscosity_bcs (Real time,
                            const int lev,
                            FArrayBox& scal_fab,
                            const Box& domain);

    void set_specific_heat_bcs (Real time,
                                const int lev,
                                FArrayBox& scal_fab,
                                const Box& domain);

    void set_thermal_conductivity_bcs (Real time,
                                       const int lev,
                                       FArrayBox& scal_fab,
                                       const Box& domain);

    void set_molecular_weight_bcs (Real time,
                                   const int lev,
                                   FArrayBox& scal_fab,
                                   const Box& domain);

    void set_mass_fractions_g_bcs (Real time,
                                   const int lev,
                                   FArrayBox& X_gk_fab,
                                   const Box& domain);

    void set_species_diffusivities_g_bcs (Real time,
                                          const int lev,
                                          FArrayBox& D_gk_fab,
                                          const Box& domain);

    void set_species_specific_heat_g_bcs (Real time,
                                          const int lev,
                                          FArrayBox& cp_gk_fab,
                                          const Box& domain);

    void set_species_enthalpy_g_bcs (Real time,
                                     const int lev,
                                     FArrayBox& h_gk_fab,
                                     const Box& domain);

    void set_vec_bcs (const int lev,
                      FArrayBox& vec_fab,
                      const Box& domain) const;

    void set_gradp_bcs (const Box& bx,
                        const int lev,
                        FArrayBox& gp_fab,
                        Box& domain);

    void set_ls_inflow (const int lev,
                        FArrayBox& ls_phi_fab,
                        const Box& domain,
                        const int* levelset_nghost,
                        const int& nref,
                        const Real* dx);

    void set_bc0 (const Box& sbx,
                  MFIter* mfi,
                  const int lev,
                  const Box& domain);

    void set_temperature_bc0 (const Box& sbx,
                              MFIter* mfi,
                              const int lev,
                              const Box& domain);

    void set_species_bc0 (const Box& sbx,
                          MFIter* mfi,
                          const int lev,
                          const Box& domain);

    void set_p0 (const Box& bx,
                 MFIter* mfi,
                 const int lev,
                 const Box& domain);

    void set_gp0 (const int lev,
                  const Box& domain);

    static int m_drag_type;

    static int m_convection_type;

    static DepositionScheme m_deposition_scheme;

    static int m_reaction_rates_type;

    static amrex::Real m_deposition_diffusion_coeff;
    static amrex::Real m_deposition_scale_factor;
    static amrex::Real m_max_solids_volume_fraction;

    template <typename F1>
    void mfix_calc_transfer_coeffs (Vector< MultiFab* > const& ep_g_in,
                                    Vector< MultiFab* > const& ro_g_in,
                                    Vector< MultiFab* > const& vel_g_in,
                                    Vector< MultiFab* > const& mu_g_in,
                                    Vector< MultiFab* > const& cp_g_in,
                                    Vector< MultiFab* > const& k_g_in,
                                    F1 DragFunc);

    template <typename F1, typename F2>
    void mfix_calc_transfer_coeffs (Vector< MultiFab* > const& ep_g_in,
                                    Vector< MultiFab* > const& ro_g_in,
                                    Vector< MultiFab* > const& vel_g_in,
                                    Vector< MultiFab* > const& mu_g_in,
                                    Vector< MultiFab* > const& cp_g_in,
                                    Vector< MultiFab* > const& k_g_in,
                                    F1 DragFunc,
                                    F2 ConvectionCoeff);

    void mfix_calc_volume_fraction (Real & sum_vol);

    void mfix_deposition_bcs (int lev, amrex::MultiFab & filled_mf);

    void mfix_redistribute_deposition (int lev,
                                       amrex::MultiFab & mf_eps,
                                       amrex::MultiFab & mf_to_redistribute,
                                       const amrex::MultiFab * volfrac,
                                       const amrex::FabArray<EBCellFlagFab>* flags,
                                       amrex::Real max_eps);

    Real volWgtSum (int lev, const MultiFab & mf, int comp, bool local=false) const;

    Real volEpsWgtSum (int lev, const MultiFab & mf, int comp, bool local=false) const;

    Real volWgtSumBox (int lev, const MultiFab & mf, int comp, const Box a_bx, bool local=false) const;

    void ReportGridStats () const;

    void mfix_set_bc_type (int lev, int nghost_bc);

    void set_bcrec_lo (const int lev, const int dir, const int l_type);
    void set_bcrec_hi (const int lev, const int dir, const int l_type);

    ////////////////////////////////////////////////////////////////////////////
    //
    // Convection
    //
    ////////////////////////////////////////////////////////////////////////////

    void set_MAC_velocity_bcs (int lev,
                               amrex::Vector< amrex::MultiFab const*> const& mac_rhs,
                               amrex::Vector< amrex::MultiFab* > const& u,
                               amrex::Vector< amrex::MultiFab* > const& v,
                               amrex::Vector< amrex::MultiFab* > const& w,
                               amrex::Real time);

  void compute_MAC_projected_velocities (amrex::Real time, const amrex::Real l_dt,
                                         amrex::Vector< amrex::MultiFab const*> const& vel_in,
                                         amrex::Vector< amrex::MultiFab*      > const& ep_u_mac,
                                         amrex::Vector< amrex::MultiFab*      > const& ep_v_mac,
                                         amrex::Vector< amrex::MultiFab*      > const& ep_w_mac,
                                         amrex::Vector< amrex::MultiFab const*> const& ep_g_in,
                                         amrex::Vector< amrex::MultiFab const*> const& ro_g_in,
                                         amrex::Vector< amrex::MultiFab      *> const& vel_forces,
                                         amrex::Vector< amrex::MultiFab const*> const& rhs_mac);


    ////////////////////////////////////////////////////////////////////////////
    //
    // Diffusion
    //
    ////////////////////////////////////////////////////////////////////////////

    std::shared_ptr<DiffusionOp> diffusion_op;

    ////////////////////////////////////////////////////////////////////////////
    //
    // Nodal Projection
    //
    ////////////////////////////////////////////////////////////////////////////

    std::shared_ptr<NodalProjector> nodal_projector;



  ////////////////////////////////////////////////////////////////////////////
  //
  // PIC solids
  //
  ////////////////////////////////////////////////////////////////////////////

  void EvolveParcels (amrex::Real dt,
                      amrex::Real time,
                      amrex::RealVect& gravity,
                      const int ls_refinement,
                      amrex::Vector< amrex::MultiFab* >& cost,
                      std::string& knapsack_weight_type,
                      const int advect_enthalpy);


  void MFIX_CalcAvgSolidsVel (Vector< MultiFab* >& avg_prop_in);

  void MFIX_CalcSolidsStress (amrex::Vector< amrex::MultiFab* >& ep_s_in,
                              amrex::Vector< amrex::MultiFab* >& avg_prop_in,
                              amrex::Vector< amrex::MultiFab* >& cost,
                              std::string& knapsack_weight_type);


  protected:
    //! Tagging cells for refinement
    virtual void ErrorEst (int lev,
                           TagBoxArray & tags,
                           Real time,
                           int ngrow) override;

    //! Make a new level using provided BoxArray and DistributionMapping and
    //! fill with interpolated coarse level data. Called by AmrCore::regrid.
    virtual void MakeNewLevelFromCoarse (int /*lev*/, Real /*time*/,
                                         const BoxArray & /*ba*/,
                                         const DistributionMapping & /*dm*/) override
    { amrex::Abort("MakeNewLevelFromCoarse: To be implemented"); }

    //! Remake an existing level using provided BoxArray and
    //! DistributionMapping and fill with existing fine and coarse data. Called
    //! by AmrCore::regrid.
    virtual void RemakeLevel (int /*lev*/, Real /*time*/,
                              const BoxArray & /*ba*/,
                              const DistributionMapping & /*dm*/) override
    { amrex::Abort("RemakeLevel: To be implemented"); }

    //! Delete level data.  Called by AmrCore::regrid.
    virtual void ClearLevel (int /*lev*/) override
    { amrex::Abort("ClearLevel: To be implemented"); }

    void mfix_init_fluid (int is_restarting, Real dt, Real stop_time);

    void mfix_set_bc0 ();

    void mfix_set_p0 ();

    void mfix_set_gp0 ();

    static RealVect gp0;
    static RealVect gravity;

    // Scalar diffusive coefficient
    Vector<Real> mu_s;

    BCList bc_list;

    Gpu::DeviceVector<Real> m_bc_u_g;
    Gpu::DeviceVector<Real> m_bc_v_g;
    Gpu::DeviceVector<Real> m_bc_w_g;
    Gpu::DeviceVector<Real> m_bc_t_g;
    Gpu::DeviceVector<Real> m_bc_ep_g;
    Gpu::DeviceVector<Real> m_bc_p_g;
    Vector< Gpu::DeviceVector<Real> > m_bc_X_gk;
    Gpu::DeviceVector<Real*> m_bc_X_gk_ptr;
    Vector<Real> m_h_bc_u_g;
    Vector<Real> m_h_bc_v_g;
    Vector<Real> m_h_bc_w_g;
    Vector<Real> m_h_bc_t_g;
    Vector<Real> m_h_bc_ep_g;
    Vector<Real> m_h_bc_p_g;
    Vector<Vector<Real> > m_h_bc_X_gk;

    std::map< std::string, Gpu::DeviceVector<int> > m_vel_g_bc_types;
    std::map< std::string, Gpu::DeviceVector<int> > m_ro_g_bc_types;
    std::map< std::string, Gpu::DeviceVector<int> > m_T_g_bc_types;
    std::map< std::string, Gpu::DeviceVector<int> > m_trac_g_bc_types;
    std::map< std::string, Gpu::DeviceVector<int> > m_X_gk_bc_types;

    //! The particles see an inflow face as a solid wall; thus instead of
    //! creating a separate EBfactory for particles, we simply modify the
    //! level-set near inflow to see it as a wall
    void mfix_set_ls_near_inflow ();

    void mfix_project_velocity ();

    void mfix_initial_iterations (Real dt, Real stop_time);

public:

    void mfix_apply_nodal_projection (Vector< amrex::MultiFab* >& a_depdt,
                                      amrex::Real a_time,
                                      amrex::Real a_dt,
                                      amrex::Real a_prev_dt,
                                      bool proj_2,
                                      amrex::Vector<amrex::MultiFab*      > const& vel_g_old_in,
                                      amrex::Vector<amrex::MultiFab*      > const& vel_g_in,
                                      amrex::Vector<amrex::MultiFab*      > const& p_g_in,
                                      amrex::Vector<amrex::MultiFab*      > const& gp_in,
                                      amrex::Vector<amrex::MultiFab*      > const& ep_g_in,
                                      amrex::Vector<amrex::MultiFab*      > const& txfr_in,
                                      amrex::Vector<amrex::MultiFab const*> const& density);

    void compute_tra_forces (amrex::Vector<amrex::MultiFab      *> const& tra_forces,
                             amrex::Vector<amrex::MultiFab const*> const& density);

    void compute_vel_forces (amrex::Vector<amrex::MultiFab*      > const& vel_forces,
                             amrex::Vector<amrex::MultiFab const*> const& velocity,
                             amrex::Vector<amrex::MultiFab const*> const& density,
                             amrex::Vector<amrex::MultiFab const*> const& txfr_in,
                             bool include_pressure_gradient = true,
                             bool include_drag_force = false);

    void compute_vel_forces_on_level ( int lev,
                                       amrex::MultiFab& vel_forces,
                                       const amrex::MultiFab& velocity,
                                       const amrex::MultiFab& density,
                                       const amrex::MultiFab& txfr_in,
                                       bool include_pressure_gradient = true,
                                       bool include_drag_force = false);

    void compute_laps (const bool update_lapT,
                       const bool update_LapTrac,
                       const bool update_LapX,
                       amrex::Vector< amrex::MultiFab*      >& lapT,
                       amrex::Vector< amrex::MultiFab*      >& lapS,
                       amrex::Vector< amrex::MultiFab*      >& lapX,
                       amrex::Vector< amrex::MultiFab*      > const& T_g_in,
                       amrex::Vector< amrex::MultiFab*      > const& trac_in,
                       amrex::Vector< amrex::MultiFab*      > const& X_gk_in,
                       amrex::Vector< amrex::MultiFab const*> const& ep_g_in,
                       amrex::Vector< amrex::MultiFab const*> const& ro_g_in);



    void apply_MAC_projection (const bool update_laplacians,
                               Vector< MultiFab* > const& lap_T,
                               Vector< MultiFab* > const& lap_X,
                               Vector< MultiFab* > const& ep_u_mac,
                               Vector< MultiFab* > const& ep_v_mac,
                               Vector< MultiFab* > const& ep_w_mac,
                               Vector< MultiFab* > const& ep_g_in,
                               Vector< MultiFab* > const& ro_g_in,
                               Vector< MultiFab* > const& MW_g_in,
                               Vector< MultiFab* > const& T_g_in,
                               Vector< MultiFab* > const& cp_g_in,
                               Vector< MultiFab* > const& k_g_in,
                               Vector< MultiFab* > const& T_g_on_eb_in,
                               Vector< MultiFab* > const& k_g_on_eb_in,
                               Vector< MultiFab* > const& X_gk_in,
                               Vector< MultiFab* > const& D_gk_in,
                               Vector< MultiFab* > const& h_gk_in,
                               Vector< MultiFab* > const& txfr_in,
                               Vector< MultiFab* > const& ro_gk_txfr_in,
                               Real time);

  protected:
    void FillPatchVel (int lev, Real time, MultiFab & mf, int icomp, int ncomp,
                       const Vector<BCRec> & bcr);

    void GetDataVel (int lev, Real time, Vector<MultiFab *> & data,
                     Vector<Real> & datatime);

    void FillPatchScalar (int lev, Real time, MultiFab & mf, int icomp, int ncomp,
                          const Vector<BCRec> & bcr);

    void FillPatchSpecies (int lev, Real time, MultiFab & mf, int icomp, int ncomp,
                           const Vector<BCRec> & bcr);

    void GetDataScalar (int lev, Real time, Vector<MultiFab *> & data, int icomp,
                        Vector<Real> & datatime);

    void GetDataSpecies (int lev, Real time, Vector<MultiFab *> & data, int icomp,
                         Vector<Real> & datatime);


    void fillpatch_force (amrex::Real time, amrex::Vector<amrex::MultiFab*> const& force, int ng);

    void mfix_calc_transfer_coeffs (Vector< MultiFab* > const& ep_g_in,
                                    Vector< MultiFab* > const& ro_g_in,
                                    Vector< MultiFab* > const& vel_g_in,
                                    Vector< MultiFab* > const& mu_g_in,
                                    Vector< MultiFab* > const& cp_g_in,
                                    Vector< MultiFab* > const& k_g_in);

    void mfix_calc_txfr_fluid (Vector<MultiFab*> const& txfr,
                               Vector< MultiFab* > const& ep_g_in,
                               Vector< MultiFab* > const& ro_g_in,
                               Vector< MultiFab* > const& vel_g_in,
                               Vector< MultiFab* > const& mu_g_in,
                               Vector< MultiFab* > const& cp_g_in,
                               Vector< MultiFab* > const& k_g_in,
                               Real time);

  protected:
    // Vector of allowed chemical reactions
    amrex::Vector< REACTIONS::ChemicalReaction* > m_chemical_reactions;

  public:
    // Initialization: read input parameters and set up reactions
    void initialize_chem_reactions ();

    void mfix_calc_chem_txfr (const Real time,
                              const Vector< MultiFab* >& ep_g_in,
                              const Vector< MultiFab* >& ro_g_in,
                              const Vector< MultiFab* >& X_gk_in);

    template <typename F1>
    void mfix_calc_chem_txfr (const Real time,
                              const Vector< MultiFab* >& ep_g_in,
                              const Vector< MultiFab* >& ro_g_in,
                              const Vector< MultiFab* >& X_gk_in,
                              F1 RRatesFunc);

  protected:
    void AllocateArrays (int lev);

    void RegridArrays (int lev);

    void RegridLevelSetArray (int a_lev);

    void avgDown (int crse_lev, const MultiFab & S_fine, MultiFab & S_crse);

  private:
    void make_eb_general ();
    void make_eb_box ();
    void make_eb_cylinder ();
    void make_eb_hopper ();
    void make_eb_cyclone ();
    void make_eb_air_reactor ();
    void make_eb_proto_clr ();
    void make_eb_hourglass ();
    void make_eb_clr ();
    void make_eb_clr_riser ();
    void make_eb_regular ();

#ifdef MFIX_GEOMETRY_CSG
  void make_eb_csg(std::string);
#endif

    void MakeBCArrays (int nghost_bc);

    void MakeNewLevelFromScratch (int lev, Real time, const BoxArray & new_grids,
                                  const DistributionMapping & new_dmap) override;

    void ReMakeNewLevelFromScratch (int lev, const BoxArray & new_grids,
                                    const DistributionMapping & new_dmap);

    void WriteCheckHeader (const std::string & name, int nstep, Real dt, Real time) const;

    void WriteJobInfo (const std::string & dir ) const;

    static void GotoNextLine (std::istream & is);

    void EvolveFluid (int nstep, Real & dt, Real & time, Real & prev_dt, Real stop_time, Real & drag_timing);

    void check_for_nans (int lev);


    void fillpatch_all (Vector< MultiFab* > const& vel_in,
                        Vector< MultiFab* > const& ro_g_in,
                        Vector< MultiFab* > const& h_g_in,
                        Vector< MultiFab* > const& trac_in,
                        Vector< MultiFab* > const& X_gk_in,
                        Real time);


    void mfix_compute_convective_term (Vector< MultiFab*      >& conv_u,
                                       Vector< MultiFab*      >& conv_s,
                                       Vector< MultiFab*      >& conv_X,
                                       Vector< MultiFab*      > const& vel_forces,
                                       Vector< MultiFab*      > const& tra_forces,
                                       Vector< MultiFab const*> const& vel_in,
                                       Vector< MultiFab const*> const& ep_g_in,
                                       Vector< MultiFab const*> const& ro_g_in,
                                       Vector< MultiFab const*> const& h_g_in,
                                       Vector< MultiFab const*> const& trac_in,
                                       Vector< MultiFab const*> const& X_gk_in,
                                       Vector< MultiFab const*> const& txfr_in,
                                       Vector< MultiFab*      > const& u_mac,
                                       Vector< MultiFab*      > const& v_mac,
                                       Vector< MultiFab*      > const& w_mac,
                                       Vector< MultiFab const*> const& rhs_mac,
                                       Vector< MultiFab      *> const& divtau_old,
                                       Real l_dt, Real time);




    int steady_state_reached (Real dt, int iter);

  public:
    void mfix_apply_predictor (Vector< MultiFab* >& conv_u_old,
                               Vector< MultiFab* >& conv_s_old,
                               Vector< MultiFab* >& conv_X_old,
                               Vector< MultiFab* >& ro_RHS_old,
                               Vector< MultiFab* >& ro_RHS,
                               Vector< MultiFab* >& lap_trac_old,
                               Vector< MultiFab* >& lap_trac,
                               Vector< MultiFab* >& enthalpy_RHS_old,
                               Vector< MultiFab* >& enthalpy_RHS,
                               Vector< MultiFab* >& lap_T_old,
                               Vector< MultiFab* >& lap_T,
                               Vector< MultiFab* >& species_RHS_old,
                               Vector< MultiFab* >& species_RHS,
                               Vector< MultiFab* >& lap_X_old,
                               Vector< MultiFab* >& lap_X,
                               Real time,
                               Real dt,
                               Real prev_dt,
                               bool proj_2);

    void mfix_density_rhs(Vector< MultiFab*      > const& rhs,
                          Vector< MultiFab const*> const& ro_gk_txfr);

    void mfix_enthalpy_rhs(Vector< MultiFab*      > const& rhs,
                           Vector< MultiFab const*> const& ep_g,
                           Vector< MultiFab const*> const& ro_g,
                           Vector< MultiFab const*> const& X_gk,
                           Vector< MultiFab const*> const& D_gk,
                           Vector< MultiFab const*> const& h_gk);

    void mfix_scalar_rhs(const bool explicit_diffusion,
                         Vector< MultiFab* > const& lap_trac,
                         Vector< MultiFab* > const& trac,
                         Vector< MultiFab* > const& ep_g,
                         Vector< MultiFab* > const& ro_g,
                         const Vector<Real>& mu_s);

    void mfix_species_X_rhs(Vector< MultiFab*      > const& rhs,
                            Vector< MultiFab const*> const& ro_gk_txfr);

    void mfix_apply_corrector (Vector< MultiFab* >& conv_u_old,
                               Vector< MultiFab* >& conv_s_old,
                               Vector< MultiFab* >& conv_X_old,
                               Vector< MultiFab* >& ro_RHS_old,
                               Vector< MultiFab* >& ro_RHS,
                               Vector< MultiFab* >& lap_trac_old,
                               Vector< MultiFab* >& lap_trac,
                               Vector< MultiFab* >& enthalpy_RHS_old,
                               Vector< MultiFab* >& enthalpy_RHS,
                               Vector< MultiFab* >& lap_T_old,
                               Vector< MultiFab* >& lap_T,
                               Vector< MultiFab* >& species_RHS_old,
                               Vector< MultiFab* >& species_RHS,
                               Vector< MultiFab* >& lap_X_old,
                               Vector< MultiFab* >& lap_X,
                               Real time,
                               Real dt,
                               Real prev_dt,
                               bool proj_2);

  public:
    void mfix_normalize_fluid_species(const Vector< MultiFab* >& X_gk);

    void mfix_update_fluid_and_species(const Vector< MultiFab* >& cp_gk,
                                        const Vector< MultiFab* >& h_gk,
                                        const Vector< MultiFab* >& MW_g,
                                        const Vector< MultiFab* >& cp_g,
                                        const Vector< MultiFab* >& h_g,
                                        const Vector< MultiFab* >& T_g,
                                        const Vector< MultiFab* >& X_gk);

    void mfix_open_system_rhs (Vector< MultiFab*       > const& rhs,
                               Vector< MultiFab*       > const& lap_T_star,
                               Vector< MultiFab*       > const& lap_X_star,
                               Vector< MultiFab const* > const& ep_g,
                               Vector< MultiFab const* > const& ro_g,
                               Vector< MultiFab const* > const& MW_g,
                               Vector< MultiFab const* > const& T_g,
                               Vector< MultiFab const* > const& cp_g,
                               Vector< MultiFab const* > const& X_gk,
                               Vector< MultiFab const* > const& D_gk,
                               Vector< MultiFab const* > const& h_gk,
                               Vector< MultiFab const* > const& txfr,
                               Vector< MultiFab const* > const& ro_gk_txfr);

    void compute_convective_term (Box const& bx, int lev, const Real l_dt, MFIter const& mfi,
                                  Array4<Real> const& dvdt, // velocity
                                  Array4<Real> const& dsdt, // density, enthalpy, tracer
                                  Array4<Real> const& dXdt, // tracer
                                  const bool l_advect_density,
                                  const bool l_advect_enthalpy,
                                  const bool l_advect_tracer,
                                  const bool l_advect_species,
                                  Array4<Real const> const& vel,
                                  Array4<Real const> const& ep_g,
                                  Array4<Real const> const& rho,
                                  Array4<Real const> const& hg,
                                  Array4<Real const> const& tra,
                                  Array4<Real const> const& Xgk,
                                  const int l_nspecies,
                                  Array4<Real const> const& divu,
                                  Array4<Real const> const& ep_umac,
                                  Array4<Real const> const& ep_vmac,
                                  Array4<Real const> const& ep_wmac,
                                  Array4<Real const> const& fvel,
                                  Array4<Real const> const& ftra);

  private:
    void mfix_print_max_vel (int lev);

    void mfix_print_max_gp (int lev);

    void mfix_print_min_epg ();

    // Make a level 0 grids covering the whole domain.  It does NOT install the new grids.
    BoxArray MakeBaseGrids () const;

    void ChopGrids (const Box& domain, BoxArray& ba, int target_size) const;

    static int nlev;
    static int ntrac;

    // Particle container
    MFIXParticleContainer* pc;

    // Debug counter used for tracking number of level-set MultiFabs used in
    // intersection/union
    int ct_ls_mf;

    // Boundary conditions types
    Vector< IArrayBox* > bc_ilo;
    Vector< IArrayBox* > bc_ihi;
    Vector< IArrayBox* > bc_jlo;
    Vector< IArrayBox* > bc_jhi;
    Vector< IArrayBox* > bc_klo;
    Vector< IArrayBox* > bc_khi;

    Vector< std::unique_ptr<LevelData> > m_leveldata;

    Vector< MultiFab* > get_ep_g () noexcept;
    Vector< MultiFab* > get_ro_g () noexcept;
    Vector< MultiFab* > get_ro_g_old () noexcept;
    Vector< MultiFab* > get_MW_g () noexcept;
    Vector< MultiFab* > get_trac () noexcept;
    Vector< MultiFab* > get_trac_old () noexcept;
    Vector< MultiFab* > get_vel_g () noexcept;
    Vector< MultiFab* > get_vel_g_old () noexcept;
    Vector< MultiFab* > get_gp () noexcept;
    Vector< MultiFab* > get_p_g () noexcept;
    Vector< MultiFab* > get_p_g_old () noexcept;
    Vector< MultiFab* > get_mu_g () noexcept;
    Vector< MultiFab* > get_T_g () noexcept;
    Vector< MultiFab* > get_T_g_old () noexcept;
    Vector< MultiFab* > get_cp_g () noexcept;
    Vector< MultiFab* > get_k_g () noexcept;
    Vector< MultiFab* > get_h_g () noexcept;
    Vector< MultiFab* > get_h_g_old () noexcept;
    Vector< MultiFab* > get_T_g_on_eb () noexcept;
    Vector< MultiFab* > get_k_g_on_eb () noexcept;
    Vector< MultiFab* > get_X_gk () noexcept;
    Vector< MultiFab* > get_X_gk_old () noexcept;
    Vector< MultiFab* > get_D_gk () noexcept;
    Vector< MultiFab* > get_cp_gk () noexcept;
    Vector< MultiFab* > get_h_gk () noexcept;
    Vector< MultiFab* > get_txfr () noexcept;
    Vector< MultiFab* > get_ro_gk_txfr () noexcept;
    Vector< MultiFab* > get_diveu () noexcept;
    Vector< MultiFab* > get_mac_phi () noexcept;
    Vector< MultiFab* > get_divtau () noexcept;


    amrex::Vector< amrex::MultiFab const* > get_ep_g_const () const noexcept;

    amrex::Vector< amrex::MultiFab const* > get_ro_g_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_ro_g_old_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_trac_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_trac_old_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_vel_g_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_vel_g_old_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_p_g_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_p_g_old_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_T_g_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_T_g_old_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_h_g_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_h_g_old_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_X_gk_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_X_gk_old_const () const noexcept;

    amrex::Vector< MultiFab const* > get_divtau_const () const noexcept;

    amrex::Vector< amrex::MultiFab const* > get_MW_g_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_mu_g_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_cp_g_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_k_g_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_T_g_on_eb_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_k_g_on_eb_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_D_gk_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_cp_gk_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_h_gk_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_txfr_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_ro_gk_txfr_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_diveu_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_mac_phi_const () const noexcept;

  amrex::Vector<int> const& get_velocity_iconserv () const noexcept { return m_iconserv_velocity; }
  amrex::Vector<int> const& get_density_iconserv () const noexcept { return m_iconserv_density; }
  amrex::Vector<int> const& get_enthalpy_iconserv () const noexcept { return m_iconserv_enthalpy; }
  amrex::Vector<int> const& get_tracer_iconserv () const noexcept { return m_iconserv_tracer; }
  amrex::Vector<int> const& get_species_iconserv () const noexcept { return m_iconserv_species; }
  //
  int const* get_velocity_iconserv_device_ptr () const noexcept {
    return m_iconserv_velocity_d.data(); }
  int const* get_density_iconserv_device_ptr () const noexcept {
    return m_iconserv_density_d.data(); }
  int const* get_enthalpy_iconserv_device_ptr () const noexcept {
    return m_iconserv_enthalpy_d.data(); }
  int const* get_tracer_iconserv_device_ptr () const noexcept {
    return m_iconserv_tracer_d.data(); }
  int const* get_species_iconserv_device_ptr () const noexcept {
    return m_iconserv_species_d.data(); }


  amrex::Vector<amrex::BCRec> const& get_velocity_bcrec () const noexcept { return m_bcrec_velocity; }
  amrex::Vector<amrex::BCRec> const& get_density_bcrec () const noexcept { return m_bcrec_density; }
  amrex::Vector<amrex::BCRec> const& get_enthalpy_bcrec () const noexcept { return m_bcrec_enthalpy; }
  amrex::Vector<amrex::BCRec> const& get_tracer_bcrec () const noexcept { return m_bcrec_tracer; }
  amrex::Vector<amrex::BCRec> const& get_species_bcrec () const noexcept { return m_bcrec_species; }
  amrex::Vector<amrex::BCRec> const& get_force_bcrec () const noexcept { return m_bcrec_force; }
  //
  amrex::BCRec const* get_velocity_bcrec_device_ptr () const noexcept {
    return m_bcrec_velocity_d.data(); }
  amrex::BCRec const* get_density_bcrec_device_ptr () const noexcept {
    return m_bcrec_density_d.data(); }
  amrex::BCRec const* get_enthalpy_bcrec_device_ptr () const noexcept {
    return m_bcrec_enthalpy_d.data(); }
  amrex::BCRec const* get_tracer_bcrec_device_ptr () const noexcept {
    return m_bcrec_tracer_d.data(); }
  amrex::BCRec const* get_species_bcrec_device_ptr () const noexcept {
    return m_bcrec_species_d.data(); }
  amrex::BCRec const* get_force_bcrec_device_ptr () const noexcept {
    return m_bcrec_force_d.data(); }


  // Although we do nothing special for different scalars,
  // we assume they are stored in the following order:
  // 0: density, 1: enthalpy, 2: tracers
  // 3+ntrac: ep_g,  4+ntrac: mu_g
  // 5+ntrac: T_g,   6+ntrac: h_g


    // used if load_balance_type == "KnapSack"
    amrex::Vector< amrex::MultiFab* > particle_cost;
    amrex::Vector< amrex::MultiFab* > fluid_cost;

    // Level-Set Data => used for particle-wall collisions and fluid
    // reconstruction in particle drag calculation. NOTE: this has at least 2
    // levels: 0) fluid reconstruction, 1..N) particle-wall (refined) collisions
    Vector< std::unique_ptr<MultiFab> > level_sets;

    int m_steady_state = 0;

    bool advect_enthalpy        = false;
    bool advect_density         = false;
    bool advect_tracer          = false;
    bool advect_fluid_species   = false;
    bool solve_reactions        = false;
    bool open_system_constraint = false; // by default use cold-flow constraint

    bool test_tracer_conservation = false;

    enum struct AdvectionType {Invalid, MOL, Godunov};
    AdvectionType m_advection_type = AdvectionType::Invalid;

    //std::string m_advection_type = "MOL";
    std::string m_redistribution_type = "FluxRedist";

    // include d(ep_g)/dt in constraint
    bool m_use_depdt_constraint = false;
    bool m_use_drag_in_projection = false;

    // If using Godunov, default to PPM (not PLM)
    bool m_godunov_ppm         = true;

    // If using Godunov, default to including viscous/diffusive
    //    terms in the forcing terms
    bool m_godunov_include_diff_in_forcing = true;

    // If using Godunov, when to add the forcing terms in
    //    the construction of the "trans" velocities
    bool m_godunov_use_forces_in_trans = false;

    // If using Godunov, if true then don't include the pressure gradient
    //    in the forcing term passed into the godunov routine; instead
    //    use gradient of mac phi which contains the full pressure
    bool m_use_mac_phi_in_godunov = false;

    // If using Godunov, if true then include the explicit fluid
    // drag force, -beta(u_g - u_s), with the pressure term.
    bool m_use_drag_in_godunov = false;

    enum struct DiffusionType {
            Invalid, Explicit, Crank_Nicolson, Implicit
    };
    // MOL: Explicit predictor / Crank_Nicolson corrector
    // Godunov: Implicit
    DiffusionType m_predictor_diff_type = DiffusionType::Invalid;
    DiffusionType m_corrector_diff_type = DiffusionType::Invalid;

    amrex::Real enthalpy_source = 0;


    // Options to control mfix verbosity level
    int m_verbose = 0;

    // ***************************************************************
    // Nodal solve
    // ***************************************************************

    // Verbosity and MLMG parameters are now ParmParse with "nodal_proj" in the
    // inputs file
    // Examples: nodal_proj.verbose = 1
    //           nodal_proj.bottom_verbose = 1
    //           nodal_proj.maxiter
    //           nodal_proj.bottom_maxiter
    //           nodal_proj.bottom_rtol
    //           nodal_proj.bottom_atol
    //           nodal_proj.bottom_solver
    // More info at "amrex/Src/LinearSolvers/Projections/AMReX_NodalProjector.cpp"
    // Verbosity

    // Control MLMG behavior
    Real nodal_mg_rtol = 1.0e-11;
    Real nodal_mg_atol = 1.0e-14;

    // Max coarsening level
    int nodal_mg_max_coarsening_level = 100;

    // Linear operator
    std::shared_ptr<MLNodeLaplacian> nodal_matrix;

    // Solver
    std::shared_ptr<MLMG> nodal_solver;

    // ***************************************************************
    // MAC solve
    // ***************************************************************
    // Verbosity and MLMG parameters are now ParmParse with "mac_proj" in the
    // inputs file
    // Examples: mac_proj.verbose = 1
    //           mac_proj.bottom_verbose = 1
    //           mac_proj.maxiter
    //           mac_proj.bottom_maxiter
    //           mac_proj.bottom_rtol
    //           mac_proj.bottom_atol
    //           mac_proj.bottom_solver
    // More info at "amrex/Src/LinearSolvers/Projections/AMReX_MacProjector.cpp"

    std::unique_ptr<amrex::MacProjector> macproj;

    // Control MLMG behavior
    Real mac_mg_rtol = 1.0e-11;
    Real mac_mg_atol = 1.0e-14;

    // Max coarsening level
    int mac_mg_max_coarsening_level = 100;

    // ***************************************************************

    // Tolerance to check for steady state --
    //    this must be specified in the inputs file
    Real steady_state_tol;

    // Maximum number of iterations to steady state --
    //    this default may be over-written in the inputs file
    int steady_state_maxiter = 100000000;

    int call_udf;

    bool dual_grid = false;

    // If true then we are restarting from a file with no thermal variables
    bool restart_from_cold_flow = false;

    // If true then print the name of the routine we are in
    bool ooo_debug = false;

    // Parameters to control load balancing
    static std::string particle_init_type;
    static std::string load_balance_type;
    static std::string knapsack_weight_type;
    static int load_balance_fluid;
    static int knapsack_nmax;

    // Options to control time stepping
    Real m_cfl = 0.5;
    Real fixed_dt;
    Real dt_min = 0.;    // Min dt allowed
    Real dt_max = 1.e14; // Max dt allowed

    int particle_max_grid_size_x = -1;
    int particle_max_grid_size_y = -1;
    int particle_max_grid_size_z = -1;

    int removeOutOfRange = 1; // Flag to remove out-of-range particles

    static EBSupport m_eb_support_level;

    Periodicity p0_periodicity;

    Real sum_vol_orig;

    Real covered_val = 1.e40;

    //
    // Options to control the computing of particle eulerian velocities
    //
    amrex::Vector<Real> avg_p_g;
    amrex::Vector<Real> avg_ep_g;
    amrex::Vector<Real> avg_vel_g;
    amrex::Vector<Real> avg_T_g;

    amrex::Vector<Real> avg_vel_p;
    amrex::Vector<Real> avg_T_p;

    amrex::Vector<Real> avg_region_x_e;
    amrex::Vector<Real> avg_region_x_w;
    amrex::Vector<Real> avg_region_y_s;
    amrex::Vector<Real> avg_region_y_n;
    amrex::Vector<Real> avg_region_z_b;
    amrex::Vector<Real> avg_region_z_t;

    // Last step at which we wrote a plotfile
    int last_plt = -1;

    // Flags for saving fluid data in plot files
    int plt_vel_g   = 1;
    int plt_ep_g    = 1;
    int plt_p_g     = 0;
    int plt_ro_g    = 0;
    int plt_MW_g    = 0;
    int plt_trac    = 0;
    int plt_cp_g    = 0;
    int plt_T_g     = 0;
    int plt_h_g     = 0;
    int plt_k_g     = 0;
    int plt_mu_g    = 0;
    int plt_X_gk    = 0;
    int plt_D_gk    = 0;
    int plt_cp_gk   = 0;
    int plt_h_gk    = 0;
    int plt_diveu   = 0;
    int plt_volfrac = 0;
    int plt_gradp_g = 0;
    int plt_vort    = 0;
    int plt_ro_txfr = 0;

    // Total number of variables to write in plot file
    int pltVarCount = 0;

    // Flags for saving particle data. By default, we have all flags on,
    // we turn off what we don't want in the init IO routine. This is somewhat
    // different from what we do with the fluid.
    Vector<int> write_real_comp = Vector<int>(AoSrealData::count+SoArealData::count,1);
    Vector<int> write_int_comp = Vector<int>(AoSintData::count+SoAintData::count,1);

    //
    // Here we set the number of ghost nodes of the field arrays.
    // Since this is  a class attribute, it will be propagated to
    // all the member functions.
    //
    // DO NOT USE  A LOCAL VARIABLE WITH THE SAME NAME
    //

    // Four ghost cells are required when using EB
    // const int nghost = 4;
    // const int m_eb_basic_grow_cells = nghost;
    // const int m_eb_volume_grow_cells = nghost;
    // const int m_eb_full_grow_cells = nghost;

    int nghost_eb_basic () const
    {
        if (advection_type() == AdvectionType::Godunov && m_redistribution_type == "StateRedist")
            return 6;
        else if (advection_type() == AdvectionType::Godunov && m_redistribution_type != "StateRedist")
            return 5;
        else if (advection_type() == AdvectionType::MOL && m_redistribution_type == "StateRedist")
            return 5;
        else
            return 4;
    }

    int nghost_eb_volume () const
    {
        if (advection_type() == AdvectionType::Godunov && m_redistribution_type == "StateRedist")
            return 6;
        else if (advection_type() == AdvectionType::Godunov && m_redistribution_type != "StateRedist")
            return 5;
        else if (advection_type() == AdvectionType::MOL && m_redistribution_type == "StateRedist")
            return 5;
        else
            return 4;
    }

    int nghost_eb_full () const {
        if (advection_type() == AdvectionType::Godunov && m_redistribution_type == "StateRedist")
            return 5;
        else
            return 4;
    }

    IntVect particle_sorting_bin = IntVect::TheZeroVector();

    bool do_initial_proj    = true;
    int  initial_iterations = 3;

    // Use tensor solve (as opposed to separate solves for each velocity component)
    bool use_tensor_solve      = true;

    // Use separate solves for each velocity component with lagged tensor terms
    bool use_tensor_correction = false;

     /***************************************************************************
     *                                                                         *
     * I/O stuff                                                               *
     *                                                                         *
     **************************************************************************/

    // Variables to simplify checkpoint IO
    // Vector< Vector< MultiFab* > > vectorVars;
    Vector< std::string > vecVarsName;

    Vector< Vector< MultiFab* > > chkScalarVars;
    Vector< std::string > chkscaVarsName;

    Vector< Vector< MultiFab* > > chkTVars;
    Vector< std::string > chkTVarsName;

    Vector< Vector< MultiFab* > > chkSpeciesVars;
    Vector< std::string > chkSpeciesVarsName;

    Vector< Vector< MultiFab* > > chkSpeciesTVars;
    Vector< std::string > chkSpeciesTVarsName;

    /***************************************************************************
     *                                                                         *
     * Used for fillpatching ...                                               *
     *                                                                         *
     **************************************************************************/

    Vector<Real> t_old;
    Vector<Real> t_new;

    Vector<BCRec> bcs_ls; // used by the level-set fill-patch

    Vector<BCRec> bcs_u; // For velocity components
    Vector<BCRec> bcs_s; // For scalars
    Vector<BCRec> bcs_X; // For fluid species
    Vector<BCRec> bcs_f; // For generic first-order extrapolation

    amrex::Vector<amrex::BCRec> m_bcrec_velocity;
    amrex::Vector<amrex::BCRec> m_bcrec_density;
    amrex::Vector<amrex::BCRec> m_bcrec_enthalpy;
    amrex::Vector<amrex::BCRec> m_bcrec_tracer;
    amrex::Vector<amrex::BCRec> m_bcrec_species;
    amrex::Vector<amrex::BCRec> m_bcrec_force;

    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_velocity_d;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_density_d;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_enthalpy_d;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_tracer_d;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_species_d;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_force_d;


    /***************************************************************************
     *                                                                         *
     * EB Data (including level-set data)                                      *
     *                                                                         *
     **************************************************************************/

    //! EB levels representing fluid boundary conditions
    Vector<const EB2::Level*> eb_levels;

    //! EB levels representing particle boundary conditions (same as
    //! `mfix::eb_levels` but might include additional walls at MI BCs).
    Vector<const EB2::Level*> particle_eb_levels;

    //! EB factory that lives on the fluid grids
    Vector<std::unique_ptr<EBFArrayBoxFactory> > ebfactory;

    amrex::EBFArrayBoxFactory const&
    EBFactory (int lev) const noexcept {
        return static_cast<amrex::EBFArrayBoxFactory const&>(*ebfactory[lev]);
    }


    //! EB factory that lives on the particle grids
    Vector<std::unique_ptr<EBFArrayBoxFactory> > particle_ebfactory;

    bool mfix_update_ebfactory (int a_lev);

  // Number of ghost cells for field arrays.
  int nghost_state () const {
  //if (!EBFactory(0).isAllRegular())
      {
        if (advection_type() == AdvectionType::Godunov && m_redistribution_type == "StateRedist")
          return 6;
        else if (advection_type() == AdvectionType::Godunov && m_redistribution_type != "StateRedist")
          return 5;
        else if (advection_type() == AdvectionType::MOL && m_redistribution_type == "StateRedist")
          return 5;
        else
          return 4;
      }
  //{
  //  return (advection_type() != AdvectionType::MOL) ? 3 : 2;
  //}
  }

    // For Godunov, we need 1 ghost cell in addition to the Box we are filling
    // For MOL    , we need 0 ghost cells
    int nghost_force () const
    {
      if (advection_type() == AdvectionType::MOL)
        return 0;
      else if (m_redistribution_type == "StateRedist")
        return 4;
      else
        return 3;
    }

    int nghost_mac () const {
      if (!EBFactory(0).isAllRegular()) return (advection_type() == AdvectionType::Godunov) ? 4 : 3;
      return (advection_type() == AdvectionType::Godunov) ? 1 : 0;
    }

    AdvectionType advection_type () const{return m_advection_type;}

   // Invalid, Explicit, Crank_Nicolson, Implicit
    DiffusionType predictor_diff_type () const{
      return (advection_type() == AdvectionType::MOL) ?
        DiffusionType::Explicit:
        DiffusionType::Implicit;
    }

    DiffusionType corrector_diff_type () const{
      return (advection_type() == AdvectionType::MOL) ?
        DiffusionType::Crank_Nicolson:
        DiffusionType::Invalid;
    }



    bool need_divtau () const {
      return ( (m_godunov_include_diff_in_forcing ||
                predictor_diff_type() != DiffusionType::Implicit) );
    }

    amrex::Vector<int> m_iconserv_velocity;
    amrex::Gpu::DeviceVector<int> m_iconserv_velocity_d;
    amrex::Vector<int> m_iconserv_density;
    amrex::Gpu::DeviceVector<int> m_iconserv_density_d;
    amrex::Vector<int> m_iconserv_enthalpy;
    amrex::Gpu::DeviceVector<int> m_iconserv_enthalpy_d;
    amrex::Vector<int> m_iconserv_tracer;
    amrex::Gpu::DeviceVector<int> m_iconserv_tracer_d;
    amrex::Vector<int> m_iconserv_species;
    amrex::Gpu::DeviceVector<int> m_iconserv_species_d;


};

#endif
