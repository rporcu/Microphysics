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
#include <hydro_NodalProjector.H>
#include <AMReX_PhysBCFunct.H>
#include <hydro_MacProjector.H>

#include <mfix_diffusion_op.H>
#include <mfix_deposition_op.H>
#include <mfix_pc.H>
#include <mfix_ic.H>
#include <mfix_bc.H>
#include <mfix_eb_if.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_fluid.H>
#include <mfix_reactions.H>
#include <mfix_leveldata.H>
#include <mfix_mlmg_options.H>
#include <mfix_run_on.H>
#include <mfix_rw.H>
#include <mfix_timer.H>

// Mfix run type
struct RunType { enum {Standard, Restart, PIC2DEM}; };

// Drag type enumeration
struct DragType { enum {Invalid=-1, WenYu, Gidaspow, BVK2, SyamOBrien, UserDrag}; };

// Convection type enumeration
struct ConvectionType { enum {Invalid=-1, RanzMarshall, Gunn, NullConvection}; };

// Advection Type
struct AdvectionType { enum {Invalid=-1, MOL, Godunov}; };

// Deposition scheme enumeration
enum DepositionScheme {trilinear, square_dpvm, true_dpvm, centroid};

// Reaction rates type enumeration
struct ReactionRatesType { enum {RRatesUser}; };

enum struct ScalarToFill { Density, Tracer, Enthalpy, VolFrac, Temperature};

class mfix : public AmrCore
{
  public:
    static constexpr Real covered_val = 1.e40;

    // Constructor
    mfix ();

    // Destructor
    ~mfix ();

    // Initial Conditions
    MFIXInitialConditions m_initial_conditions;

    // Boundary conditions
    MFIXBoundaryConditions m_boundary_conditions;

    // Embedded boundaries
    MFIXEmbeddedBoundaries m_embedded_boundaries;

    // Timer
    const MFIXTimer& timer () const { return m_timer; }
    MFIXTimer& timer () { return m_timer; }

    // Regions
    MFIXRegions regions;

    // Chemical species
    MFIXSpecies species;

    // Chemical reactions
    MFIXReactions reactions;

    // Fluid phase
    MFIXFluidPhase fluid;

    // Solids phase
    MFIXSolidsPhase solids;

    // DEM parameters
    MFIXDEM m_dem;

    // PIC parameters
    MFIXPIC m_pic;

    // IO Class
    MFIXReadWrite* m_rw;

    // Face-based coefficients b in MAC projection and implicit diffusion solve
    Vector< Array<MultiFab*,3> > bcoeff;

    void InitParams ();

    bool IsSteadyState ()
    { return (m_steady_state > 0); }

    void Init (Real time,
               const bool init_fluid_grids = true);

    void InitLevelData (Real time);

    void init_advection ();

    void mfix_init_solvers ();

    void mfix_setup_solvers ();

    void PostInit (amrex::Real& dt,
                   const amrex::Real time,
                   int restart_flag,
                   amrex::Real stop_time);

    void ResizeArrays ();

    void InitialRedistribution (amrex::Real l_time);

    void PostProjectionRedistribution (amrex::Real l_time, amrex::Real l_dt,
                                       const amrex::Vector<amrex::MultiFab*>& a_sigma);

    void PreProjectionRedistribution (amrex::Real l_time);

    void PostProjectionDiagnostics (amrex::Real l_time,
                                    amrex::Vector< amrex::MultiFab*> const& epu,
                                    amrex::Vector< amrex::MultiFab*> const& vel_g_in,
                                    amrex::Vector< amrex::MultiFab*> const& p_g_in,
                                    amrex::Vector< amrex::MultiFab*> const& gp_in,
                                    amrex::Vector< amrex::MultiFab*> const& ep_g_in,
                                    amrex::Vector< amrex::MultiFab*> const& S_cc,
                                    bool proj_for_small_dt);

    void Restart (std::string& restart_chkfile,
                  int& nstep,
                  amrex::Real& dt,
                  amrex::Real& time,
                  amrex::IntVect& Nrep);

    void Regrid ();

    void Evolve (int nstep,
                 amrex::Real & dt,
                 amrex::Real & prev_dt,
                 const amrex::Real time,
                 amrex::Real stop_time);

    void restart_pc_from_new_chkpt (const std::string& restart_file,
                                    const std::string& name);

    void mfix_usr0 () const;
    static void mfix_usr1 (Real time);
    void mfix_usr2 () const;
    void mfix_usr3 () const;

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
          if (levelset_refinement == 1) {
              eb_levels[1] = eb_levels[0];
              particle_eb_levels[1] = eb_levels[1];
          } else {
              Geometry geom_ls = amrex::refine(geom[0],levelset_refinement);
              EB2::Build(gshop, geom_ls, 0, 100);
              eb_levels[1] = &(EB2::IndexSpace::top().getLevel(geom_ls));
              particle_eb_levels[1] = eb_levels[1];
          }
      }
    }

    void build_eb_levels_from_chkpt_file ();

    const Geometry get_geom_ref (int lev)
    { return geom[lev]; }

    std::shared_ptr<UnionListIF<EB2::PlaneIF>> get_walls (bool & has_walls);

    std::shared_ptr<UnionListIF<EB2::PlaneIF>> get_real_walls (bool & has_real_walls);

    std::shared_ptr<EB2::TranslationIF<EB2::PolynomialIF>>
      get_poly (int max_order, std::string field_prefix);

    std::shared_ptr<EB2::IntersectionIF<EB2::PlaneIF,EB2::PlaneIF,EB2::PlaneIF>>
      make_wall (int dir, // direction (long edge) of wall
                 Real position, Real height, Real width);

    static inline std::string get_load_balance_type ()
    { return load_balance_type; }

    // flag enabling level-set restart (i.e. prevent make_eb_* from rebuilding
    // the level-set data).
    bool levelset_restart = false;

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

    void mfix_calc_txfr_particle (Real time,
                                  Vector< MultiFab* > const& ep_g_in,
                                  Vector< MultiFab* > const& ro_g_in,
                                  Vector< MultiFab* > const& vel_g_in,
                                  Vector< MultiFab* > const& T_g_in,
                                  Vector< MultiFab* > const& X_gk_in,
                                  Vector< Real* > const& pressure_g_in,
                                  Vector< MultiFab* > const& gp_in);

    template <typename F1>
    void mfix_calc_txfr_particle (Real time,
                                  Vector< MultiFab* > const& ep_g_in,
                                  Vector< MultiFab* > const& ro_g_in,
                                  Vector< MultiFab* > const& vel_g_in,
                                  Vector< MultiFab* > const& T_g_in,
                                  Vector< MultiFab* > const& X_gk_in,
                                  Vector< Real* > const& pressure_g_in,
                                  Vector< MultiFab* > const& gp_in,
                                  F1 HeterogeneousRatesFunc);

    void mfix_add_vel_txfr_explicit (Real dt,
                                     amrex::Vector<amrex::MultiFab*      > const& vel_in,
                                     amrex::Vector<amrex::MultiFab const*> const& txfr_in,
                                     amrex::Vector<amrex::MultiFab const*> const& rho_in,
                                     amrex::Vector<amrex::MultiFab const*> const& ep_g_in);

    void mfix_add_enthalpy_txfr_explicit (Real dt,
                                          amrex::Vector<amrex::MultiFab*      > const& h_g_in,
                                          amrex::Vector<amrex::MultiFab*      > const& T_g_in,
                                          amrex::Vector<amrex::MultiFab const*> const& X_gk_in,
                                          amrex::Vector<amrex::MultiFab const*> const& txfr_in,
                                          amrex::Vector<amrex::MultiFab const*> const& rho_in,
                                          amrex::Vector<amrex::MultiFab const*> const& ep_g_in);

    void mfix_add_vel_txfr_implicit (amrex::Real dt,
                                     amrex::Vector<amrex::MultiFab*      > const& vel_in,
                                     amrex::Vector<amrex::MultiFab const*> const& txfr_in,
                                     amrex::Vector<amrex::MultiFab const*> const& rho_in,
                                     amrex::Vector<amrex::MultiFab const*> const& ep_g_in);

    void mfix_add_enthalpy_txfr_implicit (amrex::Real dt,
                                          amrex::Vector<amrex::MultiFab*      > const& h_g_in,
                                          amrex::Vector<amrex::MultiFab*      > const& T_g_in,
                                          amrex::Vector<amrex::MultiFab const*> const& X_gk_in,
                                          amrex::Vector<amrex::MultiFab const*> const& txfr_in,
                                          amrex::Vector<amrex::MultiFab const*> const& rho_in,
                                          amrex::Vector<amrex::MultiFab const*> const& ep_g_in);

    void mfix_correct_small_cells (Vector<MultiFab*      > const& vel_in,
                                   Vector<MultiFab const*> const& ep_u_mac,
                                   Vector<MultiFab const*> const& ep_v_mac,
                                   Vector<MultiFab const*> const& ep_w_mac,
                                   Vector<MultiFab const*> const& eb_vel);

    void set_ls_inflow (const int lev,
                        FArrayBox& ls_phi_fab,
                        const Box& domain,
                        const int* levelset_nghost,
                        const int& nref,
                        const Real* dx);

    void set_p0 (const Box& bx,
                 MFIter* mfi,
                 const int lev,
                 const Box& domain);

    void set_gp0 (const int lev,
                  const Box& domain);

    static int m_run_type;

    static int m_drag_type;

    static amrex::Real m_SyamOBrien_coeff_c1;
    static amrex::Real m_SyamOBrien_coeff_d1;

    static int m_convection_type;

    static DepositionScheme m_deposition_scheme;

    static int m_reaction_rates_type;

    static amrex::Real m_deposition_diffusion_coeff;
    static amrex::Real m_deposition_scale_factor;
    static amrex::Real m_max_solids_volume_fraction;

    template <typename F1>
    void mfix_calc_transfer_coeffs (const amrex::Real time,
                                    amrex::Vector< amrex::MultiFab* > const& ep_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& ro_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& vel_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& T_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& X_gk_in,
                                    amrex::Vector< amrex::Real* > const& pressure_g_in,
                                    F1 DragFunc);

    template <typename F1, typename F2>
    void mfix_calc_transfer_coeffs (const amrex::Real time,
                                    amrex::Vector< amrex::MultiFab* > const& ep_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& ro_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& vel_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& T_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& X_gk_in,
                                    amrex::Vector< amrex::Real* > const& pressure_g_in,
                                    F1 DragFunc,
                                    F2 ConvectionCoeff);

    template <typename F1, typename F2, typename F3>
    void mfix_calc_transfer_coeffs (const amrex::Real time,
                                    amrex::Vector< amrex::MultiFab* > const& ep_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& ro_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& vel_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& T_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& X_gk_in,
                                    amrex::Vector< amrex::Real* > const& pressure_g_in,
                                    F1 DragFunc,
                                    F2 ConvectionCoeff,
                                    F3 HeterogeneousRRates);

    void mfix_calc_volume_fraction (const amrex::Real time,
                                    amrex::Real & sum_vol);

    void mfix_deposition_bcs (int lev,
                              amrex::MultiFab & filled_mf) const;

    void mfix_redistribute_deposition (int lev,
                                       amrex::MultiFab & mf_eps,
                                       amrex::MultiFab & mf_to_redistribute,
                                       const amrex::MultiFab * volfrac,
                                       const amrex::FabArray<EBCellFlagFab>* flags,
                                       amrex::Real max_eps) const;

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

    void compute_MAC_projected_velocities (amrex::Real time,
                                           const amrex::Real l_dt,
                                           amrex::Vector< amrex::MultiFab const*> const& vel_in,
                                           amrex::Vector< amrex::MultiFab*      > const& ep_u_mac,
                                           amrex::Vector< amrex::MultiFab*      > const& ep_v_mac,
                                           amrex::Vector< amrex::MultiFab*      > const& ep_w_mac,
                                           amrex::Vector< amrex::MultiFab const*> const& ep_g_in,
                                           amrex::Vector< amrex::MultiFab const*> const& ro_g_in,
                                           amrex::Vector< amrex::MultiFab const*> const& txfr_in,
                                           amrex::Vector< amrex::MultiFab const*> const& eb_vel,
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

    std::shared_ptr<Hydro::NodalProjector> nodal_projector;

    ////////////////////////////////////////////////////////////////////////////
    //
    // Deposition
    //
    ////////////////////////////////////////////////////////////////////////////

    MFIXSolidsVolume* m_solids_volume_deposition;
    MFIXInterphaseTxfr* m_interphase_txfr_deposition;

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
                        int compute_mass_balance);

    void pic_iteration (const bool apply_forces,
                        const bool update_parcels,
                        const bool use_taylor_approx,
                        const Real advance_vel_p,
                        amrex::Real dt,
                        amrex::RealVect& gravity,
                        amrex::Vector< amrex::Array<amrex::MultiFab*,3> >& vel_s_in,
                        amrex::Vector< amrex::MultiFab* >& ep_s_out,
                        amrex::EBFArrayBoxFactory* ebfactory,
                        const int ls_refinement,
                        const MultiFab* ls_phi);

    void MFIX_CalcAvgSolidsVel (Vector< Array<MultiFab*,3> >& vel_s,
                                const bool do_deposition=false);

    void MFIX_CalcSolidsStress (amrex::Vector< amrex::MultiFab* >& ep_s_in,
                                amrex::RealVect& gravity_in,
                                amrex::Vector< amrex::MultiFab* >& cost,
                                std::string& knapsack_weight_type);


    void mfix_set_pic_velocity_bcs (Vector< MultiFab* > const& vel) const;

    void set_pic_velocity_bcs (const int lev,
                               FArrayBox& vel_fab,
                               const Box& domain) const;

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

    void mfix_init_fluid (int is_restarting,
                          const amrex::Real time,
                          amrex::Real dt,
                          amrex::Real stop_time);

    void mfix_set_bc0 ();

    void mfix_set_p0 ();

    void mfix_set_gp0 ();

    static RealVect gp0;
    static RealVect gravity;

    // Scalar diffusive coefficient
    Vector<Real> mu_s;

    BCList bc_list;

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
                                      amrex::Vector<amrex::MultiFab const*> const& density,
                                      amrex::Vector< amrex::MultiFab const*> const& eb_vel);

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
                       const bool update_flux,
                       amrex::Vector< amrex::MultiFab*      > const& lapT,
                       amrex::Vector< amrex::MultiFab*      > const& lapS,
                       amrex::Vector< amrex::Array< amrex::MultiFab*, AMREX_SPACEDIM> > const& J_gk,
                       amrex::Vector< amrex::MultiFab*      > const& T_g_in,
                       amrex::Vector< amrex::MultiFab*      > const& trac_in,
                       amrex::Vector< amrex::MultiFab*      > const& X_gk_in,
                       amrex::Vector< amrex::MultiFab const*> const& ep_g_in,
                       amrex::Vector< amrex::MultiFab const*> const& ro_g_in);

  protected:
    void FillPatchVel (int lev, Real time, MultiFab & mf, int icomp, int ncomp,
                       const Vector<BCRec> & bcr);

    void GetDataVel (int lev, Real time, Vector<MultiFab *> & data,
                     Vector<Real> & datatime);

    void FillPatchScalar (int lev, Real time, MultiFab & mf, ScalarToFill scalar,
                        const Real* bc_scalar, const Vector<BCRec> & bcr);

    void FillPatchSpecies (int lev, Real time, MultiFab & mf, int icomp, int ncomp,
                           const Vector<BCRec> & bcr);

    void GetDataScalar (int lev, Real time, Vector<MultiFab *> & data,
                        ScalarToFill scalar_id, Vector<Real> & datatime);

    void GetDataSpecies (int lev, Real time, Vector<MultiFab *> & data, int icomp,
                         Vector<Real> & datatime);


    void fillpatch_force (amrex::Real time, amrex::Vector<amrex::MultiFab*> const& force, int ng);

    void mfix_calc_transfer_coeffs (const amrex::Real time,
                                    amrex::Vector< amrex::MultiFab* > const& ep_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& ro_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& vel_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& T_g_in,
                                    amrex::Vector< amrex::MultiFab* > const& X_gk_in,
                                    amrex::Vector< amrex::Real* > const& pressure_g_in);

    void mfix_calc_txfr_fluid (Vector< MultiFab* > const& txfr_out,
                               Vector< MultiFab* > const& ep_g_in,
                               Vector< MultiFab* > const& ro_g_in,
                               Vector< MultiFab* > const& vel_g_in,
                               Vector< MultiFab* > const& T_g_in,
                               Vector< MultiFab* > const& X_gk_in,
                               Vector< Real* > const& pressure_g_in,
                               Real time);

  public:
    void mfix_deposit_particles (MFIXDepositionOp* deposition_op,
                                 Vector< MultiFab* > const& txfr_out,
                                 const Real time) const;

    void AllocateArrays (int lev);

  protected:
    void RegridArrays (int lev);

    void RegridLevelSetArray (int a_lev);

    void avgDown (int crse_lev, const MultiFab & S_fine, MultiFab & S_crse);

  private:
    void make_eb_box ();
    void make_eb_cylinder ();
    void make_eb_hopper ();
    void make_eb_regular ();
    void make_eb_generic ();

#ifdef CSG_EB
  void make_eb_csg(std::string);
#endif

    void MakeNewLevelFromScratch (int lev, Real time, const BoxArray & new_grids,
                                  const DistributionMapping & new_dmap) override;

    void ReMakeNewLevelFromScratch (int lev, const BoxArray & new_grids,
                                    const DistributionMapping & new_dmap);

    void EvolveFluid (int nstep,
                      amrex::Real & dt,
                      amrex::Real & prev_dt,
                      const amrex::Real time,
                      amrex::Real stop_time,
                      amrex::Real & drag_timing);

    void check_for_nans (int lev);

    void fillpatch_all (Vector< MultiFab* > const& vel_in,
                        Vector< MultiFab* > const& ro_g_in,
                        Vector< MultiFab* > const& h_g_in,
                        Vector< MultiFab* > const& trac_in,
                        Vector< MultiFab* > const& X_gk_in,
                        Real time);

    int steady_state_reached (Real dt, int iter);

  public:
    void mfix_compute_convective_term (Vector< MultiFab*      >& conv_u,
                                       Vector< MultiFab*      >& conv_s,
                                       Vector< MultiFab*      >& conv_X,
                                       Vector< MultiFab*      > const& vel_forces,
                                       Vector< MultiFab*      > const& tra_forces,
                                       Vector< MultiFab const*> const& vel_in,
                                       Vector< MultiFab*      > const& ep_g_in,
                                       Vector< MultiFab const*> const& ro_g_in,
                                       Vector< MultiFab const*> const& h_g_in,
                                       Vector< MultiFab const*> const& trac_in,
                                       Vector< MultiFab const*> const& X_gk_in,
                                       Vector< MultiFab const*> const& txfr_in,
                                       Vector< MultiFab const*> const& eb_vel,
                                       Vector< MultiFab const*> const& eb_scalars,
                                       Vector< MultiFab const*> const& eb_species,
                                       Vector< MultiFab*      > const& u_mac,
                                       Vector< MultiFab*      > const& v_mac,
                                       Vector< MultiFab*      > const& w_mac,
                                       Real l_dt, Real time);

    void mfix_apply_predictor (Vector< MultiFab* >& conv_u_old,
                               Vector< MultiFab* >& conv_s_old,
                               Vector< MultiFab* >& conv_X_old,
                               Vector< MultiFab* >& ro_RHS_old,
                               Vector< MultiFab* >& lap_trac_old,
                               Vector< MultiFab* >& enthalpy_RHS_old,
                               Vector< MultiFab* >& lap_T_old,
                               Vector< MultiFab* >& species_RHS_old,
                               Vector< MultiFab* >& vel_RHS_old,
                               Vector< MultiFab* >& div_J_old,
                               Vector< MultiFab* >& div_hJ_old,
                               Vector< Real >& rhs_pressure_g_old,
                               Vector< Real >& rhs_pressure_g,
                               Vector< MultiFab* > eb_flow_vel,
                               Vector< MultiFab* > eb_flow_scalars,
                               Vector< MultiFab* > eb_flow_species,
                               Real time,
                               Real dt,
                               Real prev_dt,
                               bool proj_2,
                               Real& coupling_timing);

    void mfix_density_rhs (Vector< MultiFab*       > const& rhs,
                           Vector< MultiFab const* > const& txfr);

    void mfix_enthalpy_rhs(Vector< MultiFab*      > const& rhs,
                           Vector< MultiFab const*> const& ep_g,
                           Vector< MultiFab const*> const& ro_g,
                           Vector< MultiFab*      > const& X_gk,
                           Vector< MultiFab const*> const& T_g,
                           Vector< MultiFab const* > const& txfr);

    void mfix_species_X_rhs (Vector< MultiFab*       > const& rhs,
                             Vector< MultiFab const* > const& txfr);

    void mfix_momentum_rhs (Vector< MultiFab*       > const& rhs,
                            Vector< MultiFab const* > const& ep_g,
                            Vector< MultiFab const* > const& txfr);

    void mfix_apply_corrector (Vector< MultiFab* >& conv_u_old,
                               Vector< MultiFab* >& conv_s_old,
                               Vector< MultiFab* >& conv_X_old,
                               Vector< MultiFab* >& ro_RHS_old,
                               Vector< MultiFab* >& lap_trac_old,
                               Vector< MultiFab* >& enthalpy_RHS_old,
                               Vector< MultiFab* >& lap_T_old,
                               Vector< MultiFab* >& species_RHS_old,
                               Vector< MultiFab* >& vel_RHS_old,
                               Vector< MultiFab* >& div_J_old,
                               Vector< MultiFab* >& div_hJ_old,
                               Vector< Real >& rhs_pressure_g_old,
                               Vector< Real >& rhs_pressure_g,
                               Vector< MultiFab* > eb_flow_vel,
                               Vector< MultiFab* > eb_flow_scalars,
                               Vector< MultiFab* > eb_flow_species,
                               Real time,
                               Real dt,
                               Real prev_dt,
                               bool proj_2,
                               Real& coupling_timing);

  public:
    void mfix_incompressible_fluid_rhs (Vector< MultiFab* > const& rhs);

    void mfix_idealgas_opensystem_rhs (Vector< MultiFab*       > const& rhs,
                                       Vector< MultiFab const* > const& enthalpy_rhs,
                                       Vector< MultiFab const* > const& species_rhs,
                                       Vector< MultiFab const* > const& ro_g,
                                       Vector< MultiFab const* > const& T_g,
                                       Vector< MultiFab const* > const& X_gk);

    void mfix_idealgas_closedsystem_rhs (Vector< MultiFab*       > const& rhs,
                                         Vector< MultiFab const* > const& enthalpy_rhs,
                                         Vector< MultiFab const* > const& species_rhs,
                                         Vector< MultiFab const* > const& ep_g,
                                         Vector< MultiFab const* > const& ro_g,
                                         Vector< MultiFab const* > const& T_g,
                                         Vector< MultiFab const* > const& X_gk,
                                         Vector< Real const* > const& pressure_g,
                                         Vector< Real >& avgSigma,
                                         Vector< Real >& avgTheta);

    void compute_convective_term (Box const& bx, int lev, const Real l_dt, MFIter const& mfi,
                                  Array4<Real> const& dvdt, // velocity
                                  Array4<Real> const& dsdt, // density, enthalpy, tracer
                                  Array4<Real> const& dXdt, // tracer
                                  Array4<Real const> const& vel,
                                  Array4<Real const> const& ep_g,
                                  Array4<Real const> const& rho,
                                  Array4<Real const> const& hg,
                                  Array4<Real const> const& tra,
                                  Array4<Real const> const& Xgk,
                                  Array4<Real const> const& divu,
                                  Array4<Real const> const& ep_umac,
                                  Array4<Real const> const& ep_vmac,
                                  Array4<Real const> const& ep_wmac,
                                  Array4<Real const> const& fvel,
                                  Array4<Real const> const& ftra);

    amrex::IntVect mfix_locate_max_eps ( Vector< MultiFab*>& ep_s_in, const Real max_eps );

  private:
    // Make a level 0 grids covering the whole domain.  It does NOT install the new grids.
    BoxArray MakeBaseGrids (const Box& domain,
                            const IntVect& grid_sizes) const;

    void PruneBaseGrids(BoxArray &ba) const;

    void ChopGrids (const Box& domain, BoxArray& ba, int target_size) const;

  public:
    static int nlev;
    static int ntrac;

    // Particle container
    MFIXParticleContainer* pc;

    Vector< std::unique_ptr<LevelData> > m_leveldata;

    void PIC_to_PIC(const int lev,
                    std::string& restart_file,
                    const int pic_factor);

    void PIC_to_DEM(const int lev);

  private:
    // Debug counter used for tracking number of level-set MultiFabs used in
    // intersection/union
    int ct_ls_mf;

    MFIXParticleContainer* getParticleContainer() const { return pc; }

    Vector< MultiFab* > get_ep_g () noexcept;
    Vector< MultiFab* > get_ro_g () noexcept;
    Vector< MultiFab* > get_ro_g_old () noexcept;
    Vector< MultiFab* > get_trac () noexcept;
    Vector< MultiFab* > get_trac_old () noexcept;
    Vector< MultiFab* > get_vel_g () noexcept;
    Vector< MultiFab* > get_vel_g_old () noexcept;
    Vector< MultiFab* > get_gp () noexcept;
    Vector< MultiFab* > get_p_g () noexcept;
    Vector< MultiFab* > get_p_g_old () noexcept;
    Vector< Real* > get_thermodynamic_p_g () noexcept;
    Vector< Real* > get_thermodynamic_p_g_old () noexcept;
    Vector< MultiFab* > get_T_g () noexcept;
    Vector< MultiFab* > get_T_g_old () noexcept;
    Vector< MultiFab* > get_h_g () noexcept;
    Vector< MultiFab* > get_h_g_old () noexcept;
    Vector< MultiFab* > get_T_g_on_eb () noexcept;
    Vector< MultiFab* > get_X_gk () noexcept;
    Vector< MultiFab* > get_X_gk_old () noexcept;
    Vector< MultiFab* > get_txfr () noexcept;
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

    amrex::Vector< Real const* > get_thermodynamic_p_g_const () const noexcept;
    amrex::Vector< Real const* > get_thermodynamic_p_g_old_const () const noexcept;

    amrex::Vector< MultiFab const* > get_divtau_const () const noexcept;

    amrex::Vector< amrex::MultiFab const* > get_T_g_on_eb_const () const noexcept;
    amrex::Vector< amrex::MultiFab const* > get_txfr_const () const noexcept;
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
    amrex::Vector<amrex::BCRec> const& get_hydro_velocity_bcrec () const noexcept { return m_bcrec_hydro_velocity; }
    amrex::Vector<amrex::BCRec> const& get_density_bcrec () const noexcept { return m_bcrec_density; }
    amrex::Vector<amrex::BCRec> const& get_enthalpy_bcrec () const noexcept { return m_bcrec_enthalpy; }
    amrex::Vector<amrex::BCRec> const& get_tracer_bcrec () const noexcept { return m_bcrec_tracer; }
    amrex::Vector<amrex::BCRec> const& get_species_bcrec () const noexcept { return m_bcrec_species; }
    amrex::Vector<amrex::BCRec> const& get_force_bcrec () const noexcept { return m_bcrec_force; }
    //
    amrex::BCRec const* get_velocity_bcrec_device_ptr () const noexcept { return m_bcrec_velocity_d.data(); }
    amrex::BCRec const* get_hydro_velocity_bcrec_device_ptr () const noexcept { return m_bcrec_hydro_velocity_d.data(); }
    amrex::BCRec const* get_density_bcrec_device_ptr () const noexcept { return m_bcrec_density_d.data(); }
    amrex::BCRec const* get_enthalpy_bcrec_device_ptr () const noexcept { return m_bcrec_enthalpy_d.data(); }
    amrex::BCRec const* get_tracer_bcrec_device_ptr () const noexcept { return m_bcrec_tracer_d.data(); }
    amrex::BCRec const* get_species_bcrec_device_ptr () const noexcept { return m_bcrec_species_d.data(); }
    amrex::BCRec const* get_force_bcrec_device_ptr () const noexcept { return m_bcrec_force_d.data(); }

    // used if load_balance_type == "KnapSack"
    amrex::Vector< amrex::MultiFab* > particle_cost;
    amrex::Vector< amrex::MultiFab* > particle_proc;
    amrex::Vector< amrex::MultiFab* > fluid_cost;
    amrex::Vector< amrex::MultiFab* > fluid_proc;

    // set distribution map manually, used for profiling
    amrex::Vector<int> pmap;
    amrex::Vector<int> particle_pmap;

    // Level-Set Data => used for particle-wall collisions and fluid
    // reconstruction in particle drag calculation. NOTE: this has at least 2
    // levels: 0) fluid reconstruction, 1..N) particle-wall (refined) collisions
    Vector< std::unique_ptr<MultiFab> > level_sets;

    int m_steady_state = 0;

    bool test_tracer_conservation = false;

    static int m_advection_type;

    std::string m_redistribution_type = "StateRedist";

    bool m_grid_pruning = false;

    // include d(ep_g)/dt in constraint
    bool m_use_depdt_constraint = false;
    bool m_use_drag_in_projection = false;

    // Redistribute the velocity field before the nodal projection
    bool m_redistribute_before_nodal_proj = true;

    // Redistribute the velocity field after the nodal projection
    bool m_redistribute_nodal_proj = false;

    // Tolerances and max iterations for Newton solver
    amrex::Real newton_reltol = 1.e-8;
    amrex::Real newton_abstol = 1.e-8;
    int newton_maxiter = 500;

    // If using Godunov, default to PLM (not PPM)
    bool m_godunov_ppm = false;

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
    int update_mass             = 1;
    int update_momentum         = 1;
    int update_enthalpy         = 1;


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
    // More info at "AMReX-Hydro/Projections/hydro_NodalProjector.cpp"
    // Control MLMG behavior
    std::unique_ptr<MfixUtil::MLMGOptions> nodalproj_options;

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
    // More info at "AMReX-Hydro/Projections/hydro_MacProjector.cpp"
    std::unique_ptr<Hydro::MacProjector> macproj;
    // Control MLMG behavior
    std::unique_ptr<MfixUtil::MLMGOptions> macproj_options;

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
    static std::string load_balance_type;
    static std::string knapsack_weight_type;
    static int load_balance_fluid;
    static int knapsack_nmax;
    static int  greedy_dir;
    static bool greedy_3d;
    static int  greedy_min_grid_size;

    // Options to control time stepping
    Real m_cfl = 0.5;
    MFIXTimer m_timer;

    int particle_max_grid_size_x = -1;
    int particle_max_grid_size_y = -1;
    int particle_max_grid_size_z = -1;
    IntVect particle_max_grid_size;

    int agg_grid_size = -1; // agglomeration for GMG

    int removeOutOfRange = 1; // Flag to remove out-of-range particles

    int sort_particle_int = -1; // # steps to sort particles

#ifdef AMREX_USE_GPU
    bool reduceGhostParticles = true;
#else
    bool reduceGhostParticles = false;
#endif

    Real overload_toler  = 1.2;
    Real underload_toler = 0.8;

    static EBSupport m_eb_support_level;

    Periodicity p0_periodicity;

    Real sum_vol_orig;

    Real init_value = 9.87654321e32;


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

    int nghost_eb_basic () const  { return 4; }

    int nghost_eb_volume () const { return 5; }

    int nghost_eb_full () const   { return 4; }

    IntVect particle_sorting_bin = IntVect::TheUnitVector();

    bool do_initial_proj    = true;
    int  initial_iterations = 3;

    // Use tensor solve (as opposed to separate solves for each velocity component)
    bool use_tensor_solve      = true;

    // Use separate solves for each velocity component with lagged tensor terms
    bool use_tensor_correction = false;


    /***************************************************************************
     *                                                                         *
     * Used for fillpatching ...                                               *
     *                                                                         *
     **************************************************************************/

    Vector<Real> t_old;
    Vector<Real> t_new;

    Vector<BCRec> bcs_ls; // used by the level-set fill-patch
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

    // This is a special set of bc's that replace foextrap by extdir
    // for *normal* velocity only at faces that are identified as
    // pressure inflow (pinf_).   We need this because otherwise the
    // hydro routines will enforce "no backflow at outflow" and will
    // zero the inflow velocity at a pressure inflow
    amrex::Vector<amrex::BCRec> m_bcrec_hydro_velocity;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_hydro_velocity_d;

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

    //! EB factory that lives on the particle grids
    Vector<std::unique_ptr<EBFArrayBoxFactory> > particle_ebfactory;

  public:
    amrex::EBFArrayBoxFactory const&
    EBFactory (int lev) const noexcept {
        return static_cast<amrex::EBFArrayBoxFactory const&>(*ebfactory[lev]);
    }

    amrex::EBFArrayBoxFactory const&
    ParticleEBFactory (int lev) const noexcept {
        return static_cast<amrex::EBFArrayBoxFactory const&>(*particle_ebfactory[lev]);
    }

  private:
    bool mfix_update_ebfactory (int a_lev);

    // Number of ghost cells for field arrays.
    int nghost_state () const { return 4; }

    // For Godunov, we need 1 ghost cell in addition to the Box we are filling
    // For MOL    , we need 0 ghost cells
    int nghost_force () const {
       if (advection_type() == AdvectionType::MOL)
           return 0;
       else
           return 1;
    }

    int nghost_mac () const {
      if (!EBFactory(0).isAllRegular()) return (advection_type() == AdvectionType::Godunov) ? 4 : 3;
      return (advection_type() == AdvectionType::Godunov) ? 1 : 0;
    }

    int advection_type () const{return m_advection_type;}

   // Invalid, Explicit, Crank_Nicolson, Implicit
    DiffusionType predictor_diff_type () const{
      return m_predictor_diff_type;
    }

    DiffusionType corrector_diff_type () const{
      return m_corrector_diff_type;
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
