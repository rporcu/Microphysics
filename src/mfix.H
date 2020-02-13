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
#include <AMReX_EB_levelset.H>
#include <AMReX_EB_LSCore.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_NodalProjector.H>
#include <AMReX_PhysBCFunct.H>

#include <DiffusionOp.H>
#include <MFIXParticleContainer.H>
#include <MFIX_BcList.H>
#include <mfix_des_F.H>
#include <mfix_eb_if.H>
#include <MFIX_DEM_Parms.H>


enum DragType {Invalid=-1, WenYu, Gidaspow, BVK2, UserDrag};

enum DepositionScheme {trilinear, square_dpvm, true_dpvm, centroid};

class mfix : public AmrCore
{
  public:
    // Constructor
    mfix ();

    // Destructor
    ~mfix ();

    // Face-based coefficients b in MAC projection and implicit diffusion solve
    Vector< Array< std::unique_ptr<MultiFab>, AMREX_SPACEDIM> > bcoeff;

    void InitParams ();

    bool IsSteadyState ()
    { return (steady_state > 0); }

    void Init (Real time);

    void InitLevelData (Real time);

    void mfix_init_solvers ();

    void mfix_setup_solvers ();

    void PostInit (Real& dt, Real time, int restart_flag, Real stop_time);

    void ResizeArrays ();

    void WriteCheckPointFile (std::string & check_file_name,
                              int nstep = 0, Real dt = 0.0, Real time = 0.0)  const;

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

    void mfix_usr1_cpp (Real time) const;

    void usr3 ();

    void output (int estatus, int finish, int nstep, Real dt, Real time)
    { pc->output(estatus, finish, nstep, dt, time); };

    void make_mg_bc (int mg_bc[]);

    void make_eb_geometry ();

    void make_eb_factories ();

    void fill_eb_levelsets ();

    void intersect_ls_walls ();

    template<class F>
    void build_particle_eb_levels (EB2::GeometryShop<F> gshop)
    {
      if(DEM::solve)
        for (int lev = 0; lev < nlev; lev ++)
          particle_eb_levels[lev] = LSCore<F>::BuildEBLevel(gshop, lev, geom, 100, 0);
    }

    //! Construct EB levels from Geometry shop. This builds each EB level and
    //! saves a pointer to each level into `eb_levels`.
    template<class F> void build_eb_levels (EB2::GeometryShop<F> gshop)
    {
      /************************************************************************
       *                                                                      *
       * Build EB levels                                                      *
       *                                                                      *
       ***********************************************************************/
      for (int lev = 0; lev < nlev; lev ++)
      {
        // NOTE: the last two parameters are the max coarsening level and
        // the required coarsening level. They are also used by the Poisson
        // solver, => changing these will result in slight changes of the
        // fluid solve.
        eb_levels[lev] = LSCore<F>::BuildEBLevel(gshop, lev, geom, 100, 0);
        particle_eb_levels[lev] = eb_levels[lev];
      }

      if (nlev == 1)
      {
        Box dom_ls = geom[0].Domain();
        dom_ls.refine(levelset_refinement);
        Geometry geom_ls(dom_ls);
        Vector<Geometry> vgeom{geom_ls};
        // NOTE: lev here is index in vgeom
        eb_levels[1] = LSCore<F>::BuildEBLevel(gshop, 0, vgeom, 100, 0);
        particle_eb_levels[1] = eb_levels[1];
      }
    }

    const Geometry get_geom_ref (int lev)
    { return geom[lev]; }

    std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_walls (bool & has_walls);

    std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_real_walls (bool & has_real_walls);

    std::unique_ptr<EB2::TranslationIF<EB2::PolynomialIF>> 
      get_poly (int max_order, std::string field_prefix);

    std::unique_ptr<EB2::IntersectionIF<EB2::PlaneIF,EB2::PlaneIF,EB2::PlaneIF>>
      make_wall (int dir, // direction (long edge) of wall
                 Real position, Real height, Real width);

    void WriteMyEBSurface ();

    void InitIOChkData ();

    void InitIOPltData ();

    // If > 0, dump plotfiles every "plot_int" steps
    static int plot_int;

    // If > 0, dump plotfiles at as close as possible to the designated period
    // *without* changing dt
    static Real plot_per_approx;

    // If > 0, dump plotfiles at exactcly the designated period by changing dt
    static Real plot_per_exact;

    static std::string get_load_balance_type ();

    // flag enabling level-set restart (i.e. prevent make_eb_* from rebuilding
    // the level-set data).
    bool levelset_restart = false;

    void mfix_compute_vort ();

    void mfix_compute_dt (int nstep, Real time, Real stop_time, Real & dt);

    ///
    /// Parameters describing the level-set grid:
    ///
    int levelset_refinement = 1;    // refinement (wrt particle grid) of the level-set's grid
    int levelset_eb_refinement = 1; // refinement of the EB facets used to generate the level-set
    int levelset_pad = 2;           // padding of the level-set grid
    int levelset_eb_pad = 2;        // padding of the EBIS used to generate the level-set (useful for
                                    // ensuring that neighbor grids are considered when filling LS)

    bool contains_ebs = false;

    void mfix_add_drag_explicit (Real dt);

    void mfix_add_drag_implicit (Real dt);

    void mfix_calc_drag_particle (Real time);

    void mfix_add_gravity_and_gp (Real dt);

    void mfix_predict_vels_on_faces (int lev,
                                     Real time,
                                     Vector< std::unique_ptr<MultiFab> >& vel,
                                     Vector< std::unique_ptr<MultiFab> >& ep_u_mac,
                                     Vector< std::unique_ptr<MultiFab> >& ep_v_mac,
                                     Vector< std::unique_ptr<MultiFab> >& ep_w_mac,
                                     Vector< std::unique_ptr<MultiFab> >& ep_g);

    void mfix_set_velocity_bcs (Real time,
                                Vector< std::unique_ptr<MultiFab> >& vel,
                                int extrap_dir_bcs) const;

    void mfix_set_epg_bcs (const amrex::Vector< std::unique_ptr<MultiFab> >& mf_eps) const;

    void mfix_set_density_bcs (Real time,
                               Vector< std::unique_ptr<MultiFab> >& ro_g);

    void mfix_set_scalar_bcs (Real time,
                              Vector< std::unique_ptr<MultiFab> >& trac,
                              Vector< std::unique_ptr<MultiFab> >& mu_g);

    void set_velocity_bcs (Real time,
                           const int lev,
                           FArrayBox& vel_fab,
                           const Box& domain,
                           const int* extrap_dir_bcs) const;

    void set_epg_bcs (const int lev,
                      FArrayBox& eps_fab,
                      const Box& domain,
                      const int* extrap_dir_bcs) const;

    void set_density_bcs (Real time,
                          const int lev,
                          FArrayBox& scal_fab,
                          const Box& domain);

    void set_scalar_bcs (Real time,
                         const int lev,
                         FArrayBox& scal_fab,
                         const int comp,
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

    void mfix_compute_slopes (int lev, Real time, MultiFab& Sborder,
                              Vector<std::unique_ptr<MultiFab>>& xslopes_in,
                              Vector<std::unique_ptr<MultiFab>>& yslopes_in,
                              Vector<std::unique_ptr<MultiFab>>& zslopes_in,
                              int slopes_comp);

    void set_bc0 (const Box& sbx,
                  MFIter* mfi,
                  const int lev,
                  const Box& domain);

    void set_p0 (const Box& bx,
                 MFIter* mfi,
                 const int lev,
                 const Box& domain);

    void set_gp0 (const int lev,
                  const Box& domain);

    static DragType m_drag_type;

    static DepositionScheme m_deposition_scheme;

    static amrex::Real m_deposition_diffusion_coeff;
    static amrex::Real m_deposition_scale_factor;
    static amrex::Real m_max_solids_volume_fraction;

    // input parameter defaulted to 50 used in file src_des/init_collision.f90
    static amrex::Real tcoll_ratio; 

    template <typename F>
    void mfix_calc_particle_beta (F DragFunc, Real time);

    void mfix_calc_volume_fraction (Real & sum_vol);

    void mfix_deposition_bcs_scalar (int lev, amrex::MultiFab & filled_mf);

    void mfix_redistribute_deposition (int lev,
                                       amrex::MultiFab & mf_eps,
                                       amrex::MultiFab & mf_to_redistribute,
                                       const amrex::MultiFab * volfrac,
                                       const amrex::FabArray<EBCellFlagFab>* flags,
                                       amrex::Real max_eps);

    Real volWgtSum (int lev, const MultiFab & mf, int comp, bool local=false);

    Real volEpsWgtSum (int lev, const MultiFab & mf, int comp, bool local=false);

    void mfix_set_bc_type (int lev);

    ////////////////////////////////////////////////////////////////////////////
    //
    // Convection
    //
    ////////////////////////////////////////////////////////////////////////////

    void set_MAC_velocity_bcs (int lev,
                               Vector< std::unique_ptr<MultiFab> >& u,
                               Vector< std::unique_ptr<MultiFab> >& v,
                               Vector< std::unique_ptr<MultiFab> >& w,
                               Real time);

    void mfix_compute_fluxes (const int lev,
                              Vector< std::unique_ptr<MultiFab> >& a_fx,
                              Vector< std::unique_ptr<MultiFab> >& a_fy,
                              Vector< std::unique_ptr<MultiFab> >& a_fz,
                              Vector< std::unique_ptr<MultiFab> >& state,
                              const int state_comp, const int n_comp,
                              Vector< std::unique_ptr<MultiFab> >& xslopes,
                              Vector< std::unique_ptr<MultiFab> >& yslopes,
                              Vector< std::unique_ptr<MultiFab> >& zslopes,
                              const int slopes_comp,
                              Vector< std::unique_ptr<MultiFab> >& ep_u_mac,
                              Vector< std::unique_ptr<MultiFab> >& ep_v_mac,
                              Vector< std::unique_ptr<MultiFab> >& ep_w_mac);

    void mfix_compute_fluxes_on_box (const int lev, Box& bx,
                                     FArrayBox& a_fx,
                                     FArrayBox& a_fy,
                                     FArrayBox& a_fz,
                                     const FArrayBox& state,
                                     const int state_comp, const int n_comp,
                                     const FArrayBox& xslopes,
                                     const FArrayBox& yslopes,
                                     const FArrayBox& zslopes,
                                     const int slopes_comp,
                                     const FArrayBox& ep_u_mac,
                                     const FArrayBox& ep_v_mac,
                                     const FArrayBox& ep_w_mac);

    void mfix_compute_eb_fluxes_on_box (const int lev, Box& bx,
                                        FArrayBox& a_fx,
                                        FArrayBox& a_fy,
                                        FArrayBox& a_fz,
                                        const FArrayBox& state,
                                        const int state_comp, const int n_comp,
                                        const FArrayBox& xslopes,
                                        const FArrayBox& yslopes,
                                        const FArrayBox& zslopes,
                                        const int slopes_comp,
                                        const FArrayBox& ep_u_mac,
                                        const FArrayBox& ep_v_mac,
                                        const FArrayBox& ep_w_mac,
                                        const FArrayBox& ax,
                                        const FArrayBox& ay,
                                        const FArrayBox& az,
                                        const FArrayBox& fc_x,
                                        const FArrayBox& fc_y,
                                        const FArrayBox& fc_z,
                                        const FArrayBox& volfrac,
                                        const FArrayBox& bndrycent,
                                        const IArrayBox& cc_mask,
                                        const EBCellFlagFab& flags);

    ////////////////////////////////////////////////////////////////////////////
    //
    // Diffusion
    //
    ////////////////////////////////////////////////////////////////////////////

    std::unique_ptr<DiffusionOp> diffusion_op;

    ////////////////////////////////////////////////////////////////////////////
    //
    // Nodal Projection
    //
    ////////////////////////////////////////////////////////////////////////////

    std::unique_ptr<NodalProjector> nodal_projector;

  protected:
    //! Tagging cells for refinement
    virtual void ErrorEst (int lev,
                           TagBoxArray & tags,
                           Real time,
                           int ngrow) override;

    //! Make a new level using provided BoxArray and DistributionMapping and
    //! fill with interpolated coarse level data. Called by AmrCore::regrid.
    virtual void MakeNewLevelFromCoarse (int lev, Real time,
                                         const BoxArray & ba,
                                         const DistributionMapping & dm) override
    { amrex::Abort("MakeNewLevelFromCoarse: To be implemented"); }

    //! Remake an existing level using provided BoxArray and
    //! DistributionMapping and fill with existing fine and coarse data. Called
    //! by AmrCore::regrid.
    virtual void RemakeLevel (int lev,
                              Real time,
                              const BoxArray & ba,
                              const DistributionMapping & dm) override
    { amrex::Abort("RemakeLevel: To be implemented"); }

    //! Delete level data.  Called by AmrCore::regrid.
    virtual void ClearLevel (int lev) override
    { amrex::Abort("ClearLevel: To be implemented"); }

    void mfix_init_fluid (int is_restarting, Real dt, Real stop_time);

    void mfix_set_bc0 ();

    void mfix_set_p0 ();

    void mfix_set_gp0 ();

    static RealVect gp0;
    static RealVect gravity;

    // Scalar diffusive coefficient
    Vector<Real> mu_s;

    BcList bc_list;

    // "mutable" keyword is just a temporary fix to allow NodalProjection to
    // work as intended. This won't be needed when the BC will be given their own
    // class
    mutable amrex::Gpu::ManagedVector<amrex::Real> m_bc_u_g;
    mutable amrex::Gpu::ManagedVector<amrex::Real> m_bc_v_g;
    mutable amrex::Gpu::ManagedVector<amrex::Real> m_bc_w_g;
    mutable amrex::Gpu::ManagedVector<amrex::Real> m_bc_t_g;
    mutable amrex::Gpu::ManagedVector<amrex::Real> m_bc_ep_g;
    mutable amrex::Gpu::ManagedVector<amrex::Real> m_bc_p_g;

    //! The particles see an inflow face as a solid wall; thus instead of
    //! creating a separate EBfactory for particles, we simply modify the
    //! level-set near inflow to see it as a wall
    void mfix_set_ls_near_inflow ();

    void mfix_project_velocity ();

    void mfix_initial_iterations (Real dt, Real stop_time);

    void mfix_apply_nodal_projection (Vector< std::unique_ptr< amrex::MultiFab > >& a_depdt,
                                      amrex::Real a_time,
                                      amrex::Real a_dt,
                                      bool proj_2 );

    void apply_MAC_projection (Vector< std::unique_ptr<MultiFab> >& ep_u_mac,
                               Vector< std::unique_ptr<MultiFab> >& ep_v_mac,
                               Vector< std::unique_ptr<MultiFab> >& ep_w_mac,
                               Vector< std::unique_ptr<MultiFab> >& ep,
                               Vector< std::unique_ptr<MultiFab> >& ro,
                               Real time);

    void FillPatchVel (int lev, Real time, MultiFab & mf, int icomp, int ncomp,
                       const Vector<BCRec> & bcr);

    void GetDataVel (int lev, Real time, Vector<MultiFab *> & data,
                     Vector<Real> & datatime);

    void FillPatchScalar (int lev, Real time, MultiFab & mf, int icomp, int ncomp,
                          const Vector<BCRec> & bcr);

    void GetDataScalar (int lev, Real time, Vector<MultiFab *> & data, int icomp,
                        Vector<Real> & datatime);

    void mfix_calc_particle_beta (Real time);

    void mfix_calc_drag_fluid (Real time);

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

    void MakeBCArrays ();

    void MakeNewLevelFromScratch (int lev, Real time, const BoxArray & new_grids,
                                  const DistributionMapping & new_dmap) override;

    void ReMakeNewLevelFromScratch (int lev, const BoxArray & new_grids,
                                    const DistributionMapping & new_dmap);

    void WriteCheckHeader (const std::string & name, int nstep, Real dt, Real time) const;

    void WriteJobInfo (const std::string & dir ) const;

    static void GotoNextLine (std::istream & is);

    void EvolveFluid (int nstep, Real & dt, Real & time, Real stop_time, Real drag_timing);

    // Projection-related methods
    void mfix_solve_poisson_equation (Vector< std::unique_ptr<MultiFab> >& phi,
                                      Vector< std::unique_ptr<MultiFab> >& rhs,
                                      Vector< std::unique_ptr<MultiFab> >& b,
                                      Vector< std::unique_ptr<MultiFab> >& fluxes);

    // Solve implicit diffusion
    void solve_diffusion_equation (Vector< Vector< std::unique_ptr<MultiFab> > >& b,
                                   Vector< std::unique_ptr<MultiFab> >& phi,
                                   Vector< std::unique_ptr<MultiFab> >& rhs,
                                   int bc_lo[], int bc_hi[],
                                   Real dt);

    void check_for_nans (int lev);

    void mfix_compute_convective_term (Vector< std::unique_ptr<MultiFab> >& conv_u,
                                       Vector< std::unique_ptr<MultiFab> >& conv_s,
                                       Vector< std::unique_ptr<MultiFab> >& vel,
                                       Vector< std::unique_ptr<MultiFab> >& ep_g,
                                       Vector< std::unique_ptr<MultiFab> >& ro_g,
                                       Vector< std::unique_ptr<MultiFab> >& trac,
                                       Real time);

    int steady_state_reached (Real dt, int iter);

    void mfix_apply_predictor (Vector< std::unique_ptr<MultiFab> >& conv_u_old,
                               Vector< std::unique_ptr<MultiFab> >& conv_s_old,
                               Vector< std::unique_ptr<MultiFab> >& divtau_old,
                               Vector< std::unique_ptr<MultiFab> >&   laps_old,
                               Real time, Real dt, bool proj_2);

    void mfix_apply_corrector (Vector< std::unique_ptr<MultiFab> >& conv_u_old,
                               Vector< std::unique_ptr<MultiFab> >& conv_s_old,
                               Vector< std::unique_ptr<MultiFab> >& divtau_old,
                               Vector< std::unique_ptr<MultiFab> >&   laps_old,
                               Real time, Real dt, bool proj_2);

    void mfix_diffuse_scalar (const amrex::Vector< std::unique_ptr<MultiFab> >& mf_eps, amrex::Real);

    void mfix_diffuse_array (const amrex::Vector< std::unique_ptr<MultiFab> >& mf_drag, amrex::Real);

    void mfix_compute_diveu (Real time);

    void mfix_print_max_vel (int lev);

    void mfix_print_max_gp (int lev);

    // Make a level 0 grids covering the whole domain.  It does NOT install the new grids.
    BoxArray MakeBaseGrids () const;

    void ChopGrids (const Box& domain, BoxArray& ba, int target_size) const;

    static int nlev;
    static int ntrac;

    // Particle container
    std::unique_ptr<MFIXParticleContainer> pc;

    // Debug counter used for tracking number of level-set MultiFabs used in
    // intersection/union
    int ct_ls_mf;

    // Boundary conditions types
    Vector< std::unique_ptr<IArrayBox> > bc_ilo;
    Vector< std::unique_ptr<IArrayBox> > bc_ihi;
    Vector< std::unique_ptr<IArrayBox> > bc_jlo;
    Vector< std::unique_ptr<IArrayBox> > bc_jhi;
    Vector< std::unique_ptr<IArrayBox> > bc_klo;
    Vector< std::unique_ptr<IArrayBox> > bc_khi;

    // Boundary conditions flags
    Vector< std::unique_ptr<iMultiFab> > flag;

    // Void fraction
    Vector< std::unique_ptr<MultiFab> > ep_g;
    Vector< MultiFab* > ep_go;

    // Gas pressure fraction
    Vector< std::unique_ptr<MultiFab> > p_g;
    Vector< std::unique_ptr<MultiFab> > p_go;

    // Gas density
    Vector< std::unique_ptr<MultiFab> > ro_g;
    Vector< std::unique_ptr<MultiFab> > ro_go;

    // Tracer in gas
    Vector< std::unique_ptr<MultiFab> > trac;
    Vector< std::unique_ptr<MultiFab> > trac_o;

    // Gas velocity
    Vector< std::unique_ptr<MultiFab> > vel_g;
    Vector< std::unique_ptr<MultiFab> > vel_go;

    // Base state pressure
    Vector< std::unique_ptr<MultiFab> > p0_g;

    // Pressure gradients
    Vector< std::unique_ptr<MultiFab> > gp;

    // Molecular viscosity
    Vector< std::unique_ptr<MultiFab> > mu_g;

    // Cell-based
    Vector< std::unique_ptr<MultiFab> > vort;

    Vector< std::unique_ptr<MultiFab> > drag;

    // Level-Set Data => used for particle-wall collisions and fluid
    // reconstruction in particle drag calculation. NOTE: this has at least 2
    // levels: 0) fluid reconstruction, 1..N) particle-wall (refined) collisions
    Vector< std::unique_ptr<MultiFab> > level_sets;

    // These are multi-component multifabs
    Vector< std::unique_ptr<MultiFab> > xslopes_u;
    Vector< std::unique_ptr<MultiFab> > yslopes_u;
    Vector< std::unique_ptr<MultiFab> > zslopes_u;
    Vector< std::unique_ptr<MultiFab> > xslopes_s;
    Vector< std::unique_ptr<MultiFab> > yslopes_s;
    Vector< std::unique_ptr<MultiFab> > zslopes_s;

    // div (ep_g * u)
    Vector< std::unique_ptr<MultiFab> > diveu;

    // RHS for MAC solve
    Vector< std::unique_ptr<MultiFab> > mac_rhs;

    // Solution for MAC projection
    Vector< std::unique_ptr<MultiFab> > mac_phi;

    // RHS for diffusive tensor solve
    Vector< std::unique_ptr<MultiFab> > diff_rhs;
    Vector< std::unique_ptr<MultiFab> > diff_rhs1;
    Vector< std::unique_ptr<MultiFab> > diff_rhs4;

    // Solution for diffusion solves
    Vector< std::unique_ptr<MultiFab> > diff_phi;
    Vector< std::unique_ptr<MultiFab> > diff_phi1;
    Vector< std::unique_ptr<MultiFab> > diff_phi4;

    // MAC velocities
    Vector< std::unique_ptr<MultiFab> > u_mac;
    Vector< std::unique_ptr<MultiFab> > v_mac;
    Vector< std::unique_ptr<MultiFab> > w_mac;

    // Pressure increment
    Vector< std::unique_ptr<MultiFab> > phi_nd;

    // used if load_balance_type == "KnapSack"
    Vector<std::unique_ptr<MultiFab> > particle_cost;
    Vector<std::unique_ptr<MultiFab> > fluid_cost;

    int steady_state = 0;

    bool advect_density = false;
    bool advect_tracer  = false;

    bool test_tracer_conservation  = false;

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
    int  nodal_mg_max_coarsening_level = 100;

    // Linear operator
    std::unique_ptr<MLNodeLaplacian> nodal_matrix;

    // Solver
    std::unique_ptr<MLMG> nodal_solver;

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

    // Control MLMG behavior
    Real mac_mg_rtol = 1.0e-11;
    Real mac_mg_atol = 1.0e-14;

    // Max coarsening level
    int mac_mg_max_coarsening_level = 100;

    // ***************************************************************
    // Diffusion solve
    // ***************************************************************

    // Verbosity
    int diff_mg_verbose = 0;
    int diff_mg_cg_verbose = 0;

    // Control MLMG behavior
    int  diff_mg_maxiter = 100;
    int  diff_mg_cg_maxiter = 100;
    Real diff_mg_rtol = 1.0e-11;
    Real diff_mg_atol = 1.0e-14;
    std::string diff_bottom_solver_type;

    // Max coarsening level
    int  diff_mg_max_coarsening_level = 100;

    // ***************************************************************

    // Tolerance to check for steady state --
    //    this must be specified in the inputs file
    Real steady_state_tol;

    // Maximum number of iterations to steady state --
    //    this default may be over-written in the inputs file
    int steady_state_maxiter = 100000000;

    int call_udf;

    bool dual_grid = false;

    // If true then print the name of the routine we are in
    bool ooo_debug = false;

    // Parameters to control load balancing
    static std::string particle_init_type;
    static std::string load_balance_type;
    static std::string knapsack_weight_type;
    static int load_balance_fluid;
    static int knapsack_nmax;

    // Options to control time stepping
    Real cfl = 0.5;
    Real fixed_dt;
    Real dt_min  = 0.0;           // Min dt allowed
    Real dt_max  = 1.e14;         // Max dt allowed

    Real prev_dt = 0.0;           // dt used in previous timestep

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

    amrex::Vector<Real> avg_vel_p;

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
    int plt_trac    = 0;
    int plt_mu_g    = 0;
    int plt_diveu   = 0;
    int plt_volfrac = 0;
    int plt_gradp_g = 0;
    int plt_vort    = 0;

    // Total number of variables to write in plot file
    int pltVarCount = 0;

    // Flags for saving particle data. By default, we have all flags on,
    // we turn off what we don't want in the init IO routine. This is somewhat
    // different from what we do with the fluid.
    Vector<int> write_real_comp = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    Vector<int> write_int_comp = {1, 1};

    //
    // Here we set the number of ghost nodes of the field arrays.
    // Since this is  a class attribute, it will be propagated to
    // all the member functions.
    //
    // DO NOT USE  A LOCAL VARIABLE WITH THE SAME NAME
    //

    // Four ghost cells are required when using EB
    const int nghost = 4;

    const int m_eb_basic_grow_cells = nghost;
    const int m_eb_volume_grow_cells = nghost;
    const int m_eb_full_grow_cells = nghost;

    bool do_initial_proj    = true;
    int  initial_iterations = 3;

    /***************************************************************************
     *                                                                         *
     * I/O stuff                                                               *
     *                                                                         *
     **************************************************************************/

    // Variables to simplify checkpoint IO
    Vector< const Vector< std::unique_ptr<MultiFab> >* > vectorVars;
    Vector< std::string > vecVarsName;

    Vector< const Vector< std::unique_ptr<MultiFab> >* > chkscalarVars;
    Vector< std::string > chkscaVarsName;

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
    Vector<BCRec> bcs_f; // For generic first-order extrapolation

    /***************************************************************************
     *                                                                         *
     * EB Data (including level-set data)                                      *
     *                                                                         *
     **************************************************************************/

    //! EB levels representing fluid boundary conditions
    Vector<const EB2::Level *> eb_levels;

    //! EB levels representing particle boundary conditions (same as
    //! `mfix::eb_levels` but might include additional walls at MI BCs).
    Vector<const EB2::Level *> particle_eb_levels;

    //! EB factory that lives on the fluid grids
    Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > ebfactory;

    //! EB factory that lives on the particle grids
    Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > particle_ebfactory;

    bool mfix_update_ebfactory (int a_lev);
};

inline std::string mfix::get_load_balance_type ()
{ return load_balance_type; }

#endif
