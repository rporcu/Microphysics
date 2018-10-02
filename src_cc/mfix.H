#ifndef MFIX_LEVEL_H_
#define MFIX_LEVEL_H_

#include <iostream>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BLProfiler.H>

#include <MFIXParticleContainer.H>

#include <AMReX_EBFabFactory.H>

#include <AMReX_EBFArrayBox.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <mfix_eb_if.H>

#include <AMReX_EB_levelset.H>

#include <MacProjection.H>

class mfix
    : public AmrCore
{
public:

    mfix ();
    ~mfix ();

    void InitParams(int solve_fluid, int solve_dem, int call_udf );

    void Init(Real dt, Real time);

    void InitLevelData(Real dt, Real time);

    void PostInit(Real dt, Real time, int nstep, int restart_flag, Real stop_time,
                  int steady_state);

    void ResizeArrays();

    void WriteCheckPointFile(std::string& check_file_name,
            int nstep = 0, Real dt = 0.0, Real time = 0.0)  const;

    void WritePlotFile(std::string& plot_file_name,
            int nstep = 0, Real dt = 0.0, Real time = 0.0) const;

    void Restart(std::string& restart_chkfile, int *nstep, Real *dt, Real *time,
            IntVect& Nrep);

    void WriteParticleAscii(std::string& par_ascii_file_name, int nstep = 0) const;

    void WriteUSER(Real dt, Real time) const;

    void Regrid ();

    void Evolve(int nstep, int steady_state, Real& dt, Real & prev_dt, Real time, Real stop_time);

    void usr3();

    void output(int estatus, int finish, int nstep, Real dt, Real time)
    {
        pc -> output( estatus, finish, nstep, dt, time);
    };

    void make_mg_bc(int mg_bc[]);

    void make_eb_geometry ();

    std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_walls(bool & has_walls);

    std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_real_walls(bool & has_real_walls);

    std::unique_ptr<
        EB2::TranslationIF<EB2::PolynomialIF>
    > get_poly(int max_order, std::string field_prefix);

    std::unique_ptr<EB2::IntersectionIF<EB2::PlaneIF,EB2::PlaneIF,EB2::PlaneIF>> make_wall(
            int dir, // direction (long edge) of wall
            Real position, Real height, Real width
        );

    void WriteEBSurface();

    void InitIOData();

    static std::string get_load_balance_type();

    // flag enabling level-set restart (i.e. prevent make_eb_* from rebuilding
    // the level-set data).
    bool levelset__restart = false;

    // flag switching from calc_wall_collisions_ls to calc_wall_collisions if true
    bool legacy__eb_collisions = false;

    void mfix_compute_vort();

    ///
    /// Parameters describing the level-set grid:
    ///
    int levelset__refinement = 1;    // refinement (wrt particle grid) of the level-set's grid
    int levelset__eb_refinement = 1; // refinement of the EB facets used to generate the level-set
    int levelset__pad = 2;           // padding of the level-set grid
    int levelset__eb_pad = 2;        // padding of the EBIS used to generate the level-set (useful for
                                     // ensuring that neighbor grids are considered when filling LS)
protected:

    //! Tagging cells for refinement
    virtual void ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/) override
        { amrex::Abort("ErrorEst: To be implemented"); }

    //! Make a new level using provided BoxArray and
    //! DistributionMapping and fill with interpolated coarse level
    //! data.  Called by AmrCore::regrid.
    virtual void MakeNewLevelFromCoarse(int lev, amrex::Real time, const amrex::BoxArray& ba,
            const amrex::DistributionMapping& dm) override
    {
        amrex::Abort("MakeNewLevelFromCoarse: To be implemented");
    }

    //! Remake an existing level using provided BoxArray and
    //! DistributionMapping and fill with existing fine and coarse
    //! data.  Called by AmrCore::regrid.
    virtual void RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray& ba,
            const amrex::DistributionMapping& dm) override
    {
        amrex::Abort("RemakeLevel: To be implemented");
    }

    //! Delete level data.  Called by AmrCore::regrid.
    virtual void ClearLevel(int lev) override
    {
        amrex::Abort("ClearLevel: To be implemented");
    }

    void mfix_init_fluid(int is_restarting, Real dt, Real stop_time, int steady_state);
    void mfix_set_bc0();
    void mfix_set_p0();

    void mfix_project_velocity   ();
    void mfix_initial_iterations ( amrex::Real dt, amrex::Real stop_time, int steady_state );
    void mfix_apply_projection ( amrex::Real scaling_factor, bool proj_2 );

    void mfix_set_scalar_bcs ();
    void mfix_set_velocity_bcs (int extrap_dir_bcs);

    void mfix_calc_volume_fraction(Real& sum_vol);
    void mfix_calc_drag_fluid();
    void mfix_calc_drag_particle();

    void AllocateArrays (int lev);

    void RegridArrays (int lev);

private:

    void make_eb_general  ();
    void make_eb_box      ();
    void make_eb_cylinder ();
    void make_eb_hopper   ();
    void make_eb_cyclone  ();
    void make_eb_hourglass();
    void make_eb_clr      ();
    void make_eb_clr_riser();
    void make_eb_regular  ();

    void MakeBCArrays();

    void check_data();

    void MakeNewLevelFromScratch(int lev, Real time,
            const BoxArray& new_grids, const DistributionMapping& new_dmap) override;

    void ReMakeNewLevelFromScratch(int lev,
            const BoxArray& new_grids, const DistributionMapping& new_dmap);

    void WriteHeader(const std::string& name, int nstep, Real dt, Real time, bool is_checkpoint) const;
    void WritePlotHeader(const std::string& name, int nstep, Real dt, Real time) const;
    void WriteCheckHeader(const std::string& name, int nstep, Real dt, Real time) const;

    void WriteJobInfo(const std::string& dir ) const;

    static void GotoNextLine(std::istream& is);

    void EvolveFluid(int nstep, int steady_state, Real& dt, Real& time, Real stop_time);

    // Projection-related methods
    void solve_poisson_equation(Vector< Vector< std::unique_ptr<MultiFab> > >& b,
                                Vector< std::unique_ptr<MultiFab> >& phi,
                                Vector< std::unique_ptr<MultiFab> >& rhs,
                                Vector< std::unique_ptr<MultiFab> >& fluxes,
                                int bc_lo[], int bc_hi[]);

    // Solve implicit diffusion
    void solve_diffusion_equation(Vector< Vector< std::unique_ptr<MultiFab> > >& b,
                                  Vector< std::unique_ptr<MultiFab> >& phi,
                                  Vector< std::unique_ptr<MultiFab> >& rhs,
                                  int bc_lo[], int bc_hi[],
                                  amrex::Real dt);

    void mfix_compute_bcoeff_ppe();
    void mfix_compute_bcoeff_diff();

    void check_for_nans(int lev);

    void mfix_compute_dt(amrex::Real time, amrex::Real stop_time, int steady_state, amrex::Real& dt);

    void mfix_compute_ugradu_predictor(Vector< std::unique_ptr<MultiFab> >& conv,
                                       Vector< std::unique_ptr<MultiFab> >& vel);
    
    void mfix_compute_ugradu_corrector(Vector< std::unique_ptr<MultiFab> >& conv,
                                       Vector< std::unique_ptr<MultiFab> >& vel);
    
    void mfix_compute_divtau(int lev, amrex::MultiFab& divtau, Vector< std::unique_ptr<MultiFab> >& vel);

    int  steady_state_reached(Real dt);

    void mfix_apply_predictor(Vector< std::unique_ptr<MultiFab> >& conv_old,
                              Vector< std::unique_ptr<MultiFab> >& divtau_old,
                              amrex::Real dt, bool proj_2);
    void mfix_apply_corrector(Vector< std::unique_ptr<MultiFab> >& conv_old,
                              Vector< std::unique_ptr<MultiFab> >& divtau_old,
                              amrex::Real dt, bool proj_2);

    void mfix_diffuse_velocity(amrex::Real dt);

    void mfix_compute_intermediate_velocity (int lev, amrex::Real dt);

    void mfix_apply_forcing_terms    (int lev, amrex::Real dt, Vector< std::unique_ptr<MultiFab> >& vel);
    void mfix_compute_velocity_slopes(int lev, Vector< std::unique_ptr<MultiFab> >& vel);
    void mfix_compute_velocity_at_faces (int lev, Vector< std::unique_ptr<MultiFab> >& vel);

    void mfix_compute_diveu ();

    void mfix_extrap_pressure (int lev, std::unique_ptr<amrex::MultiFab>& p);

    void mfix_print_max_vel(int lev);

    // Make a level 0 grids covering the whole domain.  It does NOT install the new grids.
    BoxArray MakeBaseGrids() const;
    void ChopGrids(const Box& domain, BoxArray& ba, int target_size) const;

    void mfix_set_bc_type(int lev);

    void fill_mf_bc(int lev,  MultiFab& mf);
#if 0
    void WritePlotFile() const;

      // which step?
    Vector<int> nsubsteps;  // how many substeps on each level?
#endif

    Real mfix_norm1 ( const Vector< std::unique_ptr<MultiFab>>& mf, int lev, int comp );
    Real mfix_norm1 ( MultiFab& mf, int lev, int comp );
    Real mfix_norm0 ( const Vector< std::unique_ptr<MultiFab>>& mf, int lev, int comp );
    Real mfix_norm0 ( MultiFab& mf, int lev, int comp );

    void mfix_average_cc_to_fc ( int lev, const MultiFab& cc,
                                 Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& fc );

    static int nlev;
   
    Vector<int> istep = {1} ;

    // Particle container
    std::unique_ptr<MFIXParticleContainer> pc;

    // Unit vectors in Cartesian space
    static amrex::IntVect e_x;
    static amrex::IntVect e_y;
    static amrex::IntVect e_z;

    // Debug counter used for tracking number of level-set MultiFabs used in
    // intersection/union
    int ct_ls_mf;

    Vector< std::unique_ptr<IArrayBox> > bc_ilo ;
    Vector< std::unique_ptr<IArrayBox> > bc_ihi ;
    Vector< std::unique_ptr<IArrayBox> > bc_jlo ;
    Vector< std::unique_ptr<IArrayBox> > bc_jhi ;
    Vector< std::unique_ptr<IArrayBox> > bc_klo ;
    Vector< std::unique_ptr<IArrayBox> > bc_khi ;

    Vector< std::unique_ptr<iMultiFab> > flag;

    // Void fraction
    Vector< std::unique_ptr<MultiFab> > ep_g ;
    Vector< std::unique_ptr<MultiFab> > ep_go;

    // Gas pressure fraction
    Vector< std::unique_ptr<MultiFab> > p_g ;
    Vector< std::unique_ptr<MultiFab> > p_go;

    // Gas density
    Vector< std::unique_ptr<MultiFab> > ro_g ;
    Vector< std::unique_ptr<MultiFab> > ro_go;

    // Gas bulk density
    Vector< std::unique_ptr<MultiFab> > rop_g ;
    Vector< std::unique_ptr<MultiFab> > rop_go;

    // Gas velocity
    Vector< std::unique_ptr<MultiFab> > vel_g ;
    Vector< std::unique_ptr<MultiFab> > vel_go;

    // Base state pressure
    Vector< std::unique_ptr<MultiFab> > p0_g;

    // Pressure correction
    Vector< std::unique_ptr<MultiFab> > pp_g;

    // Pressure gradients
    Vector< std::unique_ptr<MultiFab> > gp ;
    Vector< std::unique_ptr<MultiFab> > gp0;

    // Molecular viscosity
    Vector< std::unique_ptr<MultiFab> > mu_g ;

    // Cell-based
    Vector< std::unique_ptr<MultiFab> > lambda_g;
    Vector< std::unique_ptr<MultiFab> > trD_g;
    Vector< std::unique_ptr<MultiFab> > vort;

    Vector< std::unique_ptr<MultiFab> > drag;
    Vector< std::unique_ptr<MultiFab> > f_gds;

    // Level-Set Data
    Vector< std::unique_ptr<MultiFab> > ls;

    // These are multi-component multifabs
    Vector< std::unique_ptr<MultiFab> > xslopes;
    Vector< std::unique_ptr<MultiFab> > yslopes;
    Vector< std::unique_ptr<MultiFab> > zslopes;

    // div (ep_g * u)
    Vector< std::unique_ptr<MultiFab> > diveu;

    // RHS for implicit diffusion solve
    Vector< std::unique_ptr<MultiFab> > rhs_diff;

    // Solution for implicit diffusion solve
    Vector< std::unique_ptr<MultiFab> > phi_diff;

    //
    Vector< std::unique_ptr<MultiFab> > fp;

    // Coefficients b of div(b*grad(p))
    Vector< Vector< std::unique_ptr<MultiFab> > > bcoeff;

    // Coefficients b in implicit diffusion solve
    Vector< Vector< std::unique_ptr<MultiFab> > > bcoeff_diff;

    // Pressure increment
    Vector< std::unique_ptr<MultiFab> > phi;

    // MAC velocity
    Vector< std::unique_ptr<MultiFab> > m_u_mac;
    Vector< std::unique_ptr<MultiFab> > m_v_mac;
    Vector< std::unique_ptr<MultiFab> > m_w_mac;
    // Vector< std::unique_ptr<MultiFab> > m_uo_mac;
    // Vector< std::unique_ptr<MultiFab> > m_vo_mac;
    // Vector< std::unique_ptr<MultiFab> > m_wo_mac;

    // MAC projection object
    std::unique_ptr<MacProjection> mac_projection;

    // What solver to use as the bottom solver in the MLMG solves.
    std::string bottom_solver_type;

    // used if load_balance_type == "KnapSack"
    Vector<std::unique_ptr<MultiFab> > particle_cost;
    Vector<std::unique_ptr<MultiFab> > fluid_cost;

    int solve_fluid;
    int solve_dem;

    // Options to control MLMG behavior
    int  mg_verbose = 0;
    int  mg_cg_verbose = 0;
    int  mg_max_iter = 100;
    int  mg_cg_maxiter = 100;
    int  mg_max_fmg_iter = 0;
    Real mg_rtol = 1.0e-11;
    Real mg_atol = 1.0e-14;

    // Tolerance to check for steady state
    Real steady_state_tol = 1.0e-5  ;

    int call_udf;
    bool dual_grid = false;

    // enable output_manager during des_time_loop
    bool subdt_io;

    static std::string particle_init_type;
    static std::string load_balance_type;
    static std::string knapsack_weight_type;

    // Options to control time stepping
    Real cfl = 0.5;
    int fixed_dt = 0;
    Real dt_min=0.0;            // Min dt allowed
    Real dt_max=1.e14;          // Max dt allowed

    int particle_max_grid_size_x = -1;
    int particle_max_grid_size_y = -1;
    int particle_max_grid_size_z = -1;

    static EBSupport m_eb_support_level;

    Periodicity p0_periodicity;

    Real sum_vol_orig;

    //
    // Options to control the computing of particle eulerian velocities
    //
    int avg_vel_int = -1;
    std::string avg_vel_file {"avg_region"};
    vector<Real> avg_region_x_e;
    vector<Real> avg_region_x_w;
    vector<Real> avg_region_y_s;
    vector<Real> avg_region_y_n;
    vector<Real> avg_region_z_b;
    vector<Real> avg_region_z_t;

    //
    // Here we set the number of ghost nodes of the field arrays.
    // Since this is  a class attribute, it will be propagated to
    // all the member functions.
    //
    // DO NOT USE  A LOCAL VARIABLE WITH THE SAME NAME
    //

    // These values are required when fluid sees EB -- for now
    const int nghost = 5;

    const int m_eb_basic_grow_cells = nghost;
    const int m_eb_volume_grow_cells = nghost;
    const int m_eb_full_grow_cells = nghost;
 

    int nodal_pressure = 0;
    int explicit_diffusion = 1;

    int ugradu_type = 3;

    // Variables to simplify IO -- the vectors are the same for plt/chk
    Vector< const Vector< std::unique_ptr<MultiFab> >* > vectorVars;
    Vector< std::string > vecVarsName;

    Vector< const Vector< std::unique_ptr<MultiFab> >* > chkscalarVars;
    Vector< std::string > pltscaVarsName;

    Vector< const Vector< std::unique_ptr<MultiFab> >* > pltscalarVars;
    Vector< std::string > chkscaVarsName;


    /****************************************************************************
     *                                                                          *
     * EB Data (including level-set data)                                       *
     *                                                                          *
     ****************************************************************************/

    const EB2::Level * eb_level_fluid;
    const EB2::Level * eb_level_particles;

    Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > ebfactory;
    Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > particle_ebfactory;

    // stores eb_normals (commond index-space with cells)
    std::unique_ptr<MultiFab> eb_normals;
    std::unique_ptr<MultiFab> dummy;

    // stores level-set (nodal index space corresponding to cells, but n_refinement times finer)
    std::unique_ptr<LSFactory> level_set;


    void mfix_update_ebfactory (int a_lev);
   
};

inline
std::string
mfix::get_load_balance_type()
{
    return load_balance_type;
}

#endif
