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

#include <AMReX_GeometryShop.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBTower.H>
#include <AMReX_EBFArrayBox.H>

#include <AMReX_GeometryShop.H>
#include <AMReX_PlaneIF.H>
#include <AMReX_AllRegularService.H>
#include <AMReX_FlatPlateGeom.H>
#include <AMReX_EBISLayout.H>
#include <AMReX_EBGraph.H>
#include <AMReX_EBDebugOut.H>
#include <AMReX_EBCellFAB.H>
#include <AMReX_EBCellFactory.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_UnionIF.H>
#include <AMReX_TransformIF.H>
#include <AMReX_ComplementIF.H>
#include <AMReX_IntersectionIF.H>
#include <AMReX_LatheIF.H>
#include <AMReX_PolynomialIF.H>

#include <eb_levelset.H>

class mfix_level
    : public AmrCore
{
public:

    mfix_level ();
    ~mfix_level ();

    void InitParams(int solve_fluid, int solve_dem, int max_nit, int call_udf );

    void Init(int lev, Real dt, Real time);

    void InitLevelData(int lev, Real dt, Real time);

    void PostInit(int lev, Real dt, Real time, int nstep, int restart_flag, Real stop_time,
                  int steady_state );

    void ResizeArrays();

    void WriteCheckPointFile(std::string& check_file_name,
            int nstep = 0, Real dt = 0.0, Real time = 0.0)  const;

    void WritePlotFile(std::string& plot_file_name,
            int nstep = 0, Real dt = 0.0, Real time = 0.0) const;

    void Restart(std::string& restart_chkfile, int *nstep, Real *dt, Real *time,
            IntVect& Nrep);

    void WriteParticleAscii(std::string& par_ascii_file_name, int nstep = 0) const;

    void WriteUSER(int lev, Real dt, Real time) const;

    void Regrid(int lev, int nstep);

    void RegridOnRestart(int lev);

    void Evolve(int lev, int nstep, int set_normg, int steady_state, Real& dt, Real & prev_dt,
                Real time, Real stop_time, Real normg);

    void usr3(int lev);

    void output(int lev, int estatus, int finish, int nstep, Real dt, Real time)
    {
        pc -> output( lev, estatus, finish, nstep, dt, time);
    };

    void make_mg_bc(int mg_bc[]);

    void make_eb_geometry(int lev);
    void make_eb_hourglass(int lev);
    void make_eb_clr(int lev);
    void make_eb_clr_riser(int lev);

    std::unique_ptr<BaseIF> get_walls(int lev, bool anisotropic, bool & has_walls);
    std::unique_ptr<BaseIF> get_poly(int lev, int max_order, std::string field_prefix);

    std::unique_ptr<BaseIF> make_eb_wall(
            int dir, // direction (long edge) of wall
            Real position, Real height, Real width,
            int lev, int water_tight
        );

    std::unique_ptr<BaseIF> make_cylinder(
            int dir, // direction (orientation) of cylinder axis
            Real radius, Real length, const RealVect & translation,
            int lev, bool water_tight
        );

    std::unique_ptr<BaseIF> make_cone(int dir,
            Real radius1, Real radius2, Real length, const RealVect & translation,
            int lev, bool water_tight);

    void fill_levelset(int lev, bool use_walls, bool use_poly,
                       const BaseIF & impfunc_walls, const BaseIF & impfunc_poly,
                       int max_level, int grid_size, bool eb_verbosity);

    void WriteEBSurface(int lev);

    void InitIOData();

    // flag enabling level-set restart (i.e. prevent make_eb_* from rebuilding
    // the level-set data).
    bool levelset__restart = false;

    // flag switching from calc_wall_collisions_ls to calc_wall_collisions if true
    bool legacy__eb_collisions = false;

    ///
    /// Parameters describing the level-set grid:
    ///
    int levelset__refinement = 1;     // refinement (wrt particle grid) of the level-set's grid
    int levelset__eb_refinement = 1;  // refinement (wrt particle grid) of the EB facets used to generate the level-set
    int levelset__pad = 2;            // padding of the level-set grid
    int levelset__eb_pad = 2;         // padding of the EBIS used to generate the level-set (particularly useful for
                                      // ensuring that edge-cased occur well outside the physical domain)


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

    void mfix_calc_trd_and_tau(int lev);

    void mfix_init_fluid(int lev, int is_restarting, Real dt, Real stop_time, int steady_state);
    void mfix_set_bc0(int lev);

    void mfix_conv_rop(int lev);

    void mfix_calc_mflux(int lev);

    void mfix_solve_for_u(int lev, Real dt, Real& num_u, Real& denom_u);
    void mfix_solve_for_v(int lev, Real dt, Real& num_v, Real& denom_v);
    void mfix_solve_for_w(int lev, Real dt, Real& num_w, Real& denom_w);

    // Projection
    void mfix_project_velocity   (int lev);
    void mfix_initial_iterations (int lev, Real dt, Real stop_time, int steady_state);
    void mfix_apply_projection (int lev, amrex::Real scaling_factor);
    void mfix_set_projection_bcs (int lev);
    void mfix_set_velocity_bcs (int lev);

    // End projection

    void mfix_solve_linear_equation(int eq_id,int lev,MultiFab& sol, MultiFab& matrix, MultiFab& rhs);

    void mfix_solve_for_pp(int lev, Real dt, Real& num_p, Real& denom_p);

    void mfix_correct_0(int lev);

    int solve_bicgstab(MultiFab& sol, const MultiFab& rhs, const MultiFab& matrix,
            int sweep_type, int precond_type, int maxiter, Real eps_rel, int lev);

    void mfix_calc_volume_fraction(int lev, Real& sum_vol);
    void mfix_calc_drag_fluid(int lev);
    void mfix_calc_drag_particle(int lev);

    void AllocateArrays(int lev);
    void RegridArrays(int lev, BoxArray& new_grids, DistributionMapping& new_dmap);

    void MakeBCArrays();

private:

    void check_data(int lev);

    void MakeNewLevelFromScratch(int lev, Real time,
            const BoxArray& new_grids, const DistributionMapping& new_dmap) override;

    void ReMakeNewLevelFromScratch(int lev,
            const BoxArray& new_grids, const DistributionMapping& new_dmap);

    void WriteHeader(const std::string& name, int nstep, Real dt, Real time, bool is_checkpoint) const;
    void WritePlotHeader(const std::string& name, int nstep, Real dt, Real time) const;
    void WriteCheckHeader(const std::string& name, int nstep, Real dt, Real time) const;

    void WriteJobInfo(const std::string& dir ) const;

    static void GotoNextLine(std::istream& is);

    void EvolveFluidSimple    (int lev, int nstep, int set_normg, Real dt, Real& prev_dt, Real time, Real normg);
    void EvolveFluidProjection(int lev, int nstep, int steady_state, Real& dt, Real stop_time, Real time);

    // Projection-related methods
    void solve_poisson_equation(int lev,
            Vector< Vector< std::unique_ptr<MultiFab> > >& b,
      Vector< std::unique_ptr<MultiFab> >& phi,
      Vector< std::unique_ptr<MultiFab> >& rhs,
      int bc_lo[], int bc_hi[] );

    void mfix_compute_bcoeff(int lev);

    void check_for_nans(int lev);

    void mfix_compute_dt(int lev, amrex::Real time, amrex::Real stop_time, int steady_state, amrex::Real& dt);

    void mfix_compute_fluid_acceleration(int lev,
            Vector< std::unique_ptr<MultiFab> >& u,
            Vector< std::unique_ptr<MultiFab> >& v,
            Vector< std::unique_ptr<MultiFab> >& w );

    int  steady_state_reached(int lev, Real dt);

    void mfix_apply_predictor (int lev, amrex::Real dt);

    void mfix_apply_corrector (int lev, amrex::Real dt);


    void mfix_apply_forcing_terms(int lev, amrex::Real dt,
            Vector< std::unique_ptr<MultiFab> >& u,
            Vector< std::unique_ptr<MultiFab> >& v,
            Vector< std::unique_ptr<MultiFab> >& w );

    void mfix_compute_intermediate_velocity (int lev, amrex::Real dt);

    void mfix_add_pressure_gradient(int lev, amrex::Real coeff);

    void mfix_add_grad_phi (int lev, amrex::Real coeff);

    void mfix_compute_velocity_slopes(int lev,
            Vector< std::unique_ptr<MultiFab> >& u,
            Vector< std::unique_ptr<MultiFab> >& v,
            Vector< std::unique_ptr<MultiFab> >& w );

    void mfix_compute_diveu (int lev);

    void mfix_extrap_pressure (int lev, std::unique_ptr<amrex::MultiFab>& p);

    void mfix_print_max_vel(int lev);
// End of projection-related methods


    // Make a level 0 grids covering the whole domain.  It does NOT install the new grids.
    BoxArray MakeBaseGrids() const;
    void ChopGrids(const Box& domain, BoxArray& ba, int target_size) const;

    void mfix_set_bc_type(int lev);
    void mfix_set_bc1(int lev);

    void fill_mf_bc(int lev,  MultiFab& mf);
#if 0
    void WritePlotFile() const;

      // which step?
    Vector<int> nsubsteps;  // how many substeps on each level?
#endif

    Vector<int> istep = {1} ;

    // Particle container
    std::unique_ptr<MFIXParticleContainer> pc;

    // Unit vectors in Cartesian space
    static amrex::IntVect e_x;
    static amrex::IntVect e_y;
    static amrex::IntVect e_z;

    // Debug counter used for tracking number of leve-lset MultiFabs used in
    // intersection/union
    int ct_ls_mf;

//  Array< std::unique_ptr<iMultiFab> > bc_type;
//  Array< std::unique_ptr<iMultiFab> > bc_ptr;
    IArrayBox bc_ilo;
    IArrayBox bc_ihi;
    IArrayBox bc_jlo;
    IArrayBox bc_jhi;
    IArrayBox bc_klo;
    IArrayBox bc_khi;

    Vector< std::unique_ptr<iMultiFab> > flag;

    // Matrix and rhs vector
    Vector< std::unique_ptr<MultiFab> > A_m;
    Vector< std::unique_ptr<MultiFab> > b_m;

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

    // X-axis gas velocity
    Vector< std::unique_ptr<MultiFab> > u_g ;
    Vector< std::unique_ptr<MultiFab> > u_go;
    Vector< std::unique_ptr<MultiFab> > u_gt;

    // Y-axis gas velocity
    Vector< std::unique_ptr<MultiFab> > v_g ;
    Vector< std::unique_ptr<MultiFab> > v_go;
    Vector< std::unique_ptr<MultiFab> > v_gt;

    // Z-axis gas velocity
    Vector< std::unique_ptr<MultiFab> > w_g ;
    Vector< std::unique_ptr<MultiFab> > w_go;
    Vector< std::unique_ptr<MultiFab> > w_gt;

    // Base state pressure
    Vector< std::unique_ptr<MultiFab> > p0_g;

    // Pressure correction
    Vector< std::unique_ptr<MultiFab> > pp_g;
    Vector< std::unique_ptr<MultiFab> > d_e ;
    Vector< std::unique_ptr<MultiFab> > d_n ;
    Vector< std::unique_ptr<MultiFab> > d_t ;

    // Molecular viscosity
    Vector< std::unique_ptr<MultiFab> > mu_g ;

    // Level-Set Data
    Vector< std::unique_ptr<MultiFab> > ls;

    // Cell-based
    Vector< std::unique_ptr<MultiFab> > lambda_g;
    Vector< std::unique_ptr<MultiFab> > trD_g;
    Vector< std::unique_ptr<MultiFab> > vort;

    // Face-based
    Vector< std::unique_ptr<MultiFab> > tau_u_g;
    Vector< std::unique_ptr<MultiFab> > tau_v_g;
    Vector< std::unique_ptr<MultiFab> > tau_w_g;

    Vector< std::unique_ptr<MultiFab> > fluxX;
    Vector< std::unique_ptr<MultiFab> > fluxY;
    Vector< std::unique_ptr<MultiFab> > fluxZ;

    Vector< std::unique_ptr<MultiFab> > ropX;
    Vector< std::unique_ptr<MultiFab> > ropY;
    Vector< std::unique_ptr<MultiFab> > ropZ;

    Vector< std::unique_ptr<MultiFab> > f_gds_u;
    Vector< std::unique_ptr<MultiFab> > f_gds_v;
    Vector< std::unique_ptr<MultiFab> > f_gds_w;

    Vector< std::unique_ptr<MultiFab> > drag_u;
    Vector< std::unique_ptr<MultiFab> > drag_v;
    Vector< std::unique_ptr<MultiFab> > drag_w;

    //
    // Projection
    //

    // Fluid acceleration (convection + diffusion)
    Vector< std::unique_ptr<MultiFab> > uacc;
    Vector< std::unique_ptr<MultiFab> > vacc;
    Vector< std::unique_ptr<MultiFab> > wacc;

    // These are multi-component multifab
    Vector< std::unique_ptr<MultiFab> > slopes_u;
    Vector< std::unique_ptr<MultiFab> > slopes_v;
    Vector< std::unique_ptr<MultiFab> > slopes_w;

    // div (ep_g * u)
    Vector< std::unique_ptr<MultiFab> > diveu;

    //
    Vector< std::unique_ptr<MultiFab> > fp_x;
    Vector< std::unique_ptr<MultiFab> > fp_y;
    Vector< std::unique_ptr<MultiFab> > fp_z;

    Vector< Vector< std::unique_ptr<MultiFab> > > bcoeff; // Coefficients b of div(b*grad(p))

    Vector< std::unique_ptr<MultiFab> > phi; // Pressure increment

    // used if load_balance_type == "KnapSack"
    Vector<std::unique_ptr<MultiFab> > particle_cost;
    Vector<std::unique_ptr<MultiFab> > fluid_cost;

    int solve_fluid;
    int solve_dem;
    int use_proj_method = 0; // Whether to use projection to solve NS eqs.
    // Options to control MLMG behavior (projection only)
    int  mg_verbose = 0;
    int  mg_cg_verbose = 0;
    int  mg_max_iter = 100;
    int  mg_cg_maxiter = 100;
    int  mg_max_fmg_iter = 0;
    Real mg_rtol = 1.0e-11;
    Real mg_atol = 1.0e-14;

    // Tolerance to check for steady state (projection only)
    Real steady_state_tol = 1.0e-5  ;

    int max_nit;
    int call_udf;
    int dual_grid = 0;

    // enable output_manager during des_time_loop
    bool subdt_io;

    static std::string particle_init_type;
    static std::string load_balance_type;
    static std::string knapsack_weight_type;

    // Options to control time stepping
    Real cfl = 0.5;
    int  fixed_dt=0;

    int particle_max_grid_size_x = -1;
    int particle_max_grid_size_y = -1;
    int particle_max_grid_size_z = -1;

    static int m_eb_basic_grow_cells;
    static int m_eb_volume_grow_cells;
    static int m_eb_full_grow_cells;
    static EBSupport m_eb_support_level;

    Periodicity p0_periodicity;

    Real sum_vol_orig;

    //
    // Here we set the number of ghost nodes of the field arrays.
    // Since this is  a class attribute, it will be propagated to
    // all the member functions.
    //
    // DO NOT USE  A LOCAL VARIABLE WITH THE SAME NAME
    //
    const int nghost = 1;

    // Variables to simplify IO -- the vectors are the same for plt/chk
    Vector< const Vector< std::unique_ptr<MultiFab> >* > vectorVars;
    Vector< std::string > vecVarsName;

    Vector< const Vector< std::unique_ptr<MultiFab> >* > chkscalarVars;
    Vector< std::string > pltscaVarsName;

    Vector< const Vector< std::unique_ptr<MultiFab> >* > pltscalarVars;
    Vector< std::string > chkscaVarsName;

    std::unique_ptr<amrex::EBFArrayBoxFactory> ebfactory;
    std::unique_ptr<amrex::EBFArrayBoxFactory> particle_ebfactory;

    // stores eb_normals (commond index-space with cells)
    std::unique_ptr<MultiFab> eb_normals;
    std::unique_ptr<MultiFab> dummy;

    // stores level-set (nodal index space corresponding to cells, but n_refinement times finer)
    std::unique_ptr<LSFactory> level_set;
};

#endif
