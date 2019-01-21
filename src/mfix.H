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
#include <AMReX_PhysBCFunct.H>

#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB_levelset.H>
#include <AMReX_EB_LSCore.H>

#include <MFIXParticleContainer.H>
#include <MacProjection.H>
#include <mfix_eb_if.H>

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

    void WriteCheckPointFile(std::string & check_file_name,
                             int nstep = 0, Real dt = 0.0, Real time = 0.0)  const;

    void WritePlotFile(std::string & plot_file_name,
                       int nstep = 0, Real dt = 0.0, Real time = 0.0) const;

    //! Save variables that don't change to plot file. The idea is that they can
    //! be saved _once_ per simulations and not repeatedly every Nth time step.
    void WriteStaticPlotFile(const std::string & plot_file_name) const;

    void Restart(std::string & restart_chkfile, int * nstep, Real * dt, Real * time,
                 IntVect & Nrep);

    void WriteParticleAscii(std::string & par_ascii_file_name, int nstep = 0) const;

    void WriteUSER(Real dt, Real time) const;

    void Regrid ();

    void Evolve(int nstep, int steady_state, Real & dt, Real & prev_dt, Real time, Real stop_time);

    void usr3();

    void output(int estatus, int finish, int nstep, Real dt, Real time)
    {
        pc->output( estatus, finish, nstep, dt, time);
    };

    void make_mg_bc(int mg_bc[]);

    void make_eb_geometry ();
    void make_eb_factories ();
    void fill_eb_levelsets ();
    void intersect_ls_walls ();


    //! Construct EB levels from Geometry shop. This builds each EB level and
    //! saves a pointer to each level into `eb_levels`. Also fills the implicit
    //! function values for each level in `implicit_functions`.
    template<class F> void build_eb_levels (EB2::GeometryShop<F> gshop)
    {

        /************************************************************************
         *                                                                      *
         * Build EB levels                                                      *
         *                                                                      *
         ***********************************************************************/

        for (int lev = 0; lev < nlev; lev ++)
            eb_levels[lev] = LSCore<F>::BuildEBLevel(gshop, lev, geom);


        /************************************************************************
         *                                                                      *
         * Fill each level with implicit function values                        *
         *                                                                      *
         ***********************************************************************/

        if (nlev == 1)
        {
            //___________________________________________________________________
            // Single-level case: fill normal and refined implicit functions

            // Get BoxArray and distribution mapping from EB level as the
            // mfix::grids and mfix::dmap vectors have not been built yet
            const BoxArray & ba            = eb_levels[0]->boxArray();
            const DistributionMapping & dm = eb_levels[0]->DistributionMap();

            // Base (fluid) level
            GShopLSFactory<F> gshop_lsf_1(gshop, geom[0], ba, dm, levelset__pad);
            implicit_functions[0] = gshop_lsf_1.fill_impfunc();

            // Refined level
            GShopLSFactory<F> gshop_lsf_2(gshop, geom[0], ba, dm,
                                          levelset__pad, levelset__refinement);
            implicit_functions[1] = gshop_lsf_2.fill_impfunc();

        }
        else
        {
            //___________________________________________________________________
            // Multi-level case: fill function for each level

            for (int lev = 0; lev < nlev; lev ++)
            {
                const BoxArray & ba            = eb_levels[lev]->boxArray();
                const DistributionMapping & dm = eb_levels[lev]->DistributionMap();

                GShopLSFactory<F> gshop_lsf(gshop, geom[lev], ba, dm, levelset__pad);
                implicit_functions[lev] = gshop_lsf.fill_impfunc();
            }

        }
    };

    std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_walls(bool & has_walls);

    std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_real_walls(bool & has_real_walls);

    std::unique_ptr<
        EB2::TranslationIF<EB2::PolynomialIF>
        > get_poly(int max_order, std::string field_prefix);

    std::unique_ptr<EB2::IntersectionIF<EB2::PlaneIF,EB2::PlaneIF,EB2::PlaneIF>>
        make_wall(int dir, // direction (long edge) of wall
                  Real position, Real height, Real width );

    void WriteEBSurface();

    void InitIOData();

    static std::string get_load_balance_type();

    // flag enabling level-set restart (i.e. prevent make_eb_* from rebuilding
    // the level-set data).
    bool levelset__restart = false;

    void mfix_compute_vort();

    ///
    /// Parameters describing the level-set grid:
    ///
    int levelset__refinement = 1;    // refinement (wrt particle grid) of the level-set's grid
    int levelset__eb_refinement = 1; // refinement of the EB facets used to generate the level-set
    int levelset__pad = 2;           // padding of the level-set grid
    int levelset__eb_pad = 2;        // padding of the EBIS used to generate the level-set (useful for
                                     // ensuring that neighbor grids are considered when filling LS)

     const int * get_bc_ilo_ptr(int lev) {return bc_ilo[lev]->dataPtr();}
     const int * get_bc_ihi_ptr(int lev) {return bc_ihi[lev]->dataPtr();}
     const int * get_bc_jlo_ptr(int lev) {return bc_jlo[lev]->dataPtr();}
     const int * get_bc_jhi_ptr(int lev) {return bc_jhi[lev]->dataPtr();}
     const int * get_bc_klo_ptr(int lev) {return bc_klo[lev]->dataPtr();}
     const int * get_bc_khi_ptr(int lev) {return bc_khi[lev]->dataPtr();}

     int get_nghost() { return nghost; }

protected:

     //! Tagging cells for refinement
     virtual void ErrorEst(int lev, TagBoxArray & tags, Real time, int ngrow) override;

     //! Make a new level using provided BoxArray and DistributionMapping and
     //! fill with interpolated coarse level data. Called by AmrCore::regrid.
     virtual void MakeNewLevelFromCoarse(int lev, Real time, const BoxArray & ba,
                                         const DistributionMapping & dm) override
     {
         amrex::Abort("MakeNewLevelFromCoarse: To be implemented");
     }

     //! Remake an existing level using provided BoxArray and
     //! DistributionMapping and fill with existing fine and coarse data. Called
     //! by AmrCore::regrid.
     virtual void RemakeLevel(int lev, Real time, const BoxArray & ba,
                              const DistributionMapping & dm) override
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

     //! The particles see an inflow face as a solid wall; thus instead of
     //! creating a separate EBfactory for particles, we simply modify the
     //! level-set near inflow to see it as a wall
     void mfix_set_ls_near_inflow();

     void mfix_project_velocity   ();
     void mfix_initial_iterations ( Real dt, Real stop_time, int steady_state );
     void mfix_apply_projection   ( Real dt, Real scaling_factor, bool proj_2 );

     void mfix_set_scalar_bcs ();
     void mfix_set_velocity_bcs (Real time, int extrap_dir_bcs);

     void FillPatchVel (int lev, Real time, MultiFab & mf, int icomp, int ncomp,
                        const Vector<BCRec> & bcr);
     void GetDataVel   (int lev, Real time, Vector<MultiFab *> & data,
                        Vector<Real> & datatime);

     void mfix_calc_volume_fraction(Real & sum_vol);
     void mfix_calc_drag_fluid(Real time);
     void mfix_calc_drag_particle(Real time);

     void AllocateArrays (int lev);

     void RegridArrays (int lev);

     void RegridLevelSetArray (int a_lev);

     void avgDown (int crse_lev, const MultiFab & S_fine, MultiFab & S_crse);

private:

     void make_eb_general();
     void make_eb_box();
     void make_eb_cylinder();
     void make_eb_hopper();
     void make_eb_cyclone();
     void make_eb_hourglass();
     void make_eb_clr();
     void make_eb_clr_riser();
     void make_eb_regular();

     void MakeBCArrays();

     void check_data();

     void MakeNewLevelFromScratch(int lev, Real time, const BoxArray & new_grids,
                                  const DistributionMapping & new_dmap) override;

     void ReMakeNewLevelFromScratch(int lev, const BoxArray & new_grids,
                                    const DistributionMapping & new_dmap);

     void WriteHeader(const std::string & name, int nstep, Real dt, Real time,
                      bool is_checkpoint) const;
     void WritePlotHeader(const std::string & name, int nstep, Real dt, Real time) const;
     void WriteCheckHeader(const std::string & name, int nstep, Real dt, Real time) const;

     void WriteJobInfo(const std::string & dir ) const;

     static void GotoNextLine(std::istream & is);

     void EvolveFluid(int nstep, int steady_state, Real & dt, Real & time, Real stop_time);

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
                                   Real dt);

     void mfix_compute_bcoeff_ppe();
     void mfix_compute_bcoeff_diff();

     void check_for_nans(int lev);

     void mfix_compute_dt(Real time, Real stop_time, int steady_state, Real & dt);

     void mfix_compute_ugradu_predictor(Vector< std::unique_ptr<MultiFab> >& conv,
                                        Vector< std::unique_ptr<MultiFab> >& vel,
                                        Real time);

     void mfix_compute_ugradu_corrector(Vector< std::unique_ptr<MultiFab> >& conv,
                                        Vector< std::unique_ptr<MultiFab> >& vel,
                                        Real time);

     void mfix_compute_divtau(int lev, MultiFab& divtau,
                              Vector< std::unique_ptr<MultiFab> >& vel);

     int  steady_state_reached(Real dt);

     void mfix_apply_predictor(Vector< std::unique_ptr<MultiFab> >& conv_old,
                               Vector< std::unique_ptr<MultiFab> >& divtau_old,
                               Real time, Real dt, bool proj_2);
     void mfix_apply_corrector(Vector< std::unique_ptr<MultiFab> >& conv_old,
                               Vector< std::unique_ptr<MultiFab> >& divtau_old,
                               Real time, Real dt, bool proj_2);

     void mfix_diffuse_velocity(Real time, Real dt);

     void mfix_add_drag_terms(Real dt);

     void mfix_add_gravity(int lev, Real dt, Vector< std::unique_ptr<MultiFab> >& vel);
#if 0
     void mfix_compute_velocity_slopes(int lev, Vector< std::unique_ptr<MultiFab> >& vel);
     void mfix_compute_velocity_at_faces(int lev, Vector< std::unique_ptr<MultiFab> >& vel);
#endif

     void mfix_compute_velocity_slopes(int lev, Real time, MultiFab& vel);
     void mfix_compute_MAC_velocity_at_faces(Real time, Vector< std::unique_ptr<MultiFab> >& vel);

     void mfix_compute_diveu(Real time);

     void mfix_print_max_vel(int lev);
     void mfix_print_max_gp (int lev);

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

     Real volWgtSum (int lev, const MultiFab & mf, int comp, bool local=false);

     static int nlev;

     Vector<int> istep = {1} ;

     // Particle container
     std::unique_ptr<MFIXParticleContainer> pc;

     // Unit vectors in Cartesian space
     static IntVect e_x;
     static IntVect e_y;
     static IntVect e_z;

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

     // Level-Set Data => used for particle-wall collisions and fluid
     // reconstruction in particle drag calculation. NOTE: this has at least 2
     // levels: 0) fluid reconstruction, 1..N) particle-wall (refined) collisions
     Vector< std::unique_ptr<MultiFab> > level_sets;
     Vector< std::unique_ptr<MultiFab> > implicit_functions;

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

     // Options to control mfix verbosity level
     int m_verbose = 0;

     // Options to control MLMG behavior
     int  mg_max_iter = 100;
     int  mg_cg_maxiter = 100;
     int  mg_max_fmg_iter = 0;
     Real mg_rtol = 1.0e-11;
     Real mg_atol = 1.0e-14;

     // Verbosity for the nodal solve
     int  mg_verbose = 0;
     int  mg_cg_verbose = 0;

     // Verbosity for the diffusion solve
     int  diff_mg_verbose = 0;
     int  diff_mg_cg_verbose = 0;

     // Tolerance to check for steady state
     Real steady_state_tol = 1.0e-5  ;

     int call_udf;
     bool dual_grid = false;

     // enable output_manager during des_time_loop
     bool subdt_io;

     // Parameters to control load balancing
     static std::string particle_init_type;
     static std::string load_balance_type;
     static std::string knapsack_weight_type;
     static int load_balance_fluid;
     static int knapsack_nmax;

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
     Vector<Real> avg_region_x_e;
     Vector<Real> avg_region_x_w;
     Vector<Real> avg_region_y_s;
     Vector<Real> avg_region_y_n;
     Vector<Real> avg_region_z_b;
     Vector<Real> avg_region_z_t;

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

     int explicit_diffusion = 1;

     bool do_initial_proj    = true;
     int  initial_iterations = 3;

     /***************************************************************************
      *                                                                         *
      * I/O stuff                                                               *
      *                                                                         *
      **************************************************************************/

     // Variables to simplify IO -- the vectors are the same for plt/chk
     Vector< const Vector< std::unique_ptr<MultiFab> >* > vectorVars;
     Vector< std::string > vecVarsName;

     Vector< const Vector< std::unique_ptr<MultiFab> >* > chkscalarVars;
     Vector< std::string > pltscaVarsName;

     Vector< const Vector< std::unique_ptr<MultiFab> >* > pltscalarVars;
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

     //! EB levels
     Vector<const EB2::Level *> eb_levels;

     //! EB factory that lives on the fluid grids
     Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > ebfactory;
     //! EB factory that lives on the particle grids
     Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> > particle_ebfactory;

     // Max level at which to solve the fluid equations
     int amr_max_level = 0;

     bool mfix_update_ebfactory (int a_lev);

};

inline std::string mfix::get_load_balance_type()
{
    return load_balance_type;
}

#endif
