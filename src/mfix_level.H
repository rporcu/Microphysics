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

#include <ParticleContainer.H>

class mfix_level
    : public AmrCore
{
public:

//  mfix_level (const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord);
    mfix_level ();
    ~mfix_level();

    int Verbose () const { return verbose; }

    void InitParams(int solve_fluid, int solve_dem, int cyclic_mf, int max_nit, int call_udf);

    void Init (int lev, Real dt, Real time);

    void Restart();

    void WriteCheckPointFile ( int nstep = 0, Real dt = 0.0, Real time = 0.0 )  const;

    int IsRestartEnabled () const;

    void WritePlotFile ( int nstep = 0, Real dt = 0.0, Real time = 0.0 ) const;

    void InitFromCheckpoint (int *nstep, Real *dt, Real *time) const;

    void evolve_fluid (int lev, int nstep, int set_normg,
                       Real dt, Real& prev_dt, Real time, Real normg);

    void evolve_dem (int lev, int nstep, Real dt, Real time);

    void usr3 (int lev);

    void output (int lev, int estatus, int finish, int nstep, Real dt, Real time);

    void make_mg_bc (int mg_bc[]);

protected:

    //! Tagging cells for refinement
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/) override
        { amrex::Abort("ErrorEst: To be implemented"); }

    //! Make a new level using provided BoxArray and
    //! DistributionMapping and fill with interpolated coarse level
    //! data.  Called by AmrCore::regrid.
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) override
        { amrex::Abort("MakeNewLevelFromCoarse: To be implemented"); }

    //! Remake an existing level using provided BoxArray and
    //! DistributionMapping and fill with existing fine and coarse
    //! data.  Called by AmrCore::regrid.
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm) override
        { amrex::Abort("RemakeLevel: To be implemented"); }

    //! Delete level data.  Called by AmrCore::regrid.
    virtual void ClearLevel (int lev) override
        { amrex::Abort("ClearLevel: To be implemented"); }

    void mfix_calc_coeffs (int lev, int calc_flag);

    void mfix_calc_all_coeffs (int lev);

    void mfix_calc_trd_and_tau (int lev);

    void mfix_init_fluid (int lev);

    void mfix_comp_mean_fields (int lev);

    void mfix_conv_rop (int lev, Real dt);

    void mfix_calc_mflux (int lev);

    void mfix_solve_for_vels (int lev, Real dt, Real (&residuals)[16]);

    void mfix_solve_linear_equation(int eq_id,int lev,MultiFab& sol, MultiFab& matrix, MultiFab& rhs);

    void mfix_solve_for_pp (int lev, Real dt, Real& lnormg, Real& resg, Real (&residuals)[16]);

    void mfix_physical_prop (int lev, int calc_flag);

    void mfix_correct_0 (int lev);

    int solve_bicgstab (MultiFab& sol, const MultiFab& rhs, const MultiFab& matrix,
           int sweep_type, int precond_type, int maxiter, Real eps_rel, int lev);

private:

    void ReadParameters ();

    void InitLevelData (int lev, Real dt, Real time);

    void InitIOData ();

    void MakeNewLevelFromScratch(int lev, Real time,
           const BoxArray& new_grids, const DistributionMapping& new_dmap);

    void WriteMfixHeader(const std::string& name, int nstep, Real dt, Real time) const;

    void WriteJobInfo(const std::string& dir ) const;

    static void GotoNextLine (std::istream& is);

    // Make a level 0 grids covering the whole domain.  It does NOT install the new grids.
    BoxArray MakeBaseGrids () const;

  void mfix_set_bc_type(int lev);

  void fill_mf_bc(int lev,  MultiFab& mf);
#if 0
    void WritePlotFile () const;

      // which step?
    Array<int> nsubsteps;  // how many substeps on each level?
#endif

    Array<int> istep = {1} ;

    // Particle container
    std::unique_ptr<MFIXParticleContainer> pc;

//  Array< std::unique_ptr<iMultiFab> > bc_type;
//  Array< std::unique_ptr<iMultiFab> > bc_ptr;

    IArrayBox bc_ilo;
    IArrayBox bc_ihi;
    IArrayBox bc_jlo;
    IArrayBox bc_jhi;
    IArrayBox bc_klo;
    IArrayBox bc_khi;

    Array< std::unique_ptr<iMultiFab> > flag;

    // Matrix and rhs vector
    Array< std::unique_ptr<MultiFab> > A_m;
    Array< std::unique_ptr<MultiFab> > b_m;

    // Void fraction
    Array< std::unique_ptr<MultiFab> > ep_g ;
    Array< std::unique_ptr<MultiFab> > ep_go;

    // Gas pressure fraction
    Array< std::unique_ptr<MultiFab> > p_g ;
    Array< std::unique_ptr<MultiFab> > p_go;

    // Gas density
    Array< std::unique_ptr<MultiFab> > ro_g ;
    Array< std::unique_ptr<MultiFab> > ro_go;

    // Gas bulk density
    Array< std::unique_ptr<MultiFab> > rop_g ;
    Array< std::unique_ptr<MultiFab> > rop_go;

    // X-axis gas velocity
    Array< std::unique_ptr<MultiFab> > u_g ;
    Array< std::unique_ptr<MultiFab> > u_go;
    Array< std::unique_ptr<MultiFab> > u_gt;

    // Y-axis gas velocity
    Array< std::unique_ptr<MultiFab> > v_g ;
    Array< std::unique_ptr<MultiFab> > v_go;
    Array< std::unique_ptr<MultiFab> > v_gt;

    // Z-axis gas velocity
    Array< std::unique_ptr<MultiFab> > w_g ;
    Array< std::unique_ptr<MultiFab> > w_go;
    Array< std::unique_ptr<MultiFab> > w_gt;

    // Pressure correction equation
    Array< std::unique_ptr<MultiFab> > pp_g;
    Array< std::unique_ptr<MultiFab> > d_e ;
    Array< std::unique_ptr<MultiFab> > d_n ;
    Array< std::unique_ptr<MultiFab> > d_t ;

    // Molecular viscosity
    Array< std::unique_ptr<MultiFab> > mu_g ;

    //
    Array< std::unique_ptr<MultiFab> > lambda_g ;
    Array< std::unique_ptr<MultiFab> > trD_g ;

    //
    Array< std::unique_ptr<MultiFab> > tau_u_g ;
    Array< std::unique_ptr<MultiFab> > tau_v_g ;
    Array< std::unique_ptr<MultiFab> > tau_w_g ;

    Array< std::unique_ptr<MultiFab> > fluxX ;
    Array< std::unique_ptr<MultiFab> > fluxY ;
    Array< std::unique_ptr<MultiFab> > fluxZ ;

    Array< std::unique_ptr<MultiFab> > ropX ;
    Array< std::unique_ptr<MultiFab> > ropY ;
    Array< std::unique_ptr<MultiFab> > ropZ ;

    //
    Array< std::unique_ptr<MultiFab> > f_gds  ;
    Array< std::unique_ptr<MultiFab> > drag_bm;


    Array<int> particle_state;
    Array<int> particle_phase;

    Array<Real> des_radius;
    Array<Real> ro_sol;
    Array<Real> pvol;
    Array<Real> pmass;
    Array<Real> omoi;
    Array<Real> des_pos_new;
    Array<Real> des_vel_new;
    Array<Real> des_usr_var;
    Array<Real> omega_new;
    Array<Real> des_acc_old;
    Array<Real> rot_acc_old;
    Array<Real> drag_fc;
    Array<Real> fc;
    Array<Real> tow;
    Array<int> pairs;

    // runtime parameters

    int verbose = 0;

    int solve_fluid;
    int solve_dem;
    int cyclic_mf;
    int max_nit;
    int call_udf;

    int max_step;
    Real stop_time;

    Real cfl = 0.7;

    std::string restart_chkfile {"no_restart"};

    std::string check_file {"chk"};
    std::string plot_file {"plt"};
    int check_int = -1;
    int plot_int = -1;


    // Variables to simplify IO
    Array< const Array< std::unique_ptr<MultiFab> >* > vectorVars;
    Array< const Array< std::unique_ptr<MultiFab> >* > scalarVars;
    Array< std::string > vecVarsName;
    Array< std::string > scaVarsName;

};

#endif
