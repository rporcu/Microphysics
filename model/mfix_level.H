#ifndef MFIX_LEVEL_H_
#define MFIX_LEVEL_H_

#include <iostream>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AmrCore.H>
#include <iMultiFab.H>
#include <BLProfiler.H>

#include <ParticleContainer.H>

class mfix_level
    : public AmrCore
{
public:

    mfix_level (const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord);
    mfix_level ();
    ~mfix_level();

    int Verbose () const { return verbose; }

    void InitParams(int solve_fluid, int solve_dem, int cyclic_mf, int max_nit, int call_udf);

    void Init (int lev, Real dt, Real time);

    void Restart();

    void evolve_fluid (int lev, int nstep, int set_normg, 
                       Real dt, Real& prev_dt, Real time, Real normg); 

    void evolve_dem (int lev, int nstep, Real dt, Real time);

    void usr3 (int lev);

    void output (int lev, int estatus, int finish, int nstep, Real dt, Real time);

    void make_mg_bc (int mg_bc[]);

protected:
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/) override {}

    void mfix_calc_coeffs (int lev, int calc_flag);

    void mfix_calc_all_coeffs (int lev);

    void mfix_calc_trd_and_tau (int lev);

    void mfix_init_fluid (int lev);

    void mfix_comp_mean_fields (int lev);

    void mfix_conv_rop (int lev, Real dt);

    void mfix_calc_mflux (int lev);

    void mfix_set_wall_bc (int lev);

    void mfix_set_bc1 (int lev, Real time, Real dt);

    void mfix_solve_for_vels (int lev, Real dt);

    void mfix_solve_linear_system(int eq_id,int lev,MultiFab& sol, MultiFab& matrix, MultiFab& rhs);

    void mfix_solve_for_pp (int lev, Real dt, Real& lnormg, Real& resg);

    void mfix_physical_prop (int lev, int calc_flag);

    void mfix_correct0 (int lev);

private:

    void ReadParameters ();

    void InitLevelData (int lev, Real dt, Real time);

    void MakeNewLevel (int lev, Real time,
		       const BoxArray& new_grids, const DistributionMapping& new_dmap);

    // Make a level 0 grids covering the whole domain.  It does NOT install the new grids.
    BoxArray MakeBaseGrids () const;

#if 0
    void WritePlotFile () const;

    Array<int> istep;      // which step?
    Array<int> nsubsteps;  // how many substeps on each level?
#endif

    // Particle container
    std::unique_ptr<MyParticleContainer> mypc;

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

    Array< std::unique_ptr<MultiFab> > flux_gE ;
    Array< std::unique_ptr<MultiFab> > flux_gN ;
    Array< std::unique_ptr<MultiFab> > flux_gT ;

    Array< std::unique_ptr<MultiFab> > rop_gE ;
    Array< std::unique_ptr<MultiFab> > rop_gN ;
    Array< std::unique_ptr<MultiFab> > rop_gT ;

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

    int verbose = 1;

    int solve_fluid;
    int solve_dem;
    int cyclic_mf;
    int max_nit;
    int call_udf;

//  int max_step   = std::numeric_limits<int>::max();
//  Real stop_time = std::numeric_limits<Real>::max();

    Real cfl = 0.7;

    std::string restart_chkfile;

    std::string check_file {"chk"};
    std::string plot_file {"plt"};
    int check_int = -1;
    int plot_int = -1;
};

#endif
