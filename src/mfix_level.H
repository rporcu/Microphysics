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


class mfix_level
    : public AmrCore
{
public:

    mfix_level ();
    ~mfix_level ();

    void InitParams(int solve_fluid, int solve_dem, int max_nit, int call_udf );

    void Init (int lev, Real dt, Real time);

    void InitLevelData (int lev, Real dt, Real time);

    void PostInit (int lev, Real dt, Real time, int nstep, int restart_flag);

    void WriteCheckPointFile ( std::string& check_file_name, int nstep = 0, Real dt = 0.0, Real time = 0.0 )  const;

    void WritePlotFile ( std::string& plot_file_name, int nstep = 0, Real dt = 0.0, Real time = 0.0 ) const;

    void Restart (std::string& restart_chkfile, int *nstep, Real *dt, Real *time,
                  IntVect& Nrep);

    void WriteParticleAscii (std::string& par_ascii_file_name, int nstep = 0) const;

    void Regrid (int lev, int nstep, int dual_grid);

    void RegridOnRestart (int lev);

    void Evolve(int lev, int nstep, int set_normg,
                       Real& dt, Real& prev_dt, Real time, Real normg);

    void usr3 (int lev);

    void output (int lev, int estatus, int finish, int nstep, Real dt, Real time)
       { pc -> output( lev, estatus, finish, nstep, dt, time);};

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

    void mfix_calc_trd_and_tau (int lev);

    void mfix_init_fluid (int lev, int is_restarting = 0);
    void mfix_set_bc0 (int lev);

    void mfix_conv_rop (int lev);

    void mfix_calc_mflux (int lev);

    void mfix_solve_for_u    (int lev, Real dt, Real& num_u, Real& denom_u);
    void mfix_solve_for_v    (int lev, Real dt, Real& num_v, Real& denom_v);
    void mfix_solve_for_w    (int lev, Real dt, Real& num_w, Real& denom_w);

    // Projection
    void mfix_apply_projection (int lev, amrex::Real dt);
    // End projection
    
    void mfix_solve_linear_equation(int eq_id,int lev,MultiFab& sol, MultiFab& matrix, MultiFab& rhs);

    void mfix_solve_for_pp (int lev, Real dt, Real& num_p, Real& denom_p);

    void mfix_physical_prop (int lev, int calc_flag);

    void mfix_correct_0 (int lev);

    int solve_bicgstab (MultiFab& sol, const MultiFab& rhs, const MultiFab& matrix,
           int sweep_type, int precond_type, int maxiter, Real eps_rel, int lev);

    void mfix_calc_volume_fraction(int lev, Real& sum_vol);
    void mfix_calc_drag_fluid(int lev);
    void mfix_calc_drag_particle(int lev);

    void AllocateArrays(int lev);
    void RegridArrays(int lev, BoxArray& new_grids, DistributionMapping& new_dmap);

private:

    void check_data (int lev);

    void InitIOData ();

    void MakeNewLevelFromScratch(int lev, Real time,
          const BoxArray& new_grids, const DistributionMapping& new_dmap);

    void ReMakeNewLevelFromScratch(int lev, 
          const BoxArray& new_grids, const DistributionMapping& new_dmap);

    void WriteHeader(const std::string& name, int nstep, Real dt, Real time) const;

    void WriteJobInfo(const std::string& dir ) const;

    static void GotoNextLine (std::istream& is);

    void EvolveFluid (int lev, int nstep, int set_normg,
		      Real dt, Real& prev_dt, Real time, Real normg);

    void EvolveFluidProjection  (int lev, int nstep, int set_normg,
				 Real& dt, Real& prev_dt, Real time, Real normg);

    // Projection-related methods
    //void solve_poisson_equation ( MultiFab& b, MultiFab& phi, MultiFab& rhs );
    void solve_poisson_equation ( int lev,
				  Vector< Vector< std::unique_ptr<MultiFab> > >& b,
				  Vector< std::unique_ptr<MultiFab> >& phi,
				  Vector< std::unique_ptr<MultiFab> >& rhs );

    void compute_oro_g ( int lev );
    
    void init_tests_projection (int lev);

    void check_for_nans (int lev);
    
    amrex::Real compute_dt (int lev);

    void mfix_apply_pcm_prediction (int lev, amrex::Real dt);

    void mfix_apply_pcm_correction (int lev, amrex::Real dt);

    void mfix_apply_forcing_terms (int lev, amrex::Real dt);
    
// End of projection-related methods


    
    // Make a level 0 grids covering the whole domain.  It does NOT install the new grids.
    BoxArray MakeBaseGrids () const;
    void ChopGrids (const Box& domain, BoxArray& ba, int target_size) const;

    void mfix_set_bc_type(int lev);
    void mfix_set_bc1(int lev);

    void fill_mf_bc(int lev,  MultiFab& mf);
#if 0
    void WritePlotFile () const;

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

    // Pressure correction equation
    Vector< std::unique_ptr<MultiFab> > pp_g;
    Vector< std::unique_ptr<MultiFab> > d_e ;
    Vector< std::unique_ptr<MultiFab> > d_n ;
    Vector< std::unique_ptr<MultiFab> > d_t ;

    // Molecular viscosity
    Vector< std::unique_ptr<MultiFab> > mu_g ;

    // Cell-based
    Vector< std::unique_ptr<MultiFab> > lambda_g;
    Vector< std::unique_ptr<MultiFab> > trD_g;

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


    // Projection
    Vector< std::unique_ptr<MultiFab> > ugradu_x;
    Vector< std::unique_ptr<MultiFab> > ugradu_y;
    Vector< std::unique_ptr<MultiFab> > ugradu_z;
    
    Vector< std::unique_ptr<MultiFab> > divtau_x;
    Vector< std::unique_ptr<MultiFab> > divtau_y;
    Vector< std::unique_ptr<MultiFab> > divtau_z;

    Vector< std::unique_ptr<MultiFab> > fp_x;
    Vector< std::unique_ptr<MultiFab> > fp_y;
    Vector< std::unique_ptr<MultiFab> > fp_z;

    
    Vector< Vector< std::unique_ptr<MultiFab> > > oro_g; // Over ro_g = ep_g / rop_g

    Vector< std::unique_ptr<MultiFab> > phi; // Pressure increment
    
    
    int solve_fluid;
    int solve_dem;
    int use_proj_method = 0; // Whether to use projection to solve NS eqs.
    int max_nit;
    int call_udf;

    static std::string particle_init_type;
    static std::string load_balance_type;

    Real cfl = 0.7;

    Real sum_vol_orig;

    int nghost_bc = 2;

    // Variables to simplify IO
    Vector< const Vector< std::unique_ptr<MultiFab> >* > vectorVars;
    Vector< const Vector< std::unique_ptr<MultiFab> >* > scalarVars;
    Vector< std::string > vecVarsName;
    Vector< std::string > scaVarsName;

};

#endif
