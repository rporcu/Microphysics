#ifndef MFIX_NODAL_PROJECTION_H
#define MFIX_NODAL_PROJECTION_H

#include <AMReX.H>
#include <mfix.H>

class NodalProjection
{

public:

    NodalProjection () {};
    NodalProjection (const mfix* a_mfix,
                     std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                     std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_hi )
        : m_mfix(a_mfix), m_bc_lo(a_bc_lo), m_bc_hi(a_bc_hi), m_ok(true) {};

    void define ( const mfix* a_mfix,
                  std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                  std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_hi );

    void project (      Vector< std::unique_ptr< amrex::MultiFab > >& a_vel,
                  const Vector< std::unique_ptr< amrex::MultiFab > >& a_ep,
                  const Vector< std::unique_ptr< amrex::MultiFab > >& a_ro,
                  const Vector< std::unique_ptr< amrex::MultiFab > >& a_depdt,
                        Real a_time, Real a_scale_factor = 1.0 );

    Vector< const amrex::MultiFab* > getPhi () const {return GetVecOfConstPtrs(m_phi);};

    Vector< const amrex::MultiFab* > getGradPhi () const {return GetVecOfConstPtrs(m_fluxes);};

    Vector< const amrex::MultiFab* > getRHS () const {return GetVecOfConstPtrs(m_rhs);};


private:

    bool m_ok = false;

    // Verbosity
    int  m_mg_verbose = 0;
    int  m_mg_cg_verbose = 0;

    // Control MLMG behavior
    int  m_mg_maxiter = 100;
    int  m_mg_cg_maxiter = 100;
    Real m_mg_rtol = 1.0e-11;
    Real m_mg_atol = 1.0e-14;
    std::string m_bottom_solver_type = "bicgcg";

    // Max coarsening level
    int  m_mg_max_coarsening_level = 100;

    // Boundary conditions
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_bc_lo;
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_bc_hi;

    // We need ebfactory, amrcore and some functions defined in mfix.
    // Thus we make this class a friend of mfix' and use its methods.
    // In the future we will decouple those method from mfix and just pass
    // amrcore and ebfactory to the constructor
    const mfix* m_mfix;

    // Cell-centered data
    Vector< std::unique_ptr< amrex::MultiFab > >  m_fluxes;
    Vector< std::unique_ptr< amrex::MultiFab > >  m_sigma;

    // Node-centered data
    Vector< std::unique_ptr< amrex::MultiFab > >  m_phi;
    Vector< std::unique_ptr< amrex::MultiFab > >  m_rhs;

    // Linear operator
    std::unique_ptr< amrex::MLNodeLaplacian > m_matrix;

    // Solver
    std::unique_ptr< amrex::MLMG > m_solver;


    void readParameters ();

    void setup ();

    void computeRHS (       Vector< std::unique_ptr< amrex::MultiFab > >& a_vel,
                      const Vector< std::unique_ptr< amrex::MultiFab > >& a_ep,
                      const Vector< std::unique_ptr< amrex::MultiFab > >& a_depdt,
                            Real a_time );

};


#endif
