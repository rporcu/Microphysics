#ifndef DIFFUSION_EQUATION_H_
#define DIFFUSION_EQUATION_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLEBTensorOp.H>

#include <mfix_fluid_parms.H>

//
// Solver for the implicit part of the diffusion equation:
//
// ( 1 - dt / rho * div ( eta grad ) ) u* = rhs
//
// Input "vel" holds the right hand side of this equation:
//
//  rhs = u + dt * ( conv + div ( eta (grad u)^T ) / rho + g - grad(p + p0) / rho )
//
// This is replaced by the solution of the matrix solve.
// Output "vel" therefore holds u*
//
// Note: we actually solve the above equation multiplied by the density ro.
//

class DiffusionOp
{
  public:
    // Constructor, with initialisation
    DiffusionOp (amrex::AmrCore* _amrcore,
                 amrex::Vector< const amrex::EBFArrayBoxFactory* >const& _ebfactory,
                 FluidPhase& _fluid,
                 std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_velbc_lo,
                 std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_velbc_hi,
                 std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_scalbc_lo,
                 std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_scalbc_hi,
                 std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_temperaturebc_lo,
                 std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_temperaturebc_hi,
                 std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_specieslbc_lo,
                 std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_specieslbc_hi,
                 int _nghost);

    // Destructor
    ~DiffusionOp ();

    // Actually do the setup work here
    void setup (amrex::AmrCore* _amrcore,
                amrex::Vector< const amrex::EBFArrayBoxFactory* >const& _ebfactory);

    // Read input from ParmParse (solver settings)
    void readParameters ();

    // Set user-supplied solver settings (must be done every time step)
    void setSolverSettings (amrex::MLMG& solver);

    // Tensor viscous solve
    void diffuse_velocity (const amrex::Vector< amrex::MultiFab* >& vel,
                           const amrex::Vector< amrex::MultiFab* >& ep_ro,
                           const amrex::Vector< amrex::MultiFab* >& T_g,
                           const int advect_enthalpy,
                           amrex::Real dt);

    // Temperature diffusion
    void diffuse_temperature (const amrex::Vector< amrex::MultiFab* >& T_g,
                              const amrex::Vector< amrex::MultiFab* >& ep_g,
                              const amrex::Vector< amrex::MultiFab* >& ro_g,
                              const amrex::Vector< amrex::MultiFab* >& h_g,
                              const amrex::Vector< amrex::MultiFab* >& T_g_on_eb,
                              amrex::Real dt);

    // Species mass fraction diffusion
    void diffuse_species (const amrex::Vector< amrex::MultiFab* >& X_gk,
                          const amrex::Vector< amrex::MultiFab* >& ep_ro_g,
                          const amrex::Vector< amrex::MultiFab* >& T_g,
                                amrex::Vector<amrex::BCRec> const& h_bcrec,
                          const amrex::Real dt);

    // Scalar diffusion
    void diffuse_scalar (const amrex::Vector< amrex::MultiFab* >& scal,
                         const amrex::Vector< amrex::MultiFab* >& ep_ro_g,
                         const amrex::Vector< amrex::Real > mu_s,
                               amrex::Vector<amrex::BCRec> const& h_bcrec,
                         amrex::Real dt);

    void diffuse_drag (const amrex::Vector< amrex::MultiFab* >& drag,
                       amrex::Real doceff);

    void diffuse_volfrac (const amrex::Vector< amrex::MultiFab* >& scal,
                          amrex::Real doceff);

    // Explicit tensor viscous term
    void ComputeDivTau (const amrex::Vector< amrex::MultiFab* >& divtau_out,
                        const amrex::Vector< amrex::MultiFab* >& vel_in,
                        const amrex::Vector< amrex::MultiFab* >& ep_in,
                        const amrex::Vector< amrex::MultiFab* >& T_g_in,
                        const int advect_enthalpy);

    // Explicit tensor temperature term
    void ComputeLapT (const amrex::Vector< amrex::MultiFab*      >& lapT_out,
                      const amrex::Vector< amrex::MultiFab*      >& T_g_in,
                      const amrex::Vector< amrex::MultiFab const*>& ep_g_in,
                      const amrex::Vector< amrex::MultiFab const*>& T_g_on_eb);

    // Explicit tensor species mass fractions term
    void ComputeLapX (const amrex::Vector< amrex::MultiFab*      >& lapX_out,
                      const amrex::Vector< amrex::MultiFab*      >& X_gk_in,
                      const amrex::Vector< amrex::MultiFab const*>& ro_g_in,
                      const amrex::Vector< amrex::MultiFab const*>& ep_g_in,
                      const amrex::Vector< amrex::MultiFab const*>& T_g_in);

    // Explicit tensor species mass fractions term
    void SubtractDivXGX (const amrex::Vector< amrex::MultiFab*      >& X_gk_in,
                         const amrex::Vector< amrex::MultiFab const*>& ro_g_in,
                         const amrex::Vector< amrex::MultiFab const*>& ep_g_in,
                         const amrex::Vector< amrex::MultiFab const*>& T_g_in,
                         const amrex::Real& dt);

    // Enthalpy tensor due to species mass fractions fluxes
    void ComputeLaphX (const amrex::Vector< amrex::MultiFab*       >& lap_hX_out,
                       const amrex::Vector< amrex::MultiFab*       >& X_gk_in,
                       const amrex::Vector< amrex::MultiFab const* >& ro_in,
                       const amrex::Vector< amrex::MultiFab const* >& ep_in,
                       const amrex::Vector< amrex::MultiFab const* >& T_g_in);

    // Explicit tensor viscous term // tracer
    void ComputeLapS (const amrex::Vector< amrex::MultiFab*      >& laps_out,
                      const amrex::Vector< amrex::MultiFab*      >& scal_in,
                      const amrex::Vector< amrex::MultiFab const*>& ro_in,
                      const amrex::Vector< amrex::MultiFab const*>& ep_in,
                      const amrex::Vector< amrex::Real >& mu_s);

  private:
    // AmrCore data
    amrex::AmrCore* amrcore;
    amrex::Vector< const amrex::EBFArrayBoxFactory* > ebfactory;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;

    int nghost = 0;

    // Boundary conditions
    std::array<amrex::LinOpBCType,3> m_velbc_lo;
    std::array<amrex::LinOpBCType,3> m_velbc_hi;
    std::array<amrex::LinOpBCType,3> m_scalbc_lo;
    std::array<amrex::LinOpBCType,3> m_scalbc_hi;
    std::array<amrex::LinOpBCType,3> m_temperaturebc_lo;
    std::array<amrex::LinOpBCType,3> m_temperaturebc_hi;
    std::array<amrex::LinOpBCType,3> m_speciesbc_lo;
    std::array<amrex::LinOpBCType,3> m_speciesbc_hi;

    // Internal data used in the matrix solve
    //
    // MultiLevel EmbeddedBoundary cell-centered Laplacian:
    //
    // ( alpha a - beta div ( b grad ) ) phi = rhs
    //
    std::unique_ptr<amrex::MLEBTensorOp> vel_matrix;

    amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM>> b;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> phi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> rhs;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vel_eb;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> species_phi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> species_rhs;
    //
    // MultiLevel EmbeddedBoundary cell-centered Laplacian:
    //
    // ( alpha a - beta div ( b grad ) ) phi = rhs
    //
    std::unique_ptr<amrex::MLEBABecLap> scal_matrix;
    std::unique_ptr<amrex::MLEBABecLap> temperature_matrix;
    std::unique_ptr<amrex::MLEBABecLap> species_matrix;

    amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM>> species_b;

    // DiffusionOp verbosity
    int verbose = 0;

    // Dirichlet vs Neumann bc's at EB walls
    bool eb_is_dirichlet;
    bool eb_temperature_is_dirichlet;

    // Options to control MLMG behavior
    int mg_verbose = 0;
    int mg_bottom_verbose = 0;
    int mg_maxiter = 100;
    int mg_bottom_maxiter = 100;
    int mg_max_fmg_iter = 0;
    int mg_max_coarsening_level = 100;
    int mg_agg_grid_size = -1;
    amrex::Real mg_rtol = 1.0e-11;
    amrex::Real mg_atol = 1.0e-14;
    std::string bottom_solver_type = "bicgstab";

    FluidPhase& fluid;
};

#endif
