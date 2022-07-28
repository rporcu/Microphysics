#ifndef _MFIX_BC_H_
#define _MFIX_BC_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>

#include <AMReX_Vector.H>
#include <AMReX_IArrayBox.H>

#include <AMReX_Vector.H>
#include <AMReX_LO_BCTYPES.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>

#include <mfix_eb.H>
#include <mfix_regions.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_pic.H>

#include <bits/stdc++.h>


// Boundary conditions list auxiliary structure
//
struct BCList {

  public:
    // Types
    enum {undefined=0, ig=9, pinf=10, pout=11, minf=20, nsw=100, eb=111};

    // Constructor
    BCList (const int nlev_in);

    // Destructor
    ~BCList();

    void MakeBCArrays (int nghost,
                       bool ooo_debug,
                       amrex::Vector<amrex::Geometry>& geom);

  private:
    const int nlev;

  public:
    // Boundary conditions types
    amrex::Vector<amrex::IArrayBox*> bc_ilo;
    amrex::Vector<amrex::IArrayBox*> bc_ihi;
    amrex::Vector<amrex::IArrayBox*> bc_jlo;
    amrex::Vector<amrex::IArrayBox*> bc_jhi;
    amrex::Vector<amrex::IArrayBox*> bc_klo;
    amrex::Vector<amrex::IArrayBox*> bc_khi;
};


// Forward declaration
struct BC_t;


class MFIXBoundaryConditions
{
  public:
    // Constructor
    MFIXBoundaryConditions (MFIXEmbeddedBoundaries& embedded_boundaries);

    void Initialize (amrex::Geometry& geom,
                     const MFIXRegions& regions,
                     MFIXFluidPhase& fluid,
                     const MFIXSolidsPhase& solids,
                     MFIXDEM& dem,
                     MFIXPIC& pic);

    void read_bc_density (amrex::ParmParse pp,
                          FLUID_t* fluid);

    void read_bc_velocity (amrex::ParmParse pp,
                           FLUID_t* fluid);

    void read_bc_volflow (amrex::ParmParse pp,
                          FLUID_t* fluid);

    void read_bc_temperature (amrex::ParmParse pp,
                              FLUID_t* fluid);

    void read_bc_species (amrex::ParmParse pp,
                          const amrex::Vector<std::string> species,
                          FLUID_t* fluid);

    const amrex::Vector<int>& bc_xlo () const { return m_bc_lo[0]; }
    const amrex::Vector<int>& bc_xhi () const { return m_bc_hi[0]; }
    const amrex::Vector<int>& bc_ylo () const { return m_bc_lo[1]; }
    const amrex::Vector<int>& bc_yhi () const { return m_bc_hi[1]; }
    const amrex::Vector<int>& bc_zlo () const { return m_bc_lo[2]; }
    const amrex::Vector<int>& bc_zhi () const { return m_bc_hi[2]; }

    int bc_xlo (const int i) const { return m_bc_lo[0][i]; }
    int bc_xhi (const int i) const { return m_bc_hi[0][i]; }
    int bc_ylo (const int i) const { return m_bc_lo[1][i]; }
    int bc_yhi (const int i) const { return m_bc_hi[1][i]; }
    int bc_zlo (const int i) const { return m_bc_lo[2][i]; }
    int bc_zhi (const int i) const { return m_bc_hi[2][i]; }

    // Flag to solve fluid equations
    int delp_dir () const { return m_delp_dir; }
    const amrex::RealVect& delp () const { return m_delp; }
    amrex::Real delp (const int i) const { return m_delp[i]; }
    
    void set_delp (const int i, const amrex::Real val) { m_delp[i] = val; }

    // Flags for DEM domain extent reinforcement.
    const std::array<int,6>& domain_bc () const { return m_domain_bc; }
    int domain_bc (const int i) const { return m_domain_bc[i]; }

    const std::bitset<6>& flow_plane () const { return m_flow_plane; }

    // Vectors storing EB planes for level-set creation
    const amrex::Vector<amrex::EB2::PlaneIF>& flow_planes () const { return m_flow_planes; };
    const amrex::Vector<amrex::EB2::PlaneIF>& wall_planes () const { return m_wall_planes; }

    // Use AMReX LinOpBCType to define the domain extents
    const std::array<amrex::LinOpBCType,3>& ppe_lobc () const { return m_ppe_lobc; }
    const std::array<amrex::LinOpBCType,3>& ppe_hibc () const { return m_ppe_hibc; }

    const std::array<amrex::LinOpBCType,3>& diff_vel_lobc () const { return m_diff_vel_lobc; }
    const std::array<amrex::LinOpBCType,3>& diff_vel_hibc () const { return m_diff_vel_hibc; }

    const std::array<amrex::LinOpBCType,3>& diff_scal_lobc () const { return m_diff_scal_lobc; }
    const std::array<amrex::LinOpBCType,3>& diff_scal_hibc () const { return m_diff_scal_hibc; }

    const std::array<amrex::LinOpBCType,3>& diff_temperature_lobc () const { return m_diff_temperature_lobc; }
    const std::array<amrex::LinOpBCType,3>& diff_temperature_hibc () const { return m_diff_temperature_hibc; }

    const std::array<amrex::LinOpBCType,3>& diff_species_lobc () const { return m_diff_species_lobc; }
    const std::array<amrex::LinOpBCType,3>& diff_species_hibc () const { return m_diff_species_hibc; }

    // Data structure for BC data
    const amrex::Vector<BC_t>& bc () const { return m_bc; }
    BC_t& bc (const int bcv) { return m_bc[bcv]; }
    const BC_t& bc (const int bcv) const { return m_bc[bcv]; }

    amrex::Real& h_bc_u_g (const int bcv) { return m_h_bc_u_g[bcv]; };
    amrex::Gpu::HostVector<amrex::Real>& h_bc_u_g () { return m_h_bc_u_g; };
    amrex::Gpu::DeviceVector<amrex::Real>& bc_u_g () { return m_bc_u_g; };

    amrex::Real& h_bc_v_g (const int bcv) { return m_h_bc_v_g[bcv]; };
    amrex::Gpu::HostVector<amrex::Real>& h_bc_v_g () { return m_h_bc_v_g; };
    amrex::Gpu::DeviceVector<amrex::Real>& bc_v_g () { return m_bc_v_g; };

    amrex::Real& h_bc_w_g (const int bcv) { return m_h_bc_w_g[bcv]; };
    amrex::Gpu::HostVector<amrex::Real>& h_bc_w_g () { return m_h_bc_w_g; };
    amrex::Gpu::DeviceVector<amrex::Real>& bc_w_g () { return m_bc_w_g; };

    amrex::Real& h_bc_t_g (const int bcv) { return m_h_bc_t_g[bcv]; };
    amrex::Gpu::HostVector<amrex::Real>& h_bc_t_g () { return m_h_bc_t_g; };
    amrex::Gpu::DeviceVector<amrex::Real>& bc_t_g () { return m_bc_t_g; };

    amrex::Real& h_bc_h_g (const int bcv) { return m_h_bc_h_g[bcv]; };
    amrex::Gpu::HostVector<amrex::Real>& h_bc_h_g () { return m_h_bc_h_g; };
    amrex::Gpu::DeviceVector<amrex::Real>& bc_h_g () { return m_bc_h_g; };

    amrex::Real& h_bc_ro_g (const int bcv) { return m_h_bc_ro_g[bcv]; };
    amrex::Gpu::HostVector<amrex::Real>& h_bc_ro_g () { return m_h_bc_ro_g; };
    amrex::Gpu::DeviceVector<amrex::Real>& bc_ro_g () { return m_bc_ro_g; };

    amrex::Real& h_bc_tracer (const int bcv) { return m_h_bc_tracer[bcv]; };
    amrex::Gpu::HostVector<amrex::Real>& h_bc_tracer () { return m_h_bc_tracer; };
    amrex::Gpu::DeviceVector<amrex::Real>& bc_tracer () { return m_bc_tracer; };

    amrex::Real& h_bc_ep_g (const int bcv) { return m_h_bc_ep_g[bcv]; };
    amrex::Gpu::HostVector<amrex::Real>& h_bc_ep_g () { return m_h_bc_ep_g; };
    amrex::Gpu::DeviceVector<amrex::Real>& bc_ep_g () { return m_bc_ep_g; };

    amrex::Real& h_bc_p_g (const int bcv) { return m_h_bc_p_g[bcv]; };
    amrex::Gpu::HostVector<amrex::Real>& h_bc_p_g () { return m_h_bc_p_g; };
    amrex::Gpu::DeviceVector<amrex::Real>& bc_p_g () { return m_bc_p_g; };

    amrex::Real& h_bc_X_gk (const int n, const int bcv) { return m_h_bc_X_gk[n][bcv]; };
    amrex::Gpu::HostVector<amrex::Real>& h_bc_X_gk (const int n) { return m_h_bc_X_gk[n]; };
    amrex::Vector<amrex::Gpu::HostVector<amrex::Real>>& h_bc_X_gk () { return m_h_bc_X_gk; };
    amrex::Gpu::DeviceVector<amrex::Real>& bc_X_gk (const int n) { return m_bc_X_gk[n]; };
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>>& bc_X_gk () { return m_bc_X_gk; };
    amrex::Real* h_bc_X_gk_ptr (const int n) { return m_h_bc_X_gk_ptr[n]; };

    void set_h_bc_X_gk_ptr (const int n)
    {
      if (m_bc_X_gk[n].size() > 0)
        m_h_bc_X_gk_ptr[n] = m_bc_X_gk[n].dataPtr();
    };

    amrex::Gpu::HostVector<amrex::Real*>& h_bc_X_gk_ptr () { return m_h_bc_X_gk_ptr; };
    amrex::Real* bc_X_gk_ptr (const int n) { return m_bc_X_gk_ptr[n]; };
    amrex::Gpu::DeviceVector<amrex::Real*>& bc_X_gk_ptr () { return m_bc_X_gk_ptr; };

  private:
    // BC types for each different variable
    std::map< std::string, amrex::Gpu::DeviceVector<int> > m_vel_g_bc_types;
    std::map< std::string, amrex::Gpu::DeviceVector<int> > m_ro_g_bc_types;
    std::map< std::string, amrex::Gpu::DeviceVector<int> > m_T_g_bc_types;
    std::map< std::string, amrex::Gpu::DeviceVector<int> > m_trac_g_bc_types;
    std::map< std::string, amrex::Gpu::DeviceVector<int> > m_X_gk_bc_types;

    // BCs
    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_u_g;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_u_g;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_v_g;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_v_g;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_w_g;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_w_g;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_t_g;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_t_g;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_h_g;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_h_g;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_ro_g;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_ro_g;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_tracer;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_tracer;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_ep_g;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_ep_g;

    amrex::Gpu::HostVector<amrex::Real>   m_h_bc_p_g;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_p_g;

    amrex::Vector< amrex::Gpu::HostVector<amrex::Real> >   m_h_bc_X_gk;
    amrex::Vector< amrex::Gpu::DeviceVector<amrex::Real> > m_bc_X_gk;
    amrex::Gpu::HostVector<amrex::Real*>   m_h_bc_X_gk_ptr;
    amrex::Gpu::DeviceVector<amrex::Real*> m_bc_X_gk_ptr;

    // Lists of BCs applied to the domain extent
    std::array<amrex::Vector<int>,3> m_bc_lo;
    std::array<amrex::Vector<int>,3> m_bc_hi;

    // Flag to solve fluid equations
    int m_delp_dir;
    amrex::RealVect m_delp;

    // Flags for DEM domain extent reinforcement.
    std::array<int,6> m_domain_bc;

    std::bitset<6> m_flow_plane;

    // Vectors storing EB planes for level-set creation
    amrex::Vector<amrex::EB2::PlaneIF> m_flow_planes;
    amrex::Vector<amrex::EB2::PlaneIF> m_wall_planes;

    // Use AMReX LinOpBCType to define the domain extents
    std::array<amrex::LinOpBCType,3> m_ppe_lobc;
    std::array<amrex::LinOpBCType,3> m_ppe_hibc;

    std::array<amrex::LinOpBCType,3> m_diff_vel_lobc;
    std::array<amrex::LinOpBCType,3> m_diff_vel_hibc;

    std::array<amrex::LinOpBCType,3> m_diff_scal_lobc;
    std::array<amrex::LinOpBCType,3> m_diff_scal_hibc;

    std::array<amrex::LinOpBCType,3> m_diff_temperature_lobc;
    std::array<amrex::LinOpBCType,3> m_diff_temperature_hibc;

    std::array<amrex::LinOpBCType,3> m_diff_species_lobc;
    std::array<amrex::LinOpBCType,3> m_diff_species_hibc;

    MFIXEmbeddedBoundaries& m_embedded_boundaries;

    // Data structure for BC data
    amrex::Vector<BC_t> m_bc;
};


struct BC_t {

  // Constructor
  BC_t (MFIXFluidPhase& fluid_in) :
    fluid(fluid_in)
  {}

  // BC type (mi, pi, po, nsw, eb)
  int type;

  const amrex::RealBox* region;

  EB_t eb;

  FLUID_t fluid;

  amrex::Vector<SOLIDS_t> solids;

  const SOLIDS_t* get_solid (const int phase) const
  {
    for (int n(0); n < solids.size(); ++n)
      if (phase == solids[n].phase)
        return &(solids[n]);

    amrex::Abort("Error: solid type not found");
    return nullptr;
  }
};

#endif
