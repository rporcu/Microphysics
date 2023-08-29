#ifndef _MFIX_BC_H_
#define _MFIX_BC_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>

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
    MFIXBoundaryConditions (amrex::Vector<amrex::Geometry>& geom,
                            BCList& bc_list,
                            MFIXEmbeddedBoundaries& embedded_boundaries);

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

    void set_epg_bcs (amrex::Real time,
                      const amrex::Vector< amrex::MultiFab* >& ep_g,
                      const int dir_bc);

    void set_density_bc_values (amrex::Real time_in);

    void set_density_bcs (amrex::Real time,
                          amrex::Vector< amrex::MultiFab* > const& ro_g);

    void set_velocity_bc_values (amrex::Real time_in);

    void set_velocity_bcs (amrex::Real time,
                           amrex::Vector< amrex::MultiFab* > const& vel_g,
                           int extrap_dir_bcs);

    void set_eb_velocity_bcs (amrex::Real time,
                              amrex::Vector< amrex::MultiFab* > const& vel_g_on_EB);

    void set_tracer_bc_values (amrex::Real time_in,
                               MFIXFluidPhase& fluid);

    void set_tracer_bcs (amrex::Real time,
                         MFIXFluidPhase& fluid,
                         amrex::Vector< amrex::MultiFab* > const& trac);

    void set_temperature_bc_values (amrex::Real time_in,
                                    MFIXFluidPhase& fluid);

    void set_temperature_bcs (amrex::Real time,
                              MFIXFluidPhase& fluid,
                              amrex::Vector< amrex::MultiFab* > const& T_g);

    void set_enthalpy_bcs (amrex::Real time,
                           MFIXFluidPhase& fluid,
                           amrex::Vector< amrex::MultiFab* > const& T_g);

    void set_eb_temperature_bcs (amrex::Vector< amrex::MultiFab* > const& T_g_on_EB);

    void set_species_bc_values (amrex::Real time_in,
                                MFIXFluidPhase& fluid);

    void set_species_bcs (amrex::Real time,
                          MFIXFluidPhase& fluid,
                          amrex::Vector< amrex::MultiFab* > const& X_gk);

    void set_eb_scalar_bcs (MFIXFluidPhase& fluid,
                            amrex::Vector< amrex::MultiFab* > const& scalars_on_EB,
                            amrex::Vector< amrex::MultiFab* > const& species_on_EB);

    void set_gradp_bcs (const amrex::Box& bx,
                        const int lev,
                        amrex::FArrayBox& gp_fab,
                        amrex::Box& domain);

    void set_vec_bcs (amrex::Real time,
                      amrex::Vector< amrex::MultiFab* > const& mf_in);

    void calc_bc_areas (int const lev,
                        amrex::BoxArray            const& a_grids,
                        amrex::DistributionMapping const& a_dmap,
                        amrex::EBFArrayBoxFactory* factory);

    amrex::Real calc_bc_area (int             const       a_lev,
                              int             const       a_dir,
                              int             const       a_bcv,
                              amrex::Box      const&      a_bx,
                              BC_t            const&      a_bc,
                              amrex::MultiFab const&      a_MF,
                              amrex::EBFArrayBoxFactory*  a_factory);

    amrex::Real calc_eb_bc_area (int             const       a_lev,
                                 int             const       a_bcv,
                                 BC_t            const&      a_bc,
                                 amrex::MultiFab const&      a_MF,
                                 amrex::EBFArrayBoxFactory*  a_factory);

    template <typename F1>
    AMREX_FORCE_INLINE
    void set_bcs (amrex::Real /*time*/,
                  F1 bcs_function,
                  amrex::Vector< amrex::MultiFab* > const& mf_in)
    {
      const int nlev = mf_in.size();

      for (int lev = 0; lev < nlev; lev++) {

        const int ncomp = mf_in[lev]->nComp();

        using BctArr = amrex::Array4<const int>;

        BctArr const& bct_ilo = m_bc_list.bc_ilo[lev]->const_array();
        BctArr const& bct_jlo = m_bc_list.bc_jlo[lev]->const_array();
        BctArr const& bct_klo = m_bc_list.bc_klo[lev]->const_array();

        amrex::Array<BctArr const*, 3> bct_lo_arr = {&bct_ilo, &bct_jlo, &bct_klo};

        BctArr const& bct_ihi = m_bc_list.bc_ihi[lev]->const_array();
        BctArr const& bct_jhi = m_bc_list.bc_jhi[lev]->const_array();
        BctArr const& bct_khi = m_bc_list.bc_khi[lev]->const_array();

        amrex::Array<BctArr const*, 3> bct_hi_arr = {&bct_ihi, &bct_jhi, &bct_khi};

        amrex::Array<amrex::Array<BctArr const*, 3>*, 2> bct_vect = {&bct_lo_arr, &bct_hi_arr};

        amrex::Box domain(m_geom[lev].Domain());

        amrex::IntVect dom_lo(domain.loVect());
        amrex::IntVect dom_hi(domain.hiVect());

        amrex::Array<amrex::IntVect*, 2> dom_vect = {&dom_lo, &dom_hi};

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(*mf_in[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

          amrex::FArrayBox& mf_fab = (*mf_in[lev])[mfi];

          const amrex::EBFArrayBox& mf_EB_fab = static_cast<amrex::EBFArrayBox const&>(mf_fab);
          const amrex::EBCellFlagFab& flags_fab = mf_EB_fab.getEBCellFlagFab();

          amrex::IntVect mf_lo(mf_fab.loVect());
          amrex::IntVect mf_hi(mf_fab.hiVect());

          amrex::Array<amrex::IntVect*, 2> mf_vect = {&mf_lo, &mf_hi};

          amrex::Array4<amrex::Real> const& mf_arr = mf_fab.array();
          amrex::Array4<const amrex::EBCellFlag> const& flags_arr = flags_fab.const_array();

          for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {

            amrex::IntVect bx_lo_hi(mf_hi);
            amrex::IntVect bx_hi_lo(mf_lo);

            bx_lo_hi[dir] = dom_lo[dir]-1;
            bx_hi_lo[dir] = dom_hi[dir]+1;

            amrex::Box bx_lo(mf_lo, bx_lo_hi);
            amrex::Box bx_hi(bx_hi_lo, mf_hi);

            amrex::Array<amrex::Box*, 2> bx_vect = {&bx_lo, &bx_hi};

            for (int face(-1); face <= 1; face += 2) {

              const int face_idx = (face+1) / 2;

              const int face_flag = amrex::max(0,
                  face*((*(mf_vect[face_idx]))[dir]-((*(dom_vect[face_idx]))[dir])));

              if (face_flag > 0) {

                amrex::Box bx = *(bx_vect[face_idx]);

                amrex::Array4<const int> const& bct_arr = *((*(bct_vect[face_idx]))[dir]);
                amrex::IntVect dom = *(dom_vect[face_idx]);

                amrex::ParallelFor(bx, ncomp, [bct_arr,dom,mf_arr,flags_arr,
                    bcs_function,dir,face]
                  AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                  amrex::GpuArray<int, AMREX_SPACEDIM> ijk = {i,j,k};
                  amrex::IntVect dom_ijk(ijk.data());
                  dom_ijk[dir] = dom[dir]+face;

                  const int bct = bct_arr(dom_ijk,0);
                  const int bcv = bct_arr(dom_ijk,1);

                  dom_ijk[dir] = dom[dir];

                  bcs_function(bct, bcv, {i,j,k}, dom_ijk, n, mf_arr, flags_arr, dir);
                });
              }
            }
          }

        } // end MFIter loop

        mf_in[lev]->FillBoundary(m_geom[lev].periodicity());

        amrex::EB_set_covered(*mf_in[lev], 0, mf_in[lev]->nComp(), mf_in[lev]->nGrow(),
            m_covered_val);
      }
    }

    template <typename F1>
    AMREX_FORCE_INLINE
    void set_bcs_2D (amrex::Real /*time*/,
                     F1 bcs_function,
                     amrex::Vector< amrex::MultiFab* > const& mf_in)
    {
      const int nlev = mf_in.size();

      for (int lev = 0; lev < nlev; lev++) {

        const int ncomp = mf_in[lev]->nComp();

        using BctArr = amrex::Array4<const int>;

        BctArr const& bct_ilo = m_bc_list.bc_ilo[lev]->const_array();
        BctArr const& bct_jlo = m_bc_list.bc_jlo[lev]->const_array();
        BctArr const& bct_klo = m_bc_list.bc_klo[lev]->const_array();

        amrex::Array<BctArr const*, 3> bct_lo_arr = {&bct_ilo, &bct_jlo, &bct_klo};

        BctArr const& bct_ihi = m_bc_list.bc_ihi[lev]->const_array();
        BctArr const& bct_jhi = m_bc_list.bc_jhi[lev]->const_array();
        BctArr const& bct_khi = m_bc_list.bc_khi[lev]->const_array();

        amrex::Array<BctArr const*, 3> bct_hi_arr = {&bct_ihi, &bct_jhi, &bct_khi};

        amrex::Array<amrex::Array<BctArr const*, 3>*, 2> bct_vect = {&bct_lo_arr, &bct_hi_arr};

        amrex::Box domain(m_geom[lev].Domain());

        amrex::IntVect dom_lo(domain.loVect());
        amrex::IntVect dom_hi(domain.hiVect());

        amrex::Array<amrex::IntVect*, 2> dom_vect = {&dom_lo, &dom_hi};

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(*mf_in[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

          amrex::FArrayBox& mf_fab = (*mf_in[lev])[mfi];

          amrex::IntVect mf_lo(mf_fab.loVect());
          amrex::IntVect mf_hi(mf_fab.hiVect());

          amrex::Array<amrex::IntVect*, 2> mf_vect = {&mf_lo, &mf_hi};

          amrex::Array4<amrex::Real> const& mf_arr = mf_fab.array();

          for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {

            amrex::IntVect bx_lo_lo(mf_lo), bx_lo_hi(mf_hi);
            amrex::IntVect bx_hi_lo(mf_lo), bx_hi_hi(mf_hi);

            bx_lo_lo[dir] = dom_lo[dir]-1;
            bx_lo_hi[dir] = dom_lo[dir]-1;
            bx_hi_lo[dir] = dom_hi[dir]+1;
            bx_hi_hi[dir] = dom_hi[dir]+1;

            amrex::Box bx_lo(bx_lo_lo, bx_lo_hi);
            amrex::Box bx_hi(bx_hi_lo, bx_hi_hi);

            amrex::Array<amrex::Box*, 2> bx_vect = {&bx_lo, &bx_hi};

            for (int face(-1); face <= 1; face += 2) {

              const int face_idx = (face+1) / 2;

              const int face_flag = amrex::max(0,
                  face*((*(mf_vect[face_idx]))[dir] - (*(dom_vect[face_idx]))[dir]));

              if (face_flag > 0)
              {
                amrex::Box bx = *(bx_vect[face_idx]);

                amrex::Array4<const int> const& bct_arr = *((*(bct_vect[face_idx]))[dir]);
                amrex::IntVect dom = *(dom_vect[face_idx]);

                amrex::ParallelFor(bx, ncomp, [bct_arr,dom,mf_arr,dir,face,
                    bcs_function]
                  AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                  amrex::GpuArray<int, AMREX_SPACEDIM> ijk = {i,j,k};
                  amrex::IntVect dom_ijk(ijk.data());
                  dom_ijk[dir] = dom[dir]+face;

                  const int bct = bct_arr(dom_ijk,0);

                  dom_ijk[dir] = dom[dir];

                  amrex::IntVect near_ijk(ijk.data());
                  near_ijk[dir] -= face;

                  bcs_function(bct, {i,j,k}, dom_ijk, near_ijk, n, mf_arr);
                });
              }
            }
          }

        } // end MFIter loop

        mf_in[lev]->FillBoundary(m_geom[lev].periodicity());

        amrex::EB_set_covered(*mf_in[lev], 0, mf_in[lev]->nComp(), mf_in[lev]->nGrow(),
            m_covered_val);
      }
    }

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
    const amrex::Vector<amrex::EB2::PlaneIF>& flow_planes () const { return m_flow_planes; }
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

    amrex::Real& h_bc_u_g (const int bcv) { return m_h_bc_u_g[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_u_g () { return m_h_bc_u_g; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_u_g () { return m_bc_u_g; }

    amrex::Real& h_bc_v_g (const int bcv) { return m_h_bc_v_g[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_v_g () { return m_h_bc_v_g; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_v_g () { return m_bc_v_g; }

    amrex::Real& h_bc_w_g (const int bcv) { return m_h_bc_w_g[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_w_g () { return m_h_bc_w_g; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_w_g () { return m_bc_w_g; }

    amrex::Real& h_bc_t_g (const int bcv) { return m_h_bc_t_g[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_t_g () { return m_h_bc_t_g; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_t_g () { return m_bc_t_g; }

    amrex::Real& h_bc_h_g (const int bcv) { return m_h_bc_h_g[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_h_g () { return m_h_bc_h_g; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_h_g () { return m_bc_h_g; }

    amrex::Real& h_bc_ro_g (const int bcv) { return m_h_bc_ro_g[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_ro_g () { return m_h_bc_ro_g; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_ro_g () { return m_bc_ro_g; }

    amrex::Real& h_bc_tracer (const int bcv) { return m_h_bc_tracer[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_tracer () { return m_h_bc_tracer; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_tracer () { return m_bc_tracer; }

    amrex::Real& h_bc_ep_g (const int bcv) { return m_h_bc_ep_g[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_ep_g () { return m_h_bc_ep_g; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_ep_g () { return m_bc_ep_g; }

    amrex::Real& h_bc_p_g (const int bcv) { return m_h_bc_p_g[bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_p_g () { return m_h_bc_p_g; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_p_g () { return m_bc_p_g; }

    amrex::Real& h_bc_X_gk (const int n, const int bcv) { return m_h_bc_X_gk[n][bcv]; }
    amrex::Gpu::HostVector<amrex::Real>& h_bc_X_gk (const int n) { return m_h_bc_X_gk[n]; }
    amrex::Vector<amrex::Gpu::HostVector<amrex::Real>>& h_bc_X_gk () { return m_h_bc_X_gk; }
    amrex::Gpu::DeviceVector<amrex::Real>& bc_X_gk (const int n) { return m_bc_X_gk[n]; }
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>>& bc_X_gk () { return m_bc_X_gk; }
    amrex::Real* h_bc_X_gk_ptr (const int n) { return m_h_bc_X_gk_ptr[n]; }

    void set_h_bc_X_gk_ptr (const int n)
    {
      if (m_bc_X_gk[n].size() > 0)
        m_h_bc_X_gk_ptr[n] = m_bc_X_gk[n].dataPtr();
    }

    amrex::Gpu::HostVector<amrex::Real*>& h_bc_X_gk_ptr () { return m_h_bc_X_gk_ptr; }
    amrex::Real* bc_X_gk_ptr (const int n) { return m_bc_X_gk_ptr[n]; }
    amrex::Gpu::DeviceVector<amrex::Real*>& bc_X_gk_ptr () { return m_bc_X_gk_ptr; }

    amrex::Real& get_bc_area (const int bcv);

    amrex::Real& get_bc_fab_area (const int a_bcv, std::pair<int,int> a_index);

    int& get_dir (const int bcv) { return m_dir[bcv]; }

  private:
    amrex::Vector<amrex::Geometry>& m_geom;

    BCList& m_bc_list;

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

    // Domain face of BC
    // lo x/y/z: 0, 2, 4
    // hi x/y/z: 1: 3, 5
    amrex::Vector<int> m_dir;

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

    amrex::Vector<amrex::Real> m_area;

    amrex::Vector<std::map<std::pair<int,int>, amrex::Real>> m_fab_area;

    // Data structure for BC data
    amrex::Vector<BC_t> m_bc;

    const amrex::Real m_covered_val = 1.e+40;
};


struct BC_t {

  // Constructor
  BC_t (MFIXFluidPhase* fluid_in) :
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
