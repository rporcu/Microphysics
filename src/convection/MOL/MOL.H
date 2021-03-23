//
//  MOL.H
//

#ifndef MOL_H
#define MOL_H

#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>

#include <AMReX_MultiCutFab.H>


namespace mol {

  void predict_vels_on_faces (int lev,
                              amrex::MultiFab& ep_u_mac,
                              amrex::MultiFab& ep_v_mac,
                              amrex::MultiFab& ep_w_mac,
                              const amrex::MultiFab& vel,
                              std::map<std::string, amrex::Gpu::DeviceVector<int>>& bc_types,
                              amrex::Array4<int const> const& bc_ilo,
                              amrex::Array4<int const> const& bc_ihi,
                              amrex::Array4<int const> const& bc_jlo,
                              amrex::Array4<int const> const& bc_jhi,
                              amrex::Array4<int const> const& bc_klo,
                              amrex::Array4<int const> const& bc_khi,
                              amrex::EBFArrayBoxFactory const* ebfact,
                              amrex::Vector<amrex::Geometry> geom);

  void compute_convective_fluxes (const int lev,
                                  amrex::Vector< amrex::MultiFab* >& a_fx,
                                  amrex::Vector< amrex::MultiFab* >& a_fy,
                                  amrex::Vector< amrex::MultiFab* >& a_fz,
                                  amrex::Vector< amrex::MultiFab* > const& state,
                                  const int state_comp,
                                  const int n_comp,
                                  amrex::Vector< amrex::MultiFab* > const& ep_u_mac,
                                  amrex::Vector< amrex::MultiFab* > const& ep_v_mac,
                                  amrex::Vector< amrex::MultiFab* > const& ep_w_mac,
                                  const int nghosts,
                                  const amrex::Real covered_val,
                                  const amrex::GpuArray<int, 3> bc_types,
                                  std::map<std::string, amrex::Gpu::DeviceVector<int>>& state_bcs,
                                  amrex::Array4<int const> const& bct_ilo,
                                  amrex::Array4<int const> const& bct_ihi,
                                  amrex::Array4<int const> const& bct_jlo,
                                  amrex::Array4<int const> const& bct_jhi,
                                  amrex::Array4<int const> const& bct_klo,
                                  amrex::Array4<int const> const& bct_khi,
                                  amrex::EBFArrayBoxFactory const* ebfact,
                                  amrex::Vector<amrex::Geometry> geom);

  void compute_convective_term (const int lev,
                                amrex::Vector<amrex::BoxArray> grids,
                                amrex::DistributionMapping const& dmap,
                                amrex::Vector< amrex::MultiFab* >& conv_u_in,
                                amrex::Vector< amrex::MultiFab* >& conv_s_in,
                                amrex::Vector< amrex::MultiFab* >& conv_X_in,
                                const bool advect_density,
                                const bool advect_enthalpy,
                                const bool advect_tracer,
                                const bool advect_species,
                                amrex::Vector< amrex::MultiFab* >& fx,
                                amrex::Vector< amrex::MultiFab* >& fy,
                                amrex::Vector< amrex::MultiFab* >& fz,
                                amrex::Vector< amrex::MultiFab* > const& vel_in,
                                amrex::Vector< amrex::MultiFab* > const& ep_g_in,
                                amrex::Vector< amrex::MultiFab* > const& ep_u_mac,
                                amrex::Vector< amrex::MultiFab* > const& ep_v_mac,
                                amrex::Vector< amrex::MultiFab* > const& ep_w_mac,
                                amrex::Vector< amrex::MultiFab* > const& ro_g_in,
                                amrex::Vector< amrex::MultiFab* > const& h_g_in,
                                amrex::Vector< amrex::MultiFab* > const& trac_in,
                                amrex::Vector< amrex::MultiFab* > const& X_gk_in,
                                const amrex::GpuArray<int, 3> bc_types,
                                std::map<std::string, amrex::Gpu::DeviceVector<int>>& velocity_bcs,
                                std::map<std::string, amrex::Gpu::DeviceVector<int>>& density_bcs,
                                std::map<std::string, amrex::Gpu::DeviceVector<int>>& enthalpy_bcs,
                                std::map<std::string, amrex::Gpu::DeviceVector<int>>& tracer_bcs,
                                std::map<std::string, amrex::Gpu::DeviceVector<int>>& species_bcs,
                                amrex::Array4<int const> const& bct_ilo,
                                amrex::Array4<int const> const& bct_ihi,
                                amrex::Array4<int const> const& bct_jlo,
                                amrex::Array4<int const> const& bct_jhi,
                                amrex::Array4<int const> const& bct_klo,
                                amrex::Array4<int const> const& bct_khi,
                                amrex::EBFArrayBoxFactory const* ebfact,
                                amrex::Vector<amrex::Geometry> geom);

  void compute_convective_fluxes (amrex::Box const& bx,
                                  const int ncomp,
                                  const int state_comp,
                                  amrex::Array4<amrex::Real const> const& state,
                                  amrex::Array4<amrex::Real      > const& fx,
                                  amrex::Array4<amrex::Real      > const& fy,
                                  amrex::Array4<amrex::Real      > const& fz,
                                  amrex::Array4<amrex::Real const> const& ep_u_mac,
                                  amrex::Array4<amrex::Real const> const& ep_v_mac,
                                  amrex::Array4<amrex::Real const> const& ep_w_mac,
                                  const amrex::GpuArray<int, 3> bc_types,
                                  std::map<std::string, amrex::Gpu::DeviceVector<int>>& state_bcs,
                                  amrex::Array4<int const> const& bct_ilo,
                                  amrex::Array4<int const> const& bct_ihi,
                                  amrex::Array4<int const> const& bct_jlo,
                                  amrex::Array4<int const> const& bct_jhi,
                                  amrex::Array4<int const> const& bct_klo,
                                  amrex::Array4<int const> const& bct_khi,
                                  amrex::Geometry& geom);

  void compute_convective_fluxes_eb (amrex::Box const& bx,
                                     const int ncomp,
                                     const int state_comp,
                                     amrex::Array4<amrex::Real const> const& state,
                                     amrex::Array4<amrex::Real      > const& fx,
                                     amrex::Array4<amrex::Real      > const& fy,
                                     amrex::Array4<amrex::Real      > const& fz,
                                     amrex::Array4<amrex::Real const> const& ep_u_mac,
                                     amrex::Array4<amrex::Real const> const& ep_v_mac,
                                     amrex::Array4<amrex::Real const> const& ep_w_mac,
                                     amrex::Array4<amrex::EBCellFlag const> const& flag,
                                     amrex::Array4<amrex::Real const> const& fcx,
                                     amrex::Array4<amrex::Real const> const& fcy,
                                     amrex::Array4<amrex::Real const> const& fcz,
                                     amrex::Array4<amrex::Real const> const& ccc,
                                     const amrex::GpuArray<int, 3> bc_types,
                                     std::map<std::string, amrex::Gpu::DeviceVector<int>>& state_bcs,
                                     amrex::Array4<int const> const& bct_ilo,
                                     amrex::Array4<int const> const& bct_ihi,
                                     amrex::Array4<int const> const& bct_jlo,
                                     amrex::Array4<int const> const& bct_jhi,
                                     amrex::Array4<int const> const& bct_klo,
                                     amrex::Array4<int const> const& bct_khi,
                                     amrex::Geometry& geom);

  void compute_convective_rate ( amrex::Box const& bx, int ncomp, int scomp,
                                 amrex::Array4<amrex::Real> const& dUdt,
                                 amrex::Array4<amrex::Real const> const& fx,
                                 amrex::Array4<amrex::Real const> const& fy,
                                 amrex::Array4<amrex::Real const> const& fz,
                                 amrex::Geometry& geom);

  void compute_convective_rate_eb (amrex::Box const& bx, int ncomp, int scomp,
                                   amrex::Array4<amrex::Real> const& dUdt,
                                   amrex::Array4<amrex::Real const> const& fx,
                                   amrex::Array4<amrex::Real const> const& fy,
                                   amrex::Array4<amrex::Real const> const& fz,
                                   amrex::Array4<amrex::EBCellFlag const> const& flag,
                                   amrex::Array4<amrex::Real const> const& vfrac,
                                   amrex::Array4<amrex::Real const> const& apx,
                                   amrex::Array4<amrex::Real const> const& apy,
                                   amrex::Array4<amrex::Real const> const& apz,
                                   amrex::Geometry& geom);

} // namespace mol

#endif /* MOL_H */
