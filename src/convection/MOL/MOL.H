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
                              amrex::MultiFab& ep_face_x,
                              amrex::MultiFab& ep_face_y,
                              amrex::MultiFab& ep_face_z,
                              amrex::MultiFab& vel,
                              std::map<std::string, amrex::Gpu::DeviceVector<int>>& bc_types,
                              amrex::Array4<int const> const& bc_ilo,
                              amrex::Array4<int const> const& bc_ihi,
                              amrex::Array4<int const> const& bc_jlo,
                              amrex::Array4<int const> const& bc_jhi,
                              amrex::Array4<int const> const& bc_klo,
                              amrex::Array4<int const> const& bc_khi,
                              amrex::EBFArrayBoxFactory const* ebfact,
                              amrex::Vector<amrex::Geometry> geom);

  void mfix_compute_fluxes (const int lev,
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
                            const amrex::GpuArray<int, 2> bc_types,
                            amrex::Array4<int const> const& bct_ilo,
                            amrex::Array4<int const> const& bct_ihi,
                            amrex::Array4<int const> const& bct_jlo,
                            amrex::Array4<int const> const& bct_jhi,
                            amrex::Array4<int const> const& bct_klo,
                            amrex::Array4<int const> const& bct_khi,
                            amrex::EBFArrayBoxFactory const* ebfact,
                            amrex::Vector<amrex::Geometry> geom);

} // namespace mol

#endif /* MOL_H */
