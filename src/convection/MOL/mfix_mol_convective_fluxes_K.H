#ifndef MFIX_MOL_CONVECTIVE_FLUXES_K_H_
#define MFIX_MOL_CONVECTIVE_FLUXES_K_H_

#include <mfix_upwind.H>
#include <mfix_algorithm.H>

#include <AMReX_Slopes_K.H>

using namespace amrex;
using namespace aux;

namespace mol {

void compute_convective_fluxes (const int lev,
                                Box const& domain_box,
                                Box const& xbx,
                                Box const& ybx,
                                Box const& zbx,
                                const int ncomp,
                                const int state_comp,
                                Array4<Real const> const& state,
                                Array4<Real      > const& fx,
                                Array4<Real      > const& fy,
                                Array4<Real      > const& fz,
                                Array4<Real const> const& ep_u_mac,
                                Array4<Real const> const& ep_v_mac,
                                Array4<Real const> const& ep_w_mac,
                                const GpuArray<int, 3> bc_types,
                                std::map<std::string, Gpu::DeviceVector<int>>& state_bcs,
                                Array4<int const> const& bct_ilo,
                                Array4<int const> const& bct_ihi,
                                Array4<int const> const& bct_jlo,
                                Array4<int const> const& bct_jhi,
                                Array4<int const> const& bct_klo,
                                Array4<int const> const& bct_khi)
{

  constexpr int order = 2;
  constexpr Real small_vel = 1.e-10;

  const int* dirichlet_bcs  = (state_bcs["Dirichlet"]).data();
  const int  dirichlet_size = (state_bcs["Dirichlet"]).size();

  /*********************************************************************************
   *                                                                               *
   *                                  x-direction                                  *
   *                                                                               *
   *********************************************************************************/

  const int domain_ilo = domain_box.smallEnd(0);
  const int domain_ihi = domain_box.bigEnd(0);

  bool check_extdir_ilo = (domain_ilo >= xbx.smallEnd(0) and domain_ilo <= xbx.bigEnd(0));
  bool check_extdir_ihi = (domain_ihi >= xbx.smallEnd(0) and domain_ihi <= xbx.bigEnd(0));

  if (check_extdir_ilo or check_extdir_ihi)
  {
    ParallelFor(xbx, ncomp, [fx, ep_u_mac, state, state_comp, order, small_vel,
    check_extdir_ilo, check_extdir_ihi, domain_ilo, domain_ihi,
    bc_types, dirichlet_bcs, dirichlet_size, bct_ilo, bct_ihi]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int* bct_data = bc_types.data();
      const int  bct_size = bc_types.size();

      //
      // West face
      //
      // In the case of MINF       we are using the prescribed Dirichlet value
      // In the case of PINF, POUT we are using the upwind value

      bool valid_bc_ilo = aux::any_of(&bct_data[0], &bct_data[bct_size],
                  aux::is_equal<int>(bct_ilo(domain_ilo-1,j,k,0)));

      bool valid_bc_ihi = aux::any_of(&bct_data[0], &bct_data[bct_size],
                  aux::is_equal<int>(bct_ihi(domain_ihi+1,j,k,0)));

      int comp = n+state_comp;

      Real fx_val(0.);

      if (i <= domain_ilo && valid_bc_ilo) {
        fx_val = state(domain_ilo-1,j,k,comp);

      } else if (i >= domain_ihi+1 && valid_bc_ihi) {
        fx_val = state(domain_ihi+1,j,k,comp);

      } else {

        bool extdir_ilo = check_extdir_ilo && (i-1 <= domain_ilo) &&
          aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                      aux::is_equal<int>(bct_ilo(domain_ilo-1,j,k,0)));

        bool extdir_ihi = check_extdir_ihi && (i   >= domain_ihi) &&
          aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                      aux::is_equal<int>(bct_ihi(domain_ihi+1,j,k,0)));

        Real spls = state(i  ,j,k,comp) - 0.5 * amrex_calc_xslope_extdir
          (i  ,j,k,comp,order,state,extdir_ilo,extdir_ihi,domain_ilo,domain_ihi);

        Real smns = state(i-1,j,k,comp) + 0.5 * amrex_calc_xslope_extdir
          (i-1,j,k,comp,order,state,extdir_ilo,extdir_ihi,domain_ilo,domain_ihi);

        fx_val = upwind(smns, spls, ep_u_mac(i,j,k));

      }

      fx(i,j,k,n) = fx_val * ep_u_mac(i,j,k);

    });

  }
  else
  {
    ParallelFor(xbx, ncomp, [fx, ep_u_mac, state, state_comp, order]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      int comp = n+state_comp;

      Real spls = state(i  ,j,k,comp) - 0.5 * amrex_calc_xslope(i  ,j,k,comp,order,state);
      Real smns = state(i-1,j,k,comp) + 0.5 * amrex_calc_xslope(i-1,j,k,comp,order,state);

      fx(i,j,k,n) = upwind(smns, spls, ep_u_mac(i,j,k));

      fx(i,j,k,n) *= ep_u_mac(i,j,k);
    });
  }


  /*********************************************************************************
   *                                                                               *
   *                                  y-direction                                  *
   *                                                                               *
   *********************************************************************************/

  const int domain_jlo = domain_box.smallEnd(1);
  const int domain_jhi = domain_box.bigEnd(1);

  bool check_extdir_jlo = (domain_jlo >= ybx.smallEnd(1) and domain_jlo <= zbx.bigEnd(1));
  bool check_extdir_jhi = (domain_jhi >= ybx.smallEnd(1) and domain_jhi <= zbx.bigEnd(1));

  if (check_extdir_jlo or check_extdir_jhi)
  {
    ParallelFor(ybx, ncomp, [fy, ep_v_mac, state, state_comp, order, small_vel,
    check_extdir_jlo, check_extdir_jhi, domain_jlo, domain_jhi,
    bc_types, dirichlet_bcs, dirichlet_size, bct_jlo, bct_jhi]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int* bct_data = bc_types.data();
      const int  bct_size = bc_types.size();

      //
      // South face
      //
      // In the case of MINF       we are using the prescribed Dirichlet value
      // In the case of PINF, POUT we are using the upwind value
      bool valid_bc_jlo = aux::any_of(&bct_data[0], &bct_data[bct_size],
                  aux::is_equal<int>(bct_jlo(i,domain_jlo-1,k,0)));

      bool valid_bc_jhi = aux::any_of(&bct_data[0], &bct_data[bct_size],
                  aux::is_equal<int>(bct_jhi(i,domain_jhi+1,k,0)));

      int comp = n+state_comp;

      Real fy_val(0);

      if (j <= domain_jlo && valid_bc_jlo) {
        fy_val = state(i,domain_jlo-1,k,comp);

      } else if (j >= domain_jhi+1 && valid_bc_jhi) {
        fy_val = state(i,domain_jhi+1,k,comp);
    }
      else {

        bool extdir_jlo = check_extdir_jlo && (j-1 <= domain_jlo) &&
          aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                      aux::is_equal<int>(bct_jlo(i,domain_jlo-1,k,0)));

        bool extdir_jhi = check_extdir_jhi && (j   >= domain_jlo) &&
          aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                      aux::is_equal<int>(bct_jhi(i,domain_jhi+1,k,0)));

        Real spls = state(i,j  ,k,comp) - 0.5 * amrex_calc_yslope_extdir
          (i,j  ,k,comp,order,state,extdir_jlo,extdir_jhi,domain_jlo,domain_jhi);

        Real smns = state(i,j-1,k,comp) + 0.5 * amrex_calc_yslope_extdir
          (i,j-1,k,comp,order,state,extdir_jlo,extdir_jhi,domain_jlo,domain_jhi);

        fy_val = upwind(smns, spls, ep_v_mac(i,j,k));
      }

      fy(i,j,k,n) = fy_val * ep_v_mac(i,j,k);

    });
  }
  else
  {
    ParallelFor(ybx, ncomp, [fy, ep_v_mac, state, state_comp, order, small_vel]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      int comp = n+state_comp;

      Real spls = state(i,j  ,k,comp) - 0.5 * amrex_calc_yslope(i,j  ,k,comp,order,state);
      Real smns = state(i,j-1,k,comp) + 0.5 * amrex_calc_yslope(i,j-1,k,comp,order,state);

      fy(i,j,k,n) = upwind(smns, spls, ep_v_mac(i,j,k));

      fy(i,j,k,n) *= ep_v_mac(i,j,k);

    });
  }

  /*********************************************************************************
   *                                                                               *
   *                                  z-direction                                  *
   *                                                                               *
   *********************************************************************************/

  const int domain_klo = domain_box.smallEnd(2);
  const int domain_khi = domain_box.bigEnd(2);

  bool check_extdir_klo = (domain_klo >= xbx.smallEnd(2) and domain_klo <= xbx.bigEnd(2));
  bool check_extdir_khi = (domain_khi >= xbx.smallEnd(2) and domain_khi <= xbx.bigEnd(2));

  if (check_extdir_klo or check_extdir_khi)
  {
    ParallelFor(zbx, ncomp, [fz, ep_w_mac, state, state_comp, order, small_vel,
    check_extdir_klo, check_extdir_khi, domain_klo, domain_khi,
    bc_types, dirichlet_bcs, dirichlet_size, bct_klo, bct_khi]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int* bct_data = bc_types.data();
      const int  bct_size = bc_types.size();

      //
      // Bottom face
      //
      // In the case of MINF       we are using the prescribed Dirichlet value
      // In the case of PINF, POUT we are using the upwind value

      bool valid_bc_klo = aux::any_of(&bct_data[0], &bct_data[bct_size],
                  aux::is_equal<int>(bct_klo(i,j,domain_klo-1,0)));

      bool valid_bc_khi = aux::any_of(&bct_data[0], &bct_data[bct_size],
                  aux::is_equal<int>(bct_khi(i,j,domain_khi+1,0)));

      int comp = n+state_comp;

      Real fz_val(0);

      if (k <= domain_klo && valid_bc_klo) {
        fz_val = state(i,j,domain_klo-1,comp);

      } else if (k >= domain_khi+1 && valid_bc_khi) {
        fz_val = state(i,j,domain_khi+1,comp);
      }
      else {

        bool extdir_klo = check_extdir_klo && (k-1 <= domain_klo) &&
          aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                      aux::is_equal<int>(bct_klo(i,j,domain_klo-1,0)));

        bool extdir_khi = check_extdir_khi && (k   >= domain_khi) &&
          aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                      aux::is_equal<int>(bct_khi(i,j,domain_khi+1,0)));

        Real state_pls = state(i,j,k  ,comp) - 0.5 * amrex_calc_zslope_extdir
          (i,j,k  ,comp,order,state,extdir_klo,extdir_khi,domain_klo,domain_khi);

        Real state_mns = state(i,j,k-1,comp) + 0.5 * amrex_calc_zslope_extdir
          (i,j,k-1,comp,order,state,extdir_klo,extdir_khi,domain_klo,domain_khi);

        fz_val = upwind(state_mns, state_pls, ep_w_mac(i,j,k));
      }

      fz(i,j,k,n) = fz_val * ep_w_mac(i,j,k);

    });
  }
  else
  {
    ParallelFor(zbx, ncomp, [fz, ep_w_mac, state, state_comp, order, small_vel]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      int comp = n+state_comp;

      Real state_pls = state(i,j,k  ,comp) - 0.5 * amrex_calc_zslope(i,j,k  ,comp,order,state);
      Real state_mns = state(i,j,k-1,comp) + 0.5 * amrex_calc_zslope(i,j,k-1,comp,order,state);

      fz(i,j,k,n) = upwind(state_mns, state_pls, ep_w_mac(i,j,k));

      fz(i,j,k,n) *= ep_w_mac(i,j,k);

    });

  }

}

}//end mol namespace
#endif
