#ifndef MFIX_MOL_CONVECTIVE_FLUXES_EB_K_H_
#define MFIX_MOL_CONVECTIVE_FLUXES_EB_K_H_

#include <mfix_upwind.H>
#include <mfix_algorithm.H>
#include <AMReX_EB_slopes_K.H>

using namespace aux;
using namespace amrex;

namespace mol {

//
// Compute the three components of the convection term when we have embedded
// boundaries
//
void compute_convective_fluxes_eb (const int lev,
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
                                   Array4<EBCellFlag const> const& flag,
                                   Array4<Real const> const& fcx,
                                   Array4<Real const> const& fcy,
                                   Array4<Real const> const& fcz,
                                   Array4<Real const> const& ccc,
                                   const GpuArray<int, 3> bc_types,
                                   std::map<std::string, Gpu::DeviceVector<int>>& state_bcs,
                                   Array4<int const> const& bct_ilo,
                                   Array4<int const> const& bct_ihi,
                                   Array4<int const> const& bct_jlo,
                                   Array4<int const> const& bct_jhi,
                                   Array4<int const> const& bct_klo,
                                   Array4<int const> const& bct_khi)
{
  constexpr Real my_huge = 1.e200;

  const int* dirichlet_bcs  = (state_bcs["Dirichlet"]).data();
  const int  dirichlet_size = (state_bcs["Dirichlet"]).size();

  const int domain_ilo = domain_box.smallEnd(0);
  const int domain_ihi = domain_box.bigEnd(0);
  const int domain_jlo = domain_box.smallEnd(1);
  const int domain_jhi = domain_box.bigEnd(1);
  const int domain_klo = domain_box.smallEnd(2);
  const int domain_khi = domain_box.bigEnd(2);

  const bool check_extdir_ilo = (domain_ilo >= xbx.smallEnd(0) and domain_ilo <= xbx.bigEnd(0));
  const bool check_extdir_ihi = (domain_ihi >= xbx.smallEnd(0) and domain_ihi <= xbx.bigEnd(0));
  const bool check_extdir_jlo = (domain_jlo >= ybx.smallEnd(1) and domain_jlo <= ybx.bigEnd(1));
  const bool check_extdir_jhi = (domain_jhi >= ybx.smallEnd(1) and domain_jhi <= ybx.bigEnd(1));
  const bool check_extdir_klo = (domain_klo >= zbx.smallEnd(2) and domain_klo <= zbx.bigEnd(2));
  const bool check_extdir_khi = (domain_khi >= zbx.smallEnd(2) and domain_khi <= zbx.bigEnd(2));

  const bool check_extdir = check_extdir_ilo or check_extdir_ihi or
                            check_extdir_jlo or check_extdir_jhi or
                            check_extdir_klo or check_extdir_khi;
  /*********************************************************************************
   *                                                                               *
   *                                  x-direction                                  *
   *                                                                               *
   *********************************************************************************/

  if (check_extdir)
  {
    ParallelFor(xbx,ncomp,[fx,ep_u_mac,state,state_comp,ccc,fcx,fcy,fcz,flag,
    check_extdir_ilo, check_extdir_ihi, domain_ilo, domain_ihi, bct_ilo, bct_ihi,
    check_extdir_jlo, check_extdir_jhi, domain_jlo, domain_jhi, bct_jlo, bct_jhi,
    check_extdir_klo, check_extdir_khi, domain_klo, domain_khi, bct_klo, bct_khi,
    bc_types, dirichlet_bcs, dirichlet_size]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      Real sx;

      if( flag(i,j,k).isConnected(-1,0,0) ) {

        const int* bct_data = bc_types.data();
        const int  bct_size = bc_types.size();

        bool valid_bc_ilo = aux::any_of(&bct_data[0], &bct_data[bct_size],
                    aux::is_equal<int>(bct_ilo(domain_ilo-1,j,k,0)));

        bool valid_bc_ihi = aux::any_of(&bct_data[0], &bct_data[bct_size],
                    aux::is_equal<int>(bct_ihi(domain_ihi+1,j,k,0)));

        int comp = n+state_comp;

        if (i <= domain_ilo && valid_bc_ilo) {
          sx = state(domain_ilo-1,j,k,comp);

        } else if (i >= domain_ihi+1 && valid_bc_ihi) {
          sx = state(domain_ihi+1,j,k,comp);

        } else {

          bool extdir_ilo = check_extdir_ilo && (i-1 <= domain_ilo) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_ilo(domain_ilo-1,j,k,0)));

          bool extdir_ihi = check_extdir_ihi && (i   >= domain_ihi) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_ihi(domain_ihi+1,j,k,0)));

          bool extdir_jlo = check_extdir_jlo && (j-1 <= domain_jlo) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_jlo(i,domain_jlo-1,k,0)));

          bool extdir_jhi = check_extdir_jhi && (j   >= domain_jlo) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_jhi(i,domain_jhi+1,k,0)));

          bool extdir_klo = check_extdir_klo && (k-1 <= domain_klo) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_klo(i,j,domain_klo-1,0)));

          bool extdir_khi = check_extdir_khi && (k   >= domain_khi) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_khi(i,j,domain_khi+1,0)));

          const Real state_pls = state(i,j,k,state_comp+n);
          const Real state_mns = state(i-1,j,k,state_comp+n);

          const Real cc_umax = amrex::max(state_pls, state_mns);
          const Real cc_umin = amrex::min(state_pls, state_mns);

          const Real yf = fcx(i,j,k,0);
          const Real zf = fcx(i,j,k,1);

          Real delta_x = .5 + ccc(i,j,k,0);
          Real delta_y = yf - ccc(i,j,k,1);
          Real delta_z = zf - ccc(i,j,k,2);

          const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i  ,j,k,comp,
                                        state, ccc, fcx, fcy, fcz, flag,
                                        extdir_ilo, extdir_jlo, extdir_klo,
                                        extdir_ihi, extdir_jhi, extdir_khi,
                                        domain_ilo, domain_jlo, domain_klo,
                                        domain_ihi, domain_jhi, domain_khi);

          Real upls = state_pls - delta_x * slopes_eb_hi[0]
                                + delta_y * slopes_eb_hi[1]
                                + delta_z * slopes_eb_hi[2];

          upls = amrex::max( amrex::min(upls, cc_umax), cc_umin );

          delta_x = .5 - ccc(i-1,j,k,0);
          delta_y = yf - ccc(i-1,j,k,1);
          delta_z = zf - ccc(i-1,j,k,2);

          const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i-1,j,k,comp,
                                        state, ccc, fcx, fcy, fcz, flag,
                                        extdir_ilo, extdir_jlo, extdir_klo,
                                        extdir_ihi, extdir_jhi, extdir_khi,
                                        domain_ilo, domain_jlo, domain_klo,
                                        domain_ihi, domain_jhi, domain_khi);

          Real umns = state_mns + delta_x * slopes_eb_lo[0]
                                + delta_y * slopes_eb_lo[1]
                                + delta_z * slopes_eb_lo[2];

          umns = amrex::max( amrex::min(umns, cc_umax), cc_umin );

          sx = upwind(umns, upls, ep_u_mac(i,j,k));

        }
      }
      else {
        sx = my_huge;
      }

      fx(i,j,k,n) = ep_u_mac(i,j,k) * sx;

    });

  }
  else
  {
    ParallelFor(xbx,ncomp,[fx,ep_u_mac,state,state_comp,ccc,fcx,fcy,fcz,flag]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      Real sx;

      if( flag(i,j,k).isConnected(-1,0,0) )
      {

        int comp = n+state_comp;

        const Real state_pls = state(i,j,k,state_comp+n);
        const Real state_mns = state(i-1,j,k,state_comp+n);

        const Real cc_umax = amrex::max(state_pls, state_mns);
        const Real cc_umin = amrex::min(state_pls, state_mns);

        const Real yf = fcx(i,j,k,0);
        const Real zf = fcx(i,j,k,1);

        Real delta_x = .5 + ccc(i,j,k,0);
        Real delta_y = yf - ccc(i,j,k,1);
        Real delta_z = zf - ccc(i,j,k,2);

        const auto& slopes_eb_hi = amrex_lim_slopes_eb(i  ,j,k,comp, state, ccc, fcx, fcy, fcz, flag);

        Real upls = state_pls - delta_x * slopes_eb_hi[0]
                              + delta_y * slopes_eb_hi[1]
                              + delta_z * slopes_eb_hi[2];

        upls = amrex::max( amrex::min(upls, cc_umax), cc_umin );

        delta_x = .5 - ccc(i-1,j,k,0);
        delta_y = yf - ccc(i-1,j,k,1);
        delta_z = zf - ccc(i-1,j,k,2);

        const auto& slopes_eb_lo = amrex_lim_slopes_eb(i-1,j,k,comp, state, ccc, fcx, fcy, fcz, flag);

        Real umns = state_mns + delta_x * slopes_eb_lo[0]
                              + delta_y * slopes_eb_lo[1]
                              + delta_z * slopes_eb_lo[2];

        umns = amrex::max( amrex::min(umns, cc_umax), cc_umin );

        sx = upwind(umns, upls, ep_u_mac(i,j,k));

      }
      else {
        sx = my_huge;
      }

      fx(i,j,k,n) = ep_u_mac(i,j,k) * sx;

    });

  }



  /*********************************************************************************
   *                                                                               *
   *                                  y-direction                                  *
   *                                                                               *
   *********************************************************************************/


  if (check_extdir)
  {
    ParallelFor(ybx,ncomp,[fy,ep_v_mac,state,state_comp,ccc,fcx,fcy,fcz,flag,
    check_extdir_ilo, check_extdir_ihi, domain_ilo, domain_ihi, bct_ilo, bct_ihi,
    check_extdir_jlo, check_extdir_jhi, domain_jlo, domain_jhi, bct_jlo, bct_jhi,
    check_extdir_klo, check_extdir_khi, domain_klo, domain_khi, bct_klo, bct_khi,
    bc_types, dirichlet_bcs, dirichlet_size]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {

      Real sy;

      if( flag(i,j,k).isConnected(0,-1,0) ) {

        const int* bct_data = bc_types.data();
        const int  bct_size = bc_types.size();

        bool valid_bc_jlo = aux::any_of(&bct_data[0], &bct_data[bct_size],
                    aux::is_equal<int>(bct_jlo(i,domain_jlo-1,k,0)));

        bool valid_bc_jhi = aux::any_of(&bct_data[0], &bct_data[bct_size],
                    aux::is_equal<int>(bct_jhi(i,domain_jhi+1,k,0)));

        int comp = n+state_comp;

        if  (j <= domain_jlo && valid_bc_jlo) {
          sy = state(i,domain_jlo-1,k,comp);

        } else if (j >= domain_jhi+1 && valid_bc_jhi) {
          sy = state(i,domain_jhi+1,k,comp);

        } else {

          bool extdir_ilo = check_extdir_ilo && (i-1 <= domain_ilo) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_ilo(domain_ilo-1,j,k,0)));

          bool extdir_ihi = check_extdir_ihi && (i   >= domain_ihi) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_ihi(domain_ihi+1,j,k,0)));

          bool extdir_jlo = check_extdir_jlo && (j-1 <= domain_jlo) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_jlo(i,domain_jlo-1,k,0)));

          bool extdir_jhi = check_extdir_jhi && (j   >= domain_jlo) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_jhi(i,domain_jhi+1,k,0)));

          bool extdir_klo = check_extdir_klo && (k-1 <= domain_klo) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_klo(i,j,domain_klo-1,0)));

          bool extdir_khi = check_extdir_khi && (k   >= domain_khi) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_khi(i,j,domain_khi+1,0)));

          const Real state_pls = state(i,j  ,k,state_comp+n);
          const Real state_mns = state(i,j-1,k,state_comp+n);

          const Real cc_umax = amrex::max(state_pls, state_mns);
          const Real cc_umin = amrex::min(state_pls, state_mns);

          const Real xf = fcy(i,j,k,0);
          const Real zf = fcy(i,j,k,1);

          Real delta_y = .5 + ccc(i,j,k,1);
          Real delta_x = xf - ccc(i,j,k,0);
          Real delta_z = zf - ccc(i,j,k,2);

          const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j  ,k,comp,
                                        state, ccc, fcx, fcy, fcz, flag,
                                        extdir_ilo, extdir_jlo, extdir_klo,
                                        extdir_ihi, extdir_jhi, extdir_khi,
                                        domain_ilo, domain_jlo, domain_klo,
                                        domain_ihi, domain_jhi, domain_khi);

          Real vpls = state_pls - delta_y * slopes_eb_hi[1]
                                + delta_x * slopes_eb_hi[0]
                                + delta_z * slopes_eb_hi[2];

          vpls = amrex::max( amrex::min(vpls, cc_umax), cc_umin );

          delta_y = .5 - ccc(i,j-1,k,1);
          delta_x = xf - ccc(i,j-1,k,0);
          delta_z = zf - ccc(i,j-1,k,2);

          const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i,j-1,k,comp,
                                        state, ccc, fcx, fcy, fcz, flag,
                                        extdir_ilo, extdir_jlo, extdir_klo,
                                        extdir_ihi, extdir_jhi, extdir_khi,
                                        domain_ilo, domain_jlo, domain_klo,
                                        domain_ihi, domain_jhi, domain_khi);

          Real vmns = state_mns + delta_y * slopes_eb_lo[1]
                                + delta_x * slopes_eb_lo[0]
                                + delta_z * slopes_eb_lo[2];

          vmns = amrex::max( amrex::min(vmns, cc_umax), cc_umin );

          sy = upwind(vmns, vpls, ep_v_mac(i,j,k));
        }
      }
      else {
        sy = my_huge;
      }

      fy(i,j,k,n) = ep_v_mac(i,j,k) * sy;

    });

  }
  else
  {
    ParallelFor(ybx,ncomp,[fy,ep_v_mac,state,state_comp,ccc,fcx,fcy,fcz,flag,my_huge]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {

      Real sy;

      if( flag(i,j,k).isConnected(0,-1,0) ) {

        int comp = n+state_comp;

        const Real state_pls = state(i,j  ,k,state_comp+n);
        const Real state_mns = state(i,j-1,k,state_comp+n);

        const Real cc_umax = amrex::max(state_pls, state_mns);
        const Real cc_umin = amrex::min(state_pls, state_mns);

        const Real xf = fcy(i,j,k,0);
        const Real zf = fcy(i,j,k,1);

        Real delta_y = .5 + ccc(i,j,k,1);
        Real delta_x = xf - ccc(i,j,k,0);
        Real delta_z = zf - ccc(i,j,k,2);

        const auto& slopes_eb_hi = amrex_lim_slopes_eb(i,j  ,k,comp, state, ccc, fcx, fcy, fcz, flag);

        Real vpls = state_pls - delta_y * slopes_eb_hi[1]
                              + delta_x * slopes_eb_hi[0]
                              + delta_z * slopes_eb_hi[2];

        vpls = amrex::max( amrex::min(vpls, cc_umax), cc_umin );

        delta_y = .5 - ccc(i,j-1,k,1);
        delta_x = xf - ccc(i,j-1,k,0);
        delta_z = zf - ccc(i,j-1,k,2);

        const auto& slopes_eb_lo = amrex_lim_slopes_eb(i,j-1,k,comp, state, ccc, fcx, fcy, fcz, flag);

        Real vmns = state_mns + delta_y * slopes_eb_lo[1]
                              + delta_x * slopes_eb_lo[0]
                              + delta_z * slopes_eb_lo[2];

        vmns = amrex::max( amrex::min(vmns, cc_umax), cc_umin );

        sy = upwind(vmns, vpls, ep_v_mac(i,j,k));
      }
      else {
        sy = my_huge;
      }

      fy(i,j,k,n) = ep_v_mac(i,j,k) * sy;

    });

  }

  /*********************************************************************************
   *                                                                               *
   *                                  z-direction                                  *
   *                                                                               *
   *********************************************************************************/
  if (check_extdir)
  {

    ParallelFor(zbx,ncomp,[fz,ep_w_mac,state,state_comp,ccc,fcx,fcy,fcz,flag,
    check_extdir_ilo, check_extdir_ihi, domain_ilo, domain_ihi, bct_ilo, bct_ihi,
    check_extdir_jlo, check_extdir_jhi, domain_jlo, domain_jhi, bct_jlo, bct_jhi,
    check_extdir_klo, check_extdir_khi, domain_klo, domain_khi, bct_klo, bct_khi,
    bc_types, dirichlet_bcs, dirichlet_size]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {

      Real sz;

      if( flag(i,j,k).isConnected(0,0,-1) ) {

        const int* bct_data = bc_types.data();
        const int  bct_size = bc_types.size();

        bool valid_bc_klo = aux::any_of(&bct_data[0], &bct_data[bct_size],
                      aux::is_equal<int>(bct_klo(i,j,domain_klo-1,0)));

        bool valid_bc_khi = aux::any_of(&bct_data[0], &bct_data[bct_size],
                      aux::is_equal<int>(bct_khi(i,j,domain_khi+1,0)));

        int comp = n+state_comp;

        if (k <= domain_klo && valid_bc_klo) {
          sz = state(i, j,domain_klo-1,comp);

        } else if (k >= domain_khi+1 && valid_bc_khi) {
          sz = state(i, j,domain_khi+1,comp);

        } else {

          bool extdir_ilo = check_extdir_ilo && (i-1 <= domain_ilo) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_ilo(domain_ilo-1,j,k,0)));

          bool extdir_ihi = check_extdir_ihi && (i   >= domain_ihi) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_ihi(domain_ihi+1,j,k,0)));

          bool extdir_jlo = check_extdir_jlo && (j-1 <= domain_jlo) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_jlo(i,domain_jlo-1,k,0)));

          bool extdir_jhi = check_extdir_jhi && (j   >= domain_jlo) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_jhi(i,domain_jhi+1,k,0)));

          bool extdir_klo = check_extdir_klo && (k-1 <= domain_klo) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_klo(i,j,domain_klo-1,0)));

          bool extdir_khi = check_extdir_khi && (k   >= domain_khi) &&
            aux::any_of(&dirichlet_bcs[0], &dirichlet_bcs[dirichlet_size],
                        aux::is_equal<int>(bct_khi(i,j,domain_khi+1,0)));

          const Real state_pls = state(i,j,k  ,state_comp+n);
          const Real state_mns = state(i,j,k-1,state_comp+n);

          const Real cc_umax = amrex::max(state_pls, state_mns);
          const Real cc_umin = amrex::min(state_pls, state_mns);

          const Real xf = fcz(i,j,k,0);
          const Real yf = fcz(i,j,k,1);

          Real delta_x = xf - ccc(i,j,k,0);
          Real delta_y = yf - ccc(i,j,k,1);
          Real delta_z = .5 + ccc(i,j,k,2);

          const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j,k  ,comp,
                                        state, ccc, fcx, fcy, fcz, flag,
                                        extdir_ilo, extdir_jlo, extdir_klo,
                                        extdir_ihi, extdir_jhi, extdir_khi,
                                        domain_ilo, domain_jlo, domain_klo,
                                        domain_ihi, domain_jhi, domain_khi);

          Real wpls = state_pls - delta_z * slopes_eb_hi[2]
                                + delta_x * slopes_eb_hi[0]
                                + delta_y * slopes_eb_hi[1];

          wpls = amrex::max( amrex::min(wpls, cc_umax), cc_umin );

          delta_x = xf - ccc(i,j,k-1,0);
          delta_y = yf - ccc(i,j,k-1,1);
          delta_z = .5 - ccc(i,j,k-1,2);

          const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i,j,k-1,comp,
                                        state, ccc, fcx, fcy, fcz, flag,
                                        extdir_ilo, extdir_jlo, extdir_klo,
                                        extdir_ihi, extdir_jhi, extdir_khi,
                                        domain_ilo, domain_jlo, domain_klo,
                                        domain_ihi, domain_jhi, domain_khi);

          Real wmns = state_mns + delta_z * slopes_eb_lo[2]
                                + delta_x * slopes_eb_lo[0]
                                + delta_y * slopes_eb_lo[1];

          wmns = amrex::max( amrex::min(wmns, cc_umax), cc_umin );

          sz = upwind(wmns, wpls, ep_w_mac(i,j,k));
        }
      }
      else {
        sz = my_huge;
      }

      fz(i,j,k,n) = ep_w_mac(i,j,k) * sz;

    });
  }
  else
  {
    ParallelFor(zbx,ncomp,[fz,ep_w_mac,state,state_comp,ccc,fcx,fcy,fcz,flag]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {

      Real sz;

      if( flag(i,j,k).isConnected(0,0,-1) ) {

        int comp = n+state_comp;

        const Real state_pls = state(i,j,k  ,state_comp+n);
        const Real state_mns = state(i,j,k-1,state_comp+n);

        const Real cc_umax = amrex::max(state_pls, state_mns);
        const Real cc_umin = amrex::min(state_pls, state_mns);

        const Real xf = fcz(i,j,k,0);
        const Real yf = fcz(i,j,k,1);

        Real delta_x = xf - ccc(i,j,k,0);
        Real delta_y = yf - ccc(i,j,k,1);
        Real delta_z = .5 + ccc(i,j,k,2);

        const auto& slopes_eb_hi = amrex_lim_slopes_eb(i,j,k  ,comp, state, ccc, fcx, fcy, fcz, flag);

        Real wpls = state_pls - delta_z * slopes_eb_hi[2]
                              + delta_x * slopes_eb_hi[0]
                              + delta_y * slopes_eb_hi[1];

        wpls = amrex::max( amrex::min(wpls, cc_umax), cc_umin );

        delta_x = xf - ccc(i,j,k-1,0);
        delta_y = yf - ccc(i,j,k-1,1);
        delta_z = .5 - ccc(i,j,k-1,2);

        const auto& slopes_eb_lo = amrex_lim_slopes_eb(i,j,k-1,comp, state, ccc, fcx, fcy, fcz, flag);

        Real wmns = state_mns + delta_z * slopes_eb_lo[2]
                              + delta_x * slopes_eb_lo[0]
                              + delta_y * slopes_eb_lo[1];

        wmns = amrex::max( amrex::min(wmns, cc_umax), cc_umin );

        sz = upwind(wmns, wpls, ep_w_mac(i,j,k));

      }
      else {
        sz = my_huge;
      }

      fz(i,j,k,n) = ep_w_mac(i,j,k) * sz;

    });

  }
}

}//end mol namespace
#endif
