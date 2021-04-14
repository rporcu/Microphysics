#ifndef MFIX_MOL_PREDICT_EB_K_H_
#define MFIX_MOL_PREDICT_EB_K_H_

#include <mfix_algorithm.H>
#include <AMReX_EB_slopes_K.H>

using namespace amrex;

namespace mol {

void predict_vels_on_faces_eb (Box const& domain_box,
                               Box const& ubx,
                               Box const& vbx,
                               Box const& wbx,
                               Array4<Real      > const& u,
                               Array4<Real      > const& v,
                               Array4<Real      > const& w,
                               Array4<Real const> const& vcc,
                               Array4<EBCellFlag const> const& flag,
                               Array4<Real const> const& fcx,
                               Array4<Real const> const& fcy,
                               Array4<Real const> const& fcz,
                               Array4<Real const> const& ccc,
                               std::map<std::string, Gpu::DeviceVector<int>>& bc_types,
                               Array4<int const> const& bc_ilo,
                               Array4<int const> const& bc_ihi,
                               Array4<int const> const& bc_jlo,
                               Array4<int const> const& bc_jhi,
                               Array4<int const> const& bc_klo,
                               Array4<int const> const& bc_khi)
{
  constexpr Real small_vel = 1.e-10;

  const int* bct_Dirichlet = (bc_types["Dirichlet"]).data();
  const int  bct_size = (bc_types["Dirichlet"]).size();

  const int domain_ilo = domain_box.smallEnd(0);
  const int domain_ihi = domain_box.bigEnd(0);
  const int domain_jlo = domain_box.smallEnd(1);
  const int domain_jhi = domain_box.bigEnd(1);
  const int domain_klo = domain_box.smallEnd(2);
  const int domain_khi = domain_box.bigEnd(2);

  const bool check_extdir_ilo = (domain_ilo >= ubx.smallEnd(0) && domain_ilo <= ubx.bigEnd(0));
  const bool check_extdir_ihi = (domain_ihi >= ubx.smallEnd(0) && domain_ihi <= ubx.bigEnd(0));
  const bool check_extdir_jlo = (domain_jlo >= vbx.smallEnd(1) && domain_jlo <= vbx.bigEnd(1));
  const bool check_extdir_jhi = (domain_jhi >= vbx.smallEnd(1) && domain_jhi <= vbx.bigEnd(1));
  const bool check_extdir_klo = (domain_klo >= wbx.smallEnd(2) && domain_klo <= wbx.bigEnd(2));
  const bool check_extdir_khi = (domain_khi >= wbx.smallEnd(2) && domain_khi <= wbx.bigEnd(2));

  const bool check_extdir = check_extdir_ilo || check_extdir_ihi ||
                            check_extdir_jlo || check_extdir_jhi ||
                            check_extdir_klo || check_extdir_khi;

  /*********************************************************************************
   *                                                                               *
   *                                  x-direction                                  *
   *                                                                               *
   *********************************************************************************/

  if (check_extdir)
  {

    // This FAB has cut cells -- we predict from cell centroids to face
    // centroids
    ParallelFor(ubx, [u,vcc,ccc,fcx,fcy,fcz,flag,small_vel,bct_Dirichlet, bct_size,
    check_extdir_ilo, check_extdir_ihi, domain_ilo, domain_ihi, bc_ilo, bc_ihi,
    check_extdir_jlo, check_extdir_jhi, domain_jlo, domain_jhi, bc_jlo, bc_jhi,
    check_extdir_klo, check_extdir_khi, domain_klo, domain_khi, bc_klo, bc_khi]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      Real umac(0);

      // X-faces
      if (flag(i,j,k).isConnected(-1,0,0))
      {
        bool extdir_ilo = check_extdir_ilo && (i-1 <= domain_ilo) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_ilo(domain_ilo-1,j,k,0))));

        bool extdir_ihi = check_extdir_ihi && (i   >= domain_ihi) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_ihi(domain_ihi+1,j,k,0))));

        bool extdir_jlo = check_extdir_jlo && (j-1 <= domain_jlo) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_jlo(i,domain_jlo-1,k,0))));

        bool extdir_jhi = check_extdir_jhi && (j   >= domain_jhi) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_jhi(i,domain_jhi+1,k,0))));

        bool extdir_klo = check_extdir_klo && (k-1 <= domain_klo) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_klo(i,j,domain_klo-1,0))));

        bool extdir_khi = check_extdir_khi && (k   >= domain_khi) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_khi(i,j,domain_khi+1,0))));

        // local (y,z) of centroid of x-face we are extrapolating to
        Real yf = fcx(i,j,k,0);
        Real zf = fcx(i,j,k,1);

        Real delta_x = 0.5 + ccc(i,j,k,0);
        Real delta_y = yf  - ccc(i,j,k,1);
        Real delta_z = zf  - ccc(i,j,k,2);

        const Real vcc_mns = vcc(i-1,j,k,0);
        const Real vcc_pls = vcc(i  ,j,k,0);

        Real cc_umax = amrex::max(vcc_pls, vcc_mns);
        Real cc_umin = amrex::min(vcc_pls, vcc_mns);

        // Compute slopes of component "0" of vcc
        const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i  ,j,k,0,
                                      vcc, ccc, fcx, fcy, fcz, flag,
                                      extdir_ilo, extdir_jlo, extdir_klo,
                                      extdir_ihi, extdir_jhi, extdir_khi,
                                      domain_ilo, domain_jlo, domain_klo,
                                      domain_ihi, domain_jhi, domain_khi);

        Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                            + delta_y * slopes_eb_hi[1]
                            + delta_z * slopes_eb_hi[2];

        upls = amrex::min(upls, cc_umax);
        upls = amrex::max(upls, cc_umin);

        delta_x = .5 - ccc(i-1,j,k,0);
        delta_y = yf - ccc(i-1,j,k,1);
        delta_z = zf - ccc(i-1,j,k,2);

        const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i-1,j,k,0,
                                      vcc, ccc, fcx, fcy, fcz, flag,
                                      extdir_ilo, extdir_jlo, extdir_klo,
                                      extdir_ihi, extdir_jhi, extdir_khi,
                                      domain_ilo, domain_jlo, domain_klo,
                                      domain_ihi, domain_jhi, domain_khi);

        Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                            + delta_y * slopes_eb_lo[1]
                            + delta_z * slopes_eb_lo[2];

        umns = amrex::min(umns, cc_umax);
        umns = amrex::max(umns, cc_umin);


        if (umns >= 0 || upls <= 0) {
          Real avg = .5 * (upls + umns);

          if (avg >= small_vel) {
            umac = umns;
          }
          else if (avg <= -small_vel) {
            umac = upls;
          }
        }
      }

      u(i,j,k) = umac;

    });

  }
  else
  {

    // This FAB has cut cells -- we predict from cell centroids to face
    // centroids
    ParallelFor(ubx, [u,vcc,ccc,fcx,fcy,fcz,flag,small_vel]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      Real umac(0);

      // X-faces
      if (flag(i,j,k).isConnected(-1,0,0))
      {
        // local (y,z) of centroid of x-face we are extrapolating to
        Real yf = fcx(i,j,k,0);
        Real zf = fcx(i,j,k,1);

        Real delta_x = 0.5 + ccc(i,j,k,0);
        Real delta_y = yf  - ccc(i,j,k,1);
        Real delta_z = zf  - ccc(i,j,k,2);

        const Real vcc_mns = vcc(i-1,j,k,0);
        const Real vcc_pls = vcc(i  ,j,k,0);

        Real cc_umax = amrex::max(vcc_pls, vcc_mns);
        Real cc_umin = amrex::min(vcc_pls, vcc_mns);

        // Compute slopes of component "0" of vcc
        const auto& slopes_eb_hi = amrex_lim_slopes_eb(i  ,j,k,0,vcc,ccc,fcx,fcy,fcz,flag);

        Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                            + delta_y * slopes_eb_hi[1]
                            + delta_z * slopes_eb_hi[2];

        upls = amrex::min(upls, cc_umax);
        upls = amrex::max(upls, cc_umin);

        delta_x = .5 - ccc(i-1,j,k,0);
        delta_y = yf - ccc(i-1,j,k,1);
        delta_z = zf - ccc(i-1,j,k,2);

        const auto& slopes_eb_lo = amrex_lim_slopes_eb(i-1,j,k,0,vcc,ccc,fcx,fcy,fcz,flag);

        Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                            + delta_y * slopes_eb_lo[1]
                            + delta_z * slopes_eb_lo[2];

        umns = amrex::min(umns, cc_umax);
        umns = amrex::max(umns, cc_umin);


        if (umns >= 0 || upls <= 0) {
          Real avg = .5 * (upls + umns);

          if (avg >= small_vel) {
            umac = umns;
          }
          else if (avg <= -small_vel) {
            umac = upls;
          }
        }
      }

      u(i,j,k) = umac;

    });
  }



  /*********************************************************************************
   *                                                                               *
   *                                  y-direction                                  *
   *                                                                               *
   *********************************************************************************/

  if (check_extdir)
  {

    // This FAB has cut cells -- we predict from cell centroids to face
    // centroids
    ParallelFor(vbx, [v,vcc,ccc,fcx,fcy,fcz,flag,small_vel,bct_Dirichlet, bct_size,
    check_extdir_ilo, check_extdir_ihi, domain_ilo, domain_ihi, bc_ilo, bc_ihi,
    check_extdir_jlo, check_extdir_jhi, domain_jlo, domain_jhi, bc_jlo, bc_jhi,
    check_extdir_klo, check_extdir_khi, domain_klo, domain_khi, bc_klo, bc_khi]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      Real vmac(0);

      if (flag(i,j,k).isConnected(0,-1,0))
      {
        bool extdir_ilo = check_extdir_ilo && (i-1 <= domain_ilo) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_ilo(domain_ilo-1,j,k,0))));

        bool extdir_ihi = check_extdir_ihi && (i   >= domain_ihi) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_ihi(domain_ihi+1,j,k,0))));

        bool extdir_jlo = check_extdir_jlo && (j-1 <= domain_jlo) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_jlo(i,domain_jlo-1,k,0))));

        bool extdir_jhi = check_extdir_jhi && (j   >= domain_jhi) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_jhi(i,domain_jhi+1,k,0))));

        bool extdir_klo = check_extdir_klo && (k-1 <= domain_klo) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_klo(i,j,domain_klo-1,0))));

        bool extdir_khi = check_extdir_khi && (k   >= domain_khi) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_khi(i,j,domain_khi+1,0))));

        Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
        Real zf = fcy(i,j,k,1);

        Real delta_x = xf - ccc(i,j,k,0);
        Real delta_y = .5 + ccc(i,j,k,1);
        Real delta_z = zf - ccc(i,j,k,2);

        const Real vcc_pls = vcc(i,j,k,1);
        const Real vcc_mns = vcc(i,j-1,k,1);

        Real cc_vmax = amrex::max(vcc_pls, vcc_mns);
        Real cc_vmin = amrex::min(vcc_pls, vcc_mns);

        const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j  ,k,1,
                                      vcc, ccc, fcx, fcy, fcz, flag,
                                      extdir_ilo, extdir_jlo, extdir_klo,
                                      extdir_ihi, extdir_jhi, extdir_khi,
                                      domain_ilo, domain_jlo, domain_klo,
                                      domain_ihi, domain_jhi, domain_khi);

        Real vpls = vcc_pls - delta_y * slopes_eb_hi[1]
                            + delta_x * slopes_eb_hi[0]
                            + delta_z * slopes_eb_hi[2];

        vpls = amrex::min(vpls, cc_vmax);
        vpls = amrex::max(vpls, cc_vmin);

        delta_x = xf - ccc(i,j-1,k,0);
        delta_y = .5 - ccc(i,j-1,k,1);
        delta_z = zf - ccc(i,j-1,k,2);

        const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i,j-1,k,1,
                                      vcc, ccc, fcx, fcy, fcz, flag,
                                      extdir_ilo, extdir_jlo, extdir_klo,
                                      extdir_ihi, extdir_jhi, extdir_khi,
                                      domain_ilo, domain_jlo, domain_klo,
                                      domain_ihi, domain_jhi, domain_khi);

        Real vmns = vcc_mns + delta_y * slopes_eb_lo[1]
                            + delta_x * slopes_eb_lo[0]
                            + delta_z * slopes_eb_lo[2];

        vmns = amrex::min(vmns, cc_vmax);
        vmns = amrex::max(vmns, cc_vmin);

        if (vmns >= 0 || vpls <= 0) {
          Real avg = .5 * (vpls + vmns);

          if (avg >= small_vel) {
            vmac = vmns;
          }
          else if (avg <= -small_vel) {
            vmac = vpls;
          }
        }
      }

      v(i,j,k) = vmac;

    });

  }
  else
  {

    // This FAB has cut cells -- we predict from cell centroids to face
    // centroids
    ParallelFor(vbx, [v,vcc,ccc,fcx,fcy,fcz,flag,small_vel]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      Real vmac(0);

      if (flag(i,j,k).isConnected(0,-1,0))
      {

        Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
        Real zf = fcy(i,j,k,1);

        Real delta_x = xf - ccc(i,j,k,0);
        Real delta_y = .5 + ccc(i,j,k,1);
        Real delta_z = zf - ccc(i,j,k,2);

        const Real vcc_pls = vcc(i,j,k,1);
        const Real vcc_mns = vcc(i,j-1,k,1);

        Real cc_vmax = amrex::max(vcc_pls, vcc_mns);
        Real cc_vmin = amrex::min(vcc_pls, vcc_mns);

        const auto& slopes_eb_hi = amrex_lim_slopes_eb(i,j  ,k,1,vcc,ccc,fcx,fcy,fcz,flag);

        Real vpls = vcc_pls - delta_y * slopes_eb_hi[1]
                            + delta_x * slopes_eb_hi[0]
                            + delta_z * slopes_eb_hi[2];

        vpls = amrex::min(vpls, cc_vmax);
        vpls = amrex::max(vpls, cc_vmin);

        delta_x = xf - ccc(i,j-1,k,0);
        delta_y = .5 - ccc(i,j-1,k,1);
        delta_z = zf - ccc(i,j-1,k,2);

        const auto& slopes_eb_lo = amrex_lim_slopes_eb(i,j-1,k,1,vcc,ccc,fcx,fcy,fcz,flag);

        Real vmns = vcc_mns + delta_y * slopes_eb_lo[1]
                            + delta_x * slopes_eb_lo[0]
                            + delta_z * slopes_eb_lo[2];

        vmns = amrex::min(vmns, cc_vmax);
        vmns = amrex::max(vmns, cc_vmin);

        if (vmns >= 0 || vpls <= 0) {
          Real avg = .5 * (vpls + vmns);

          if (avg >= small_vel) {
            vmac = vmns;
          }
          else if (avg <= -small_vel) {
            vmac = vpls;
          }
        }
      }

      v(i,j,k) = vmac;

    });
  }

  /*********************************************************************************
   *                                                                               *
   *                                  z-direction                                  *
   *                                                                               *
   *********************************************************************************/

  if (check_extdir)
  {

    // This FAB has cut cells -- we predict from cell centroids to face
    // centroids
    ParallelFor(wbx, [w,vcc,ccc,fcx,fcy,fcz,flag,small_vel,bct_Dirichlet, bct_size,
    check_extdir_ilo, check_extdir_ihi, domain_ilo, domain_ihi, bc_ilo, bc_ihi,
    check_extdir_jlo, check_extdir_jhi, domain_jlo, domain_jhi, bc_jlo, bc_jhi,
    check_extdir_klo, check_extdir_khi, domain_klo, domain_khi, bc_klo, bc_khi]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      Real wmac(0);

      if (flag(i,j,k).isConnected(0,0,-1))
      {
        bool extdir_ilo = check_extdir_ilo && (i-1 <= domain_ilo) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_ilo(domain_ilo-1,j,k,0))));

        bool extdir_ihi = check_extdir_ihi && (i   >= domain_ihi) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_ihi(domain_ihi+1,j,k,0))));

        bool extdir_jlo = check_extdir_jlo && (j-1 <= domain_jlo) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_jlo(i,domain_jlo-1,k,0))));

        bool extdir_jhi = check_extdir_jhi && (j   >= domain_jhi) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_jhi(i,domain_jhi+1,k,0))));

        bool extdir_klo = check_extdir_klo && (k-1 <= domain_klo) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_klo(i,j,domain_klo-1,0))));

        bool extdir_khi = check_extdir_khi && (k   >= domain_khi) &&
          (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                       aux::is_equal<int>(bc_khi(i,j,domain_khi+1,0))));

        Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
        Real yf = fcz(i,j,k,1);

        Real delta_x = xf - ccc(i,j,k,0);
        Real delta_y = yf - ccc(i,j,k,1);
        Real delta_z = .5 + ccc(i,j,k,2);

        const Real vcc_pls = vcc(i,j,k,2);
        const Real vcc_mns = vcc(i,j,k-1,2);

        Real cc_wmax = amrex::max(vcc_pls, vcc_mns);
        Real cc_wmin = amrex::min(vcc_pls, vcc_mns);

        const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j,k  ,2,
                                      vcc, ccc, fcx, fcy, fcz, flag,
                                      extdir_ilo, extdir_jlo, extdir_klo,
                                      extdir_ihi, extdir_jhi, extdir_khi,
                                      domain_ilo, domain_jlo, domain_klo,
                                      domain_ihi, domain_jhi, domain_khi);

        Real wpls = vcc_pls - delta_z * slopes_eb_hi[2]
                            + delta_x * slopes_eb_hi[0]
                            + delta_y * slopes_eb_hi[1];

        wpls = amrex::min(wpls, cc_wmax);
        wpls = amrex::max(wpls, cc_wmin);

        delta_x = xf - ccc(i,j,k-1,0);
        delta_y = yf - ccc(i,j,k-1,1);
        delta_z = 0.5 - ccc(i,j,k-1,2);

        const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i,j,k-1,2,
                                      vcc, ccc, fcx, fcy, fcz, flag,
                                      extdir_ilo, extdir_jlo, extdir_klo,
                                      extdir_ihi, extdir_jhi, extdir_khi,
                                      domain_ilo, domain_jlo, domain_klo,
                                      domain_ihi, domain_jhi, domain_khi);

        Real wmns = vcc_mns + delta_z * slopes_eb_lo[2]
                            + delta_x * slopes_eb_lo[0]
                            + delta_y * slopes_eb_lo[1];

        wmns = amrex::min(wmns, cc_wmax);
        wmns = amrex::max(wmns, cc_wmin);

        if ( wmns >= 0 || wpls <= 0) {
          Real avg = .5 * (wpls + wmns);

          if (avg >= small_vel) {
            wmac = wmns;
          }
          else if (avg <= -small_vel) {
            wmac = wpls;
          }
        }
      }

      w(i,j,k) = wmac;

    });




  }
  else
  {
    // This FAB has cut cells -- we predict from cell centroids to face
    // centroids
    ParallelFor(wbx, [w,vcc,ccc,fcx,fcy,fcz,flag,small_vel]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      Real wmac(0);

      if (flag(i,j,k).isConnected(0,0,-1))
      {

        Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
        Real yf = fcz(i,j,k,1);

        Real delta_x = xf - ccc(i,j,k,0);
        Real delta_y = yf - ccc(i,j,k,1);
        Real delta_z = .5 + ccc(i,j,k,2);

        const Real vcc_pls = vcc(i,j,k,2);
        const Real vcc_mns = vcc(i,j,k-1,2);

        Real cc_wmax = amrex::max(vcc_pls, vcc_mns);
        Real cc_wmin = amrex::min(vcc_pls, vcc_mns);

        const auto& slopes_eb_hi = amrex_lim_slopes_eb(i,j,  k,2,vcc,ccc,fcx,fcy,fcz,flag);

        Real wpls = vcc_pls - delta_z * slopes_eb_hi[2]
                            + delta_x * slopes_eb_hi[0]
                            + delta_y * slopes_eb_hi[1];

        wpls = amrex::min(wpls, cc_wmax);
        wpls = amrex::max(wpls, cc_wmin);

        delta_x = xf - ccc(i,j,k-1,0);
        delta_y = yf - ccc(i,j,k-1,1);
        delta_z = 0.5 - ccc(i,j,k-1,2);

        const auto& slopes_eb_lo = amrex_lim_slopes_eb(i,j,k-1,2,vcc,ccc,fcx,fcy,fcz,flag);

        Real wmns = vcc_mns + delta_z * slopes_eb_lo[2]
                            + delta_x * slopes_eb_lo[0]
                            + delta_y * slopes_eb_lo[1];

        wmns = amrex::min(wmns, cc_wmax);
        wmns = amrex::max(wmns, cc_wmin);

        if ( wmns >= 0 || wpls <= 0) {
          Real avg = .5 * (wpls + wmns);

          if (avg >= small_vel) {
            wmac = wmns;
          }
          else if (avg <= -small_vel) {
            wmac = wpls;
          }
        }
      }

      w(i,j,k) = wmac;

    });
    }
}


}// end mol namespace
#endif
