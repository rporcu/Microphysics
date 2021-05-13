#ifndef MFIX_MOL_PREDICT_K_H_
#define MFIX_MOL_PREDICT_K_H_

#include <mfix_algorithm.H>
#include <AMReX_Slopes_K.H>

using namespace amrex;

namespace mol {

void predict_vels_on_faces (Box const& domain_box,
                            Box const& ubx,
                            Box const& vbx,
                            Box const& wbx,
                            Array4<Real> const& u,
                            Array4<Real> const& v,
                            Array4<Real> const& w,
                            Array4<Real const> const& vcc,
                            std::map<std::string, Gpu::DeviceVector<int>>& bc_types,
                            Array4<int const> const& bc_ilo,
                            Array4<int const> const& bc_ihi,
                            Array4<int const> const& bc_jlo,
                            Array4<int const> const& bc_jhi,
                            Array4<int const> const& bc_klo,
                            Array4<int const> const& bc_khi)
{
  constexpr int order = 2;
  constexpr Real small_vel = 1.e-8;

  const int* bct_Dirichlet = (bc_types["Dirichlet"]).data();
  const int  bct_size = (bc_types["Dirichlet"]).size();

  /*********************************************************************************
   *                                                                               *
   *                                  x-direction                                  *
   *                                                                               *
   *********************************************************************************/

  const int domain_ilo = domain_box.smallEnd(0);
  const int domain_ihi = domain_box.bigEnd(0);

  bool check_extdir_ilo = (domain_ilo >= ubx.smallEnd(0) && domain_ilo <= ubx.bigEnd(0));
  bool check_extdir_ihi = (domain_ihi >= ubx.smallEnd(0) && domain_ihi <= ubx.bigEnd(0));

  if (check_extdir_ilo || check_extdir_ihi)
  {
    amrex::ParallelFor(ubx, [vcc, small_vel, u, order,
      check_extdir_ilo, check_extdir_ihi, domain_ilo, domain_ihi,
      bc_ilo, bc_ihi, bct_Dirichlet, bct_size]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      bool extdir_ilo = check_extdir_ilo && (i-1 <= domain_ilo) &&
        (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                     aux::is_equal<int>(bc_ilo(domain_ilo-1,j,k,0))));

      bool extdir_ihi = check_extdir_ihi && (i   >= domain_ihi) &&
        (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                     aux::is_equal<int>(bc_ihi(domain_ihi+1,j,k,0))));

      const Real vcc_pls = vcc(i  ,j,k,0);
      const Real vcc_mns = vcc(i-1,j,k,0);

      const Real upls = vcc_pls - 0.5 * amrex_calc_xslope_extdir(
             i  ,j,k,0,order,vcc,extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);

      const Real umns = vcc_mns + 0.5 * amrex_calc_xslope_extdir(
             i-1,j,k,0,order,vcc,extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);

      Real umac(0);

      if (umns >= 0 || upls <= 0) {
        Real avg = .5 * ( upls + umns );

        if (avg >= small_vel) {
          umac = umns;
        }
        else if(avg <= -small_vel) {
          umac = upls;
        }
      }

      u(i,j,k) = umac;

    });

  }
  else
  {
    amrex::ParallelFor(ubx, [vcc, small_vel, u, order]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const Real vcc_pls = vcc(i  ,j,k,0);
      const Real vcc_mns = vcc(i-1,j,k,0);

      const Real upls = vcc_pls - 0.5 * amrex_calc_xslope(i  ,j,k,0,order,vcc);
      const Real umns = vcc_mns + 0.5 * amrex_calc_xslope(i-1,j,k,0,order,vcc);

      Real umac(0);

      if (umns >= 0 || upls <= 0) {
        Real avg = .5 * ( upls + umns );

        if (avg >= small_vel) {
          umac = umns;
        }
        else if(avg <= -small_vel) {
          umac = upls;
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

  const int domain_jlo = domain_box.smallEnd(1);
  const int domain_jhi = domain_box.bigEnd(1);

  bool check_extdir_jlo = (domain_jlo >= vbx.smallEnd(1) && domain_jlo <= vbx.bigEnd(1));
  bool check_extdir_jhi = (domain_jhi >= vbx.smallEnd(1) && domain_jhi <= vbx.bigEnd(1));


  if(check_extdir_jlo || check_extdir_jhi)
  {

    amrex::ParallelFor(vbx, [vcc, small_vel, v, order,
    check_extdir_jlo, check_extdir_jhi, domain_jlo, domain_jhi,
    bc_jlo, bc_jhi, bct_Dirichlet, bct_size]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      bool extdir_jlo = check_extdir_jlo && (j-1 <= domain_jlo) &&
        (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                     aux::is_equal<int>(bc_jlo(i,domain_jlo-1,k,0))));

      bool extdir_jhi = check_extdir_jhi && (j   >= domain_jhi) &&
        (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                     aux::is_equal<int>(bc_jhi(i,domain_jhi+1,k,0))));

      const Real vcc_pls = vcc(i,j  ,k,1);
      const Real vcc_mns = vcc(i,j-1,k,1);

      const Real vpls = vcc_pls - 0.5 * amrex_calc_yslope_extdir(
             i,j  ,k,1,order,vcc,extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);

      const Real vmns = vcc_mns + 0.5 * amrex_calc_yslope_extdir(
             i,j-1,k,1,order,vcc,extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);

      Real vmac(0);

      if (vmns >= 0 || vpls <= 0) {
        Real avg = .5 * (vpls + vmns);

        if (avg >= small_vel) {
          vmac = vmns;
        }
        else if (avg <= -small_vel) {
          vmac = vpls;
        }
      }

      v(i,j,k) = vmac;

    });

  }
  else
  {

    amrex::ParallelFor(vbx, [vcc, small_vel, v, order]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const Real vcc_pls = vcc(i,j  ,k,1);
      const Real vcc_mns = vcc(i,j-1,k,1);

      const Real vpls = vcc_pls - 0.5 * amrex_calc_yslope(i,j  ,k,1,order,vcc);
      const Real vmns = vcc_mns + 0.5 * amrex_calc_yslope(i,j-1,k,1,order,vcc);

      Real vmac(0);

      if (vmns >= 0 || vpls <= 0) {
        Real avg = .5 * (vpls + vmns);

        if (avg >= small_vel) {
          vmac = vmns;
        }
        else if (avg <= -small_vel) {
          vmac = vpls;
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

  const int domain_klo = domain_box.smallEnd(2);
  const int domain_khi = domain_box.bigEnd(2);

  bool check_extdir_klo = (domain_klo >= wbx.smallEnd(2) && domain_klo <= wbx.bigEnd(2));
  bool check_extdir_khi = (domain_khi >= wbx.smallEnd(2) && domain_khi <= wbx.bigEnd(2));

  if(check_extdir_klo || check_extdir_khi)
  {
    amrex::ParallelFor(wbx, [vcc,small_vel, w, order,
    check_extdir_klo, check_extdir_khi, domain_klo, domain_khi,
    bc_klo, bc_khi, bct_Dirichlet, bct_size]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      bool extdir_klo = check_extdir_klo && (k-1 <= domain_klo) &&
        (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                     aux::is_equal<int>(bc_klo(i,j,domain_klo-1,0))));

      bool extdir_khi = check_extdir_khi && (k   >= domain_khi) &&
        (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                     aux::is_equal<int>(bc_khi(i,j,domain_khi+1,0))));

      const Real vcc_pls = vcc(i,j,k  ,2);
      const Real vcc_mns = vcc(i,j,k-1,2);

      const Real wpls = vcc_pls - 0.5 * amrex_calc_zslope_extdir(
             i,j,k  ,2,order,vcc,extdir_klo, extdir_khi, domain_klo, domain_khi);

      const Real wmns = vcc_mns + 0.5 * amrex_calc_zslope_extdir(
             i,j,k-1,2,order,vcc,extdir_klo, extdir_khi, domain_klo, domain_khi);

      Real wmac(0);

      if (wmns >= 0 || wpls <= 0) {
        Real avg = .5 * (wpls + wmns);

        if (avg >= small_vel) {
          wmac = wmns;
        }
        else if (avg <= -small_vel) {
          wmac = wpls;
        }
      }

      w(i,j,k) = wmac;

    });
  }
  else
  {
    amrex::ParallelFor(wbx, [vcc,small_vel,w,order]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const Real vcc_pls = vcc(i,j,k  ,2);
      const Real vcc_mns = vcc(i,j,k-1,2);

      const Real wpls = vcc_pls - 0.5 * amrex_calc_zslope(i,j,k  ,2,order,vcc);
      const Real wmns = vcc_mns + 0.5 * amrex_calc_zslope(i,j,k-1,2,order,vcc);

      Real wmac(0);

      if (wmns >= 0 || wpls <= 0) {
        Real avg = .5 * (wpls + wmns);

        if (avg >= small_vel) {
          wmac = wmns;
        }
        else if (avg <= -small_vel) {
          wmac = wpls;
        }
      }

      w(i,j,k) = wmac;

    });

  }

}

}// end mol namespace
#endif
