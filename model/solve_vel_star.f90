module solve_vel_star_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   private

   public :: solve_u_g_star
   public :: solve_v_g_star
   public :: solve_w_g_star

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOLVE_VEL_STAR                                          !
!  Author: M. Syamlal                                 Date: 25-APR-96  !
!                                                                      !
!  Purpose: Solve starred velocity components                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine solve_u_g_star(slo, shi, lo, hi, u_g, v_g, w_g, u_go, p_g, ro_g, rop_g, &
         rop_go, ep_g, tau_u_g, d_e, flux_ge, flux_gn, flux_gt ,mu_g,  &
         f_gds, a_m, b_m, drag_bm, flag, dt, dx, dy, dz)&

         bind(C, name="solve_u_g_star")

! Module procedures ..................................................//
      USE matrix, only: init_ab_m
      USE u_g_conv_dif, only: conv_dif_u_g
      USE source_u_g_module, only: source_u_g
      USE source_u_g_module, only: point_source_u_g
      USE calc_d_mod, only: calc_d
      USE adjust_a, only: adjust_a_g
      use gas_drag_module, only: gas_drag_u
      use residual, only: calc_resid_vel
      use ur_facs, only: under_relax

! Global data .........................................................//
! Fluid array bounds
! Flag for coupling fluid and dem via gas-solid drag
   use discretelement, only: des_continuum_coupled
! Flag for existence of point souces
   use ps, only: point_source
! Global data arrays for residuals
   use residual, only: resid_u
   use residual, only: resid, num_resid, den_resid
   use residual, only: max_resid, i_resid, j_resid, k_resid

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      real(c_real), intent(in   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: tau_u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: d_e&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: a_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: dt, dx, dy, dz
!.....................................................................//

! initialize matrix and vector
      call init_ab_m (a_m, b_m)

! calculate the convection-diffusion terms
      call conv_dif_u_g (slo, shi, lo, hi, a_m, mu_g, u_g, v_g, w_g, &
         flux_ge, flux_gn, flux_gt, flag, dt, dx, dy, dz)

! calculate the source terms for the gas phase u-momentum eqs
      call source_u_g(slo, shi, lo, hi, a_m, b_m, dt, p_g, ep_g, ro_g, rop_g, rop_go, &
         u_g, u_go, tau_u_g, flag, dx, dy, dz)

! add in point sources
      if(point_source) call point_source_u_g (slo, shi, lo, hi, a_m, b_m, flag, dx, dy, dz)

! calculate coefficients for the pressure correction equation
      call calc_d(slo, shi, d_e, "X", a_m, ep_g, f_gds, flag, dx, dy, dz)

! handle special case where center coefficient is zero
      call adjust_a_g ('U', slo, shi, lo, hi, a_m, b_m, rop_g, dx, dy, dz)

! add in source terms for DEM drag coupling.
      if(des_continuum_coupled) &
         call gas_drag_u(a_m, b_m, f_gds, drag_bm, flag, dx, dy, dz)

      call calc_resid_vel (slo, shi, lo, hi, &
         u_g, v_g, w_g, a_m, b_m, &
         num_resid(resid_u), den_resid(resid_u), &
         resid(resid_u), max_resid(resid_u), &
         i_resid(resid_u),j_resid(resid_u),k_resid(resid_u), flag)

      call under_relax (slo, shi, u_g, a_m, b_m, 'U', flag, 3)

      return
   end subroutine solve_u_g_star


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOLVE_V_G_STAR                                          !
!  Author: M. Syamlal                                 Date: 25-APR-96  !
!                                                                      !
!  Purpose: Solve starred velocity components                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine solve_v_g_star(slo, shi, lo, hi, u_g, v_g, w_g, v_go, p_g, ro_g, rop_g, &
      rop_go, ep_g, tau_v_g, d_n, flux_ge, flux_gn, flux_gt, mu_g,  &
      f_gds, a_m, b_m, drag_bm, flag, dt, dx, dy, dz)&
      bind(C, name="solve_v_g_star")


! Module procedures ...................................................//
      USE matrix, only: init_ab_m
      USE v_g_conv_dif, only: conv_dif_v_g
      USE source_v_g_module, only: source_v_g
      USE source_v_g_module, only: point_source_v_g
      USE calc_d_mod, only: calc_d
      USE adjust_a, only: adjust_a_g
      use gas_drag_module, only: gas_drag_v
      use residual, only: calc_resid_vel
      use ur_facs, only: under_relax

! Global data .........................................................//
! Fluid array bounds
! Flag for coupling fluid and dem via gas-solid drag
   use discretelement, only: des_continuum_coupled
! Flag for existence of point souces
   use ps, only: point_source
! Global data arrays for residuals
   use residual, only: resid_v
   use residual, only: resid, num_resid, den_resid
   use residual, only: max_resid, i_resid, j_resid, k_resid

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

! Dummy arguments ....................................................//
      real(c_real), intent(inout) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: v_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: tau_v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: d_n&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
      real(c_real), intent(inout) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      INTEGER(C_INT), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), intent(in   ) :: dt, dx, dy, dz
!.....................................................................//

! initialize matrix and vector
      call init_ab_m (a_m, b_m)

! calculate the convection-diffusion terms
      call conv_dif_v_g (slo, shi, lo, hi, a_m, mu_g, u_g, v_g, w_g, &
         flux_ge, flux_gn, flux_gt, flag, dt, dx, dy, dz)

! calculate the source terms for the gas phase u-momentum eqs
      call source_v_g(slo, shi, lo, hi, a_m, b_m, dt, p_g, ep_g, ro_g, rop_g, rop_go, &
         v_g, v_go, tau_v_g, flag, dx, dy, dz)

! add in point sources
      if(point_source) call point_source_v_g (slo, shi, lo, hi, a_m, b_m, flag, dx, dy, dz)

! calculate coefficients for the pressure correction equation
      call calc_d(slo, shi, d_n, "Y", a_m, ep_g, f_gds, flag, dx, dy, dz)

! handle special case where center coefficient is zero
      call adjust_a_g('V',slo, shi, lo, hi, a_m, b_m, rop_g, dx, dy, dz)

! add in source terms for DEM drag coupling.
      if(des_continuum_coupled) &
         call gas_drag_v(a_m, b_m, f_gds, drag_bm, flag, dx, dy, dz)

      call calc_resid_vel (slo, shi, lo, hi, &
         v_g, w_g, u_g, a_m, b_m, &
         num_resid(resid_v), den_resid(resid_v), &
         resid(resid_v), max_resid(resid_v), &
         i_resid(resid_v),j_resid(resid_v),k_resid(resid_v), flag)

      call under_relax (slo, shi, v_g, a_m, b_m, 'V', flag, 4)

   END SUBROUTINE SOLVE_V_G_STAR


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOLVE_W_G_STAR                                          !
!  Author: M. Syamlal                                 Date: 25-APR-96  !
!                                                                      !
!  Purpose: Solve starred velocity components                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine solve_w_g_star(slo, shi, lo, hi, u_g, v_g, w_g, w_go, p_g, ro_g, rop_g, &
      rop_go, ep_g, tau_w_g, d_t, flux_ge, flux_gn, flux_gt, mu_g,  &
      f_gds, a_m, b_m, drag_bm, flag, dt, dx, dy, dz)&
      bind(C, name="solve_w_g_star")

! Module procedures ..................................................//
      USE matrix, only: init_ab_m
      USE w_g_conv_dif, only: conv_dif_w_g
      USE source_w_g_module, only: source_w_g
      USE source_w_g_module, only: point_source_w_g
      USE calc_d_mod, only: calc_d
      USE adjust_a, only: adjust_a_g
      use gas_drag_module, only: gas_drag_w
      use residual, only: calc_resid_vel
      use ur_facs, only: under_relax

! Global data .........................................................//
! Fluid array bounds
! Flag for coupling fluid and dem via gas-solid drag
   use discretelement, only: des_continuum_coupled
! Flag for existence of point souces
   use ps, only: point_source
! Global data arrays for residuals
   use residual, only: resid_w
   use residual, only: resid, num_resid, den_resid
   use residual, only: max_resid, i_resid, j_resid, k_resid

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

! Dummy arguments ....................................................//
      real(c_real), intent(inout) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: tau_w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: d_t&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
      real(c_real), intent(inout) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), intent(in   ) :: dt, dx, dy, dz
!.....................................................................//

! initialize matrix and vector
      call init_ab_m(a_m, b_m)

! calculate the convection-diffusion terms
      call conv_dif_w_g(slo, shi, lo, hi, a_m, mu_g, u_g, v_g, w_g, &
                        flux_ge, flux_gn, flux_gt, flag, dt, dx, dy, dz)

! calculate the source terms for the gas phase u-momentum eqs
      call source_w_g(slo, shi, lo, hi, a_m, b_m, dt, p_g, ep_g, ro_g, rop_g, rop_go, &
         w_g, w_go, tau_w_g, flag, dx, dy, dz)

! add in point sources
      if(point_source) call point_source_w_g (slo, shi, lo, hi, a_m, b_m, flag, dx, dy, dz)

! calculate coefficients for the pressure correction equation
      call calc_d(slo, shi, d_t, "Z", a_m, ep_g, f_gds, flag, dx, dy, dz)

! handle special case where center coefficient is zero
      call adjust_a_g('W',slo, shi, lo, hi, a_m, b_m, rop_g, dx, dy, dz)

! add in source terms for DEM drag coupling.
      if(des_continuum_coupled) &
         call gas_drag_w(a_m, b_m, f_gds, drag_bm, flag, dx, dy, dz)

      call calc_resid_vel (slo, shi, lo, hi, &
         w_g, u_g, v_g, a_m, b_m, &
         num_resid(resid_w), den_resid(resid_w), &
         resid(resid_w), max_resid(resid_w), &
         i_resid(resid_w),j_resid(resid_w),k_resid(resid_w),flag)

      call under_relax (slo, shi, w_g, a_m, b_m, 'W', flag, 5)

   END SUBROUTINE SOLVE_W_G_STAR

end module solve_vel_star_module
