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
      subroutine solve_u_g_star(&
         slo, shi, lo, hi, ulo, uhi, vlo, vhi, wlo, whi, &
         u_g, v_g, w_g, u_go, p_g, ro_g, rop_g, &
         rop_go, ep_g, tau_u_g, d_e, flux_ge, flux_gn, flux_gt ,mu_g,  &
         f_gds, a_m, b_m, drag_bm, flag, dt, dx, dy, dz)&

         bind(C, name="solve_u_g_star")

! Module procedures ..................................................//
      USE u_g_conv_dif, only: conv_dif_u_g
      USE source_u_g_module, only: source_u_g, source_u_g_bc
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

      integer(c_int)     , intent(in   ) :: slo(3),shi(3)
      integer(c_int)     , intent(in   ) ::  lo(3), hi(3)
      integer(c_int)     , intent(in   ) :: ulo(3),uhi(3)
      integer(c_int)     , intent(in   ) :: vlo(3),vhi(3)
      integer(c_int)     , intent(in   ) :: wlo(3),whi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: u_go&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: flux_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: tau_u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

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
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(  out) :: a_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
      real(c_real), intent(  out) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: d_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
!.....................................................................//

! Initialize a_m and b_m
      a_m(:,:,:,:) =  0.0d0
      a_m(:,:,:,0) = -1.0d0
      b_m(:,:,:)   =  0.0d0
      d_e(:,:,:)   =  0.0d0

! calculate the convection-diffusion terms
      call conv_dif_u_g (slo, shi, lo, hi, ulo, uhi, vlo, vhi, wlo, whi,&
         a_m, mu_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt, &
         flag, dt, dx, dy, dz)

! calculate the source terms for the gas phase u-momentum eqs
      call source_u_g(slo, shi, lo, hi, a_m, b_m, dt, p_g, ep_g, ro_g, rop_g, rop_go, &
         u_go, tau_u_g, dx, dy, dz)

! modifications for bc
      call source_u_g_bc (slo, shi, A_m, b_m, flag, dy, dz)

! add in point sources
      if(point_source) call point_source_u_g (slo, shi, b_m, flag, dx, dy, dz)

! calculate coefficients for the pressure correction equation
      call calc_d(slo, shi, lo, hi, d_e, "X", a_m, ep_g, f_gds, flag, dx, dy, dz)

! handle special case where center coefficient is zero
      call adjust_a_g ('U', slo, shi, lo, hi, a_m, b_m, rop_g, dx, dy, dz)

! add in source terms for DEM drag coupling.
      if(des_continuum_coupled) &
         call gas_drag_u(slo, shi, &
                         a_m, b_m, f_gds, drag_bm, flag, dx, dy, dz)

      call calc_resid_vel (slo, shi, lo, hi, &
         u_g, v_g, w_g, a_m, b_m, &
         num_resid(resid_u), den_resid(resid_u), &
         resid(resid_u), max_resid(resid_u), &
         i_resid(resid_u),j_resid(resid_u),k_resid(resid_u), 'U')

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
   subroutine solve_v_g_star(&
      slo, shi, lo, hi, ulo, uhi, vlo, vhi, wlo, whi, &
      u_g, v_g, w_g, v_go, p_g, ro_g, rop_g, &
      rop_go, ep_g, tau_v_g, d_n, flux_ge, flux_gn, flux_gt, mu_g,  &
      f_gds, a_m, b_m, drag_bm, flag, dt, dx, dy, dz)&
      bind(C, name="solve_v_g_star")


! Module procedures ...................................................//
      USE v_g_conv_dif, only: conv_dif_v_g
      USE source_v_g_module, only: source_v_g, source_v_g_bc
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

! Dummy arguments ....................................................//
      integer(c_int)     , intent(in   ) :: slo(3),shi(3)
      integer(c_int)     , intent(in   ) ::  lo(3), hi(3)
      integer(c_int)     , intent(in   ) :: ulo(3),uhi(3)
      integer(c_int)     , intent(in   ) :: vlo(3),vhi(3)
      integer(c_int)     , intent(in   ) :: wlo(3),whi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: flux_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: v_go&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: tau_v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

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
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(  out) :: d_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
      real(c_real), intent(  out) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
!.....................................................................//

! Initialize a_m and b_m
      a_m(:,:,:,:) =  0.0d0
      a_m(:,:,:,0) = -1.0d0
      b_m(:,:,:)   =  0.0d0

! calculate the convection-diffusion terms
      call conv_dif_v_g (slo, shi, lo, hi, a_m, mu_g, u_g, v_g, w_g, &
         flux_ge, flux_gn, flux_gt, flag, dt, dx, dy, dz)

! calculate the source terms for the gas phase u-momentum eqs
      call source_v_g(slo, shi, lo, hi, a_m, b_m, dt, p_g, ep_g, ro_g, rop_g, rop_go, &
         v_go, tau_v_g, dx, dy, dz)

! modifications for bc
      call source_v_g_bc(slo, shi, A_m, b_m, flag, dx, dz)

! add in point sources
      if(point_source) call point_source_v_g (slo, shi, b_m, flag, dx, dy, dz)

! calculate coefficients for the pressure correction equation
      call calc_d(slo, shi, lo, hi, d_n, "Y", a_m, ep_g, f_gds, flag, dx, dy, dz)

! handle special case where center coefficient is zero
      call adjust_a_g('V',slo, shi, lo, hi, a_m, b_m, rop_g, dx, dy, dz)

! add in source terms for DEM drag coupling.
      if(des_continuum_coupled) &
         call gas_drag_v(slo, shi, &
                         a_m, b_m, f_gds, drag_bm, flag, dx, dy, dz)

      call calc_resid_vel (slo, shi, lo, hi, &
         v_g, w_g, u_g, a_m, b_m, &
         num_resid(resid_v), den_resid(resid_v), &
         resid(resid_v), max_resid(resid_v), &
         i_resid(resid_v),j_resid(resid_v),k_resid(resid_v), 'V')

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
   subroutine solve_w_g_star(&
      slo, shi, lo, hi, ulo, uhi, vlo, vhi, wlo, whi, &
      u_g, v_g, w_g, w_go, p_g, ro_g, rop_g, &
      rop_go, ep_g, tau_w_g, d_t, flux_ge, flux_gn, flux_gt, mu_g,  &
      f_gds, a_m, b_m, drag_bm, flag, dt, dx, dy, dz)&
      bind(C, name="solve_w_g_star")

! Module procedures ..................................................//
      use w_g_conv_dif, only: conv_dif_w_g
      use source_w_g_module, only: source_w_g, source_w_g_bc
      use source_w_g_module, only: point_source_w_g
      use calc_d_mod, only: calc_d
      use adjust_a, only: adjust_a_g
      use gas_drag_module, only: gas_drag_w
      use residual, only: calc_resid_vel
      use ur_facs, only: under_relax

! Global data .........................................................//
! Flag for coupling fluid and dem via gas-solid drag
   use discretelement, only: des_continuum_coupled
! Flag for existence of point souces
   use ps, only: point_source
! Global data arrays for residuals
   use residual, only: resid_w
   use residual, only: resid, num_resid, den_resid
   use residual, only: max_resid, i_resid, j_resid, k_resid

      IMPLICIT NONE

! Dummy arguments ....................................................//
      integer(c_int)     , intent(in   ) :: slo(3),shi(3)
      integer(c_int)     , intent(in   ) ::  lo(3), hi(3)
      integer(c_int)     , intent(in   ) :: ulo(3),uhi(3)
      integer(c_int)     , intent(in   ) :: vlo(3),vhi(3)
      integer(c_int)     , intent(in   ) :: wlo(3),whi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: flux_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), intent(inout) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: w_go&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: tau_w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

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
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(  out) :: d_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(  out) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
      real(c_real), intent(  out) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
!.....................................................................//

! Initialize a_m and b_m
      a_m(:,:,:,:) =  0.0d0
      a_m(:,:,:,0) = -1.0d0
      b_m(:,:,:)   =  0.0d0
      d_t(:,:,:)   =  0.0d0

! calculate the convection-diffusion terms
      call conv_dif_w_g(slo, shi, lo, hi, a_m, mu_g, u_g, v_g, w_g, &
                        flux_ge, flux_gn, flux_gt, flag, dt, dx, dy, dz)

! calculate the source terms for the gas phase u-momentum eqs
      call source_w_g(slo, shi, lo, hi, a_m, b_m, dt, p_g, ep_g, ro_g, rop_g, rop_go, &
         w_go, tau_w_g, dx, dy, dz)

! modifications for bc
      call source_w_g_bc (slo, shi, A_m, b_m, flag, dx, dy)

! add in point sources
      if(point_source) call point_source_w_g (slo, shi, b_m, flag, dx, dy, dz)

! calculate coefficients for the pressure correction equation
      call calc_d(slo, shi, lo, hi, d_t, "Z", a_m, ep_g, f_gds, flag, dx, dy, dz)

! handle special case where center coefficient is zero
      call adjust_a_g('W',slo, shi, lo, hi, a_m, b_m, rop_g, dx, dy, dz)

! add in source terms for DEM drag coupling.
      if(des_continuum_coupled) &
         call gas_drag_w(slo, shi, &
                         a_m, b_m, f_gds, drag_bm, flag, dx, dy, dz)

      call calc_resid_vel (slo, shi, lo, hi, &
         w_g, u_g, v_g, a_m, b_m, &
         num_resid(resid_w), den_resid(resid_w), &
         resid(resid_w), max_resid(resid_w), &
         i_resid(resid_w),j_resid(resid_w),k_resid(resid_w),'W')

      call under_relax (slo, shi, w_g, a_m, b_m, 'W', flag, 5)

   END SUBROUTINE SOLVE_W_G_STAR

end module solve_vel_star_module
