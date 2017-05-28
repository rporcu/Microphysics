module solve_vel_star_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   implicit none

   private

   public :: solve_u_g_star
   public :: solve_v_g_star
   public :: solve_w_g_star

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: solve_u_g_star                                          !
!  Author: M. Syamlal                                 Date: 25-APR-96  !
!                                                                      !
!  Purpose: Solve starred velocity components                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine solve_u_g_star(&
         slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, lo, hi, &
         u_g, v_g, w_g, u_go, p_g, ro_g, rop_g, &
         rop_go, ep_g, tau_u_g, d_e, fluxX, fluxY, fluxZ ,mu_g,  &
         f_gds, A_m, b_m, drag_bm, &
         bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
         bc_klo_type, bc_khi_type, domlo, domhi, dt, dx, dy, dz, resid) &
         bind(C, name="solve_u_g_star")

      use u_g_conv_dif, only: conv_dif_u_g
      use source_u_g_module, only: source_u_g, source_u_g_bc
      use source_u_g_module, only: point_source_u_g
      use calc_d_mod, only: calc_d_e
      use adjust_a, only: adjust_a_g
      use gas_drag_module, only: gas_drag_u
      use residual, only: calc_resid_vel
      use ur_facs, only: under_relax
      use discretelement, only: des_continuum_coupled

      ! Flag for existence of point souces
      use ps, only: point_source

      ! Global data arrays for residuals
      use residual, only: resid_u

      integer(c_int)     , intent(in   ) :: slo(3),shi(3)
      integer(c_int)     , intent(in   ) :: ulo(3),uhi(3)
      integer(c_int)     , intent(in   ) :: vlo(3),vhi(3)
      integer(c_int)     , intent(in   ) :: wlo(3),whi(3)
      integer(c_int)     , intent(in   ) :: alo(3),ahi(3)
      integer(c_int)     , intent(in   ) :: domlo(3),domhi(3)
      integer(c_int)     , intent(in   ) ::  lo(3), hi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: u_go&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: fluxX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: tau_u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: fluxY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: fluxZ&
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

      real(c_real), intent(  out) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(  out) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real), intent(  out) :: d_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(  out) :: resid(8,2)

      integer(c_int), intent(in   ) :: bc_ilo_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jlo_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_klo_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)
      integer(c_int), intent(in   ) :: bc_khi_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)

      real(c_real) :: vol
      vol = dx*dy*dz

      ! Initialize A_m, b_m -- but only on the current tile!
      A_m(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) =  0.0d0
      A_m(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),0) = -1.0d0
      b_m(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))   =  0.0d0

      ! calculate the convection-diffusion terms
      call conv_dif_u_g (lo, hi, slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
                         A_m, mu_g, fluxX, fluxY, fluxZ, dx, dy, dz)

      ! calculate the source terms for the gas phase u-momentum eqs
      call source_u_g(lo, hi, slo, shi, ulo, uhi, alo, ahi, A_m, b_m, dt, &
         p_g, ep_g, ro_g, rop_go, u_go, tau_u_g, dx, dy, dz, domlo, domhi)

      ! modifications for bc
      call source_u_g_bc (lo, hi, slo, shi, alo, ahi, A_m, b_m, &
                          bc_ilo_type, bc_ihi_type, &
                          bc_jlo_type, bc_jhi_type, &
                          bc_klo_type, bc_khi_type, &
                          domlo, domhi, dy, dz)

      ! Add in point sources
      if(point_source) call point_source_u_g (lo, hi, alo, ahi, b_m, vol)

      ! Calculate coefficients for the pressure correction equation
      call calc_d_e(lo, hi, slo, shi, ulo, uhi, alo, ahi, d_e, A_m, &
                    ep_g, f_gds, dx, dy, dz, domlo, domhi)

      ! Handle special case where center coefficient is zero
      call adjust_a_g ('U', slo, shi, alo, ahi, lo, hi, A_m, b_m, rop_g, dx, dy, dz)

      ! Add in source terms for DEM drag coupling.
      if (des_continuum_coupled) &
         call gas_drag_u(lo, hi, slo, shi, alo, ahi, A_m, b_m, f_gds, drag_bm, vol)

      call calc_resid_vel (lo, hi, alo, ahi, &
         ulo, uhi, vlo, vhi, wlo, whi, &
         u_g, v_g, w_g, A_m, b_m, &
         resid_u, resid(resid_u,1), resid(resid_u,2), domlo, domhi)

     call under_relax (lo, hi, u_g, ulo, uhi, A_m, b_m, alo, ahi, resid_u)

   end subroutine solve_u_g_star

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: solve_v_g_star                                          !
!  Author: M. Syamlal                                 Date: 25-APR-96  !
!                                                                      !
!  Purpose: Solve starred velocity components                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine solve_v_g_star(&
      slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, lo, hi, &
      u_g, v_g, w_g, v_go, p_g, ro_g, rop_g, &
      rop_go, ep_g, tau_v_g, d_n, fluxX, fluxY, fluxZ, mu_g,  &
      f_gds, A_m, b_m, drag_bm, &
      bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
      bc_klo_type, bc_khi_type, domlo, domhi, dt, dx, dy, dz, resid) &
      bind(C, name="solve_v_g_star")


! Module procedures ...................................................//
      use v_g_conv_dif, only: conv_dif_v_g
      use source_v_g_module, only: source_v_g, source_v_g_bc
      use source_v_g_module, only: point_source_v_g
      use calc_d_mod, only: calc_d_n
      use adjust_a, only: adjust_a_g
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

! Dummy arguments ....................................................//
      integer(c_int)     , intent(in   ) :: slo(3),shi(3)
      integer(c_int)     , intent(in   ) :: ulo(3),uhi(3)
      integer(c_int)     , intent(in   ) :: vlo(3),vhi(3)
      integer(c_int)     , intent(in   ) :: wlo(3),whi(3)
      integer(c_int)     , intent(in   ) :: alo(3),ahi(3)
      integer(c_int)     , intent(in   ) :: domlo(3),domhi(3)
      integer(c_int)     , intent(in   ) ::  lo(3), hi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: fluxX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: v_go&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: fluxY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: tau_v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: fluxZ&
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

      real(c_real), intent(  out) :: d_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(  out) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real), intent(  out) :: resid(8,2)

      integer(c_int), intent(in   ) :: bc_ilo_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jlo_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_klo_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)
      integer(c_int), intent(in   ) :: bc_khi_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)
!.....................................................................//
      real(c_real) :: vol

      vol = dx*dy*dz

      ! Initialize A_m, b_m -- but only on the current tile!
      A_m(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) =  0.0d0
      A_m(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),0) = -1.0d0
      b_m(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))   =  0.0d0

! calculate the convection-diffusion terms
      call conv_dif_v_g (lo, hi, slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
         A_m, mu_g, fluxX, fluxY, fluxZ, dx, dy, dz)

! calculate the source terms for the gas phase u-momentum eqs
      call source_v_g(lo, hi, slo, shi, vlo, vhi, alo, ahi, A_m,&
         b_m, dt, p_g, ep_g, ro_g, rop_go, &
         v_go, tau_v_g, dx, dy, dz, domlo, domhi)

! modifications for bc
      call source_v_g_bc(lo, hi, slo, shi, alo, ahi, A_m, b_m, &
                         bc_ilo_type, bc_ihi_type, &
                         bc_jlo_type, bc_jhi_type, &
                         bc_klo_type, bc_khi_type, &
                         domlo, domhi, dx, dz)

      ! Add in point sources
      if(point_source) call point_source_v_g (lo, hi, alo, ahi, b_m, vol)

      ! Calculate coefficients for the pressure correction equation
      call calc_d_n(lo, hi, slo, shi, vlo, vhi, alo, ahi, d_n, A_m, &
         ep_g, f_gds, dx, dy, dz, domlo, domhi)

      ! Handle special case where center coefficient is zero
      call adjust_a_g('V',slo, shi, alo, ahi, lo, hi, A_m, b_m, rop_g, dx, dy, dz)

      ! Add in source terms for DEM drag coupling.
      if(des_continuum_coupled) &
         call gas_drag_v(lo, hi, slo, shi, alo, ahi, A_m, b_m, f_gds, drag_bm, vol)

      call calc_resid_vel (lo, hi, alo, ahi, &
         vlo, vhi, wlo, whi, ulo, uhi, &
         v_g, w_g, u_g, A_m, b_m, &
         resid_v, resid(resid_v,1), resid(resid_v,2), domlo, domhi)

      call under_relax (lo, hi, v_g, vlo, vhi, A_m, b_m, alo, ahi, resid_v)

   end subroutine solve_v_g_star


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: solve_w_g_star                                          !
!  Author: M. Syamlal                                 Date: 25-APR-96  !
!                                                                      !
!  Purpose: Solve starred velocity components                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine solve_w_g_star(&
      slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, lo, hi, &
      u_g, v_g, w_g, w_go, p_g, ro_g, rop_g, &
      rop_go, ep_g, tau_w_g, d_t, fluxX, fluxY, fluxZ, mu_g,  &
      f_gds, A_m, b_m, drag_bm, &
      bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
      bc_klo_type, bc_khi_type, domlo, domhi, dt, dx, dy, dz, resid) &
      bind(C, name="solve_w_g_star")

! Module procedures ..................................................//
      use w_g_conv_dif, only: conv_dif_w_g
      use source_w_g_module, only: source_w_g, source_w_g_bc
      use source_w_g_module, only: point_source_w_g
      use calc_d_mod, only: calc_d_t
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

! Dummy arguments ....................................................//
      integer(c_int)     , intent(in   ) :: slo(3),shi(3)
      integer(c_int)     , intent(in   ) :: ulo(3),uhi(3)
      integer(c_int)     , intent(in   ) :: vlo(3),vhi(3)
      integer(c_int)     , intent(in   ) :: wlo(3),whi(3)
      integer(c_int)     , intent(in   ) :: alo(3),ahi(3)
      integer(c_int)     , intent(in   ) :: domlo(3),domhi(3)
      integer(c_int)     , intent(in   ) ::  lo(3), hi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: fluxX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: fluxY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), intent(inout) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: w_go&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: fluxZ&
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

      real(c_real), intent(  out) :: d_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(  out) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(  out) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real), intent(  out) :: resid(8,2)

      integer(c_int), intent(in   ) :: bc_ilo_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jlo_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_klo_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)
      integer(c_int), intent(in   ) :: bc_khi_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)
!.....................................................................//
      real(c_real) :: vol

      vol = dx*dy*dz

      ! Initialize A_m, b_m -- but only on the current tile!
      A_m(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) =  0.0d0
      A_m(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),0) = -1.0d0
      b_m(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))   =  0.0d0

      ! calculate the convection-diffusion terms
      call conv_dif_w_g (lo, hi, slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
                         A_m, mu_g, fluxX, fluxY, fluxZ, dx, dy, dz)

      ! calculate the source terms for the gas phase u-momentum eqs
      call source_w_g(lo, hi, slo, shi, wlo, whi, alo, ahi, A_m, b_m, dt, &
         p_g, ep_g, ro_g, rop_go, w_go, tau_w_g, dx, dy, dz, domlo, domhi)

      ! modifications for bc
      call source_w_g_bc (lo, hi, slo, shi, alo, ahi, A_m, b_m, &
                          bc_ilo_type, bc_ihi_type, &
                          bc_jlo_type, bc_jhi_type, &
                          bc_klo_type, bc_khi_type, &
                          domlo, domhi, dx, dy)

      ! Add in point sources
      if(point_source) call point_source_w_g (lo, hi, alo, ahi, b_m, vol)

      ! calculate coefficients for the pressure correction equation
      call calc_d_t(lo, hi, slo, shi, wlo, whi, alo, ahi, &
         d_t, A_m, ep_g, f_gds, dx, dy, dz, domlo, domhi)

      ! handle special case where center coefficient is zero
      call adjust_a_g('W',slo, shi, alo, ahi, lo, hi, A_m, b_m, rop_g, dx, dy, dz)

      ! add in source terms for DEM drag coupling.
      if(des_continuum_coupled) &
         call gas_drag_w(lo, hi, slo, shi, alo, ahi, A_m, b_m, f_gds, drag_bm, vol)

      call calc_resid_vel (lo, hi, alo, ahi, &
         wlo, whi, ulo, uhi, vlo, vhi, &
         w_g, u_g, v_g, A_m, b_m, &
         resid_w, resid(resid_w,1), resid(resid_w,2), domlo, domhi)

      call under_relax (lo, hi, w_g, wlo, whi, A_m, b_m, alo, ahi, resid_w)

   end subroutine solve_w_g_star

end module solve_vel_star_module
