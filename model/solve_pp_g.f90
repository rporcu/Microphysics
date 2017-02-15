module solve_pp_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   implicit none

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: solve_pp_g                                              !
!  Purpose: Solve fluid pressure correction equation                   !
!                                                                      !
!  Author: M. Syamlal                                 Date: 19-JUN-96  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine solve_pp_g(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, lo, hi, &
      u_g, v_g, w_g, p_g, ep_g, rop_g, rop_go, &
      ro_g, rop_ge, rop_gn, rop_gt, d_e,d_n, d_t, A_m, b_m, b_mmax, &
      bc_ilo_type, bc_ihi_type, bc_jlo_type, &
      bc_jhi_type, bc_klo_type, bc_khi_type, &
      dt, dx, dy, dz)&
      bind(C, name="solve_pp_g")

! Module procedures ..................................................//
      use conv_pp_g_module, only: conv_pp_g
      use source_pp_module, only: source_pp_g
      use source_pp_module, only: source_pp_g_bc

! Global data .........................................................//
! Fluid array bounds

      ! Flag for existence of point sources
      use ps, only: point_source

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: alo(3),ahi(3)
      real(c_real)  , intent(in   ) :: dt, dx, dy, dz

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: rop_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: rop_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: rop_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: d_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: d_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: d_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer(c_int), intent(in   ) :: bc_ilo_type&
         (slo(2):shi(2),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
         (slo(2):shi(2),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_jlo_type&
         (slo(1):shi(1),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
         (slo(1):shi(1),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_klo_type&
         (slo(1):shi(1),slo(2):shi(2),2)
      integer(c_int), intent(in   ) :: bc_khi_type&
         (slo(1):shi(1),slo(2):shi(2),2)

      real(c_real), intent(  out) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(  out) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real), intent(  out) :: b_mmax&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      integer :: i,j,k

      ! Initialize A_m and b_m
      A_m(:,:,:,:)  =  0.0d0
      A_m(:,:,:,0)  = -1.0d0
      b_m(:,:,:)    =  0.0d0

      ! Forming the sparse matrix equation.
      call conv_pp_g (ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
         A_m, rop_ge, rop_gn, rop_gt, dx, dy, dz)

      call source_pp_g(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, lo, hi, &
         A_m, b_m, b_mmax, dt, u_g, v_g, w_g, p_g, ep_g,&
         rop_g, rop_go, ro_g, d_e, d_n, d_t, dx, dy, dz)

      call source_pp_g_bc(slo, shi, alo, ahi, A_m, bc_ilo_type, bc_ihi_type, &
         bc_jlo_type, bc_jhi_type, bc_klo_type, bc_khi_type)

      if (point_source) call point_source_pp_g (alo, ahi, b_m, b_mmax, dx, dy, dz)

      end subroutine solve_pp_g

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_Pp_g                                       C
!  Purpose: Adds point sources to the Pressure correction equation.    C
!                                                                      C
!  Notes: The off-diagonal coefficients are positive. The center       C
!         coefficient and the source vector are negative. See          C
!         conv_Pp_g                                                    C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_PP_G(alo, ahi, b_m, B_mmax, dx, dy, dz)

      use param1  , only: small_number
      use ps, only: dimension_ps, ps_defined, ps_massflow_g, ps_volume,&
         ps_i_w, ps_j_s, ps_k_b, ps_i_e, ps_j_n, ps_k_t

      integer(c_int), intent(in   ) :: alo(3),ahi(3)
      real(c_real)  , intent(in   ) :: dx,dy,dz

      ! Vector b_m
      real(c_real), intent(INOUT) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      ! maximum term in b_m expression
      real(c_real), intent(INOUT) :: b_mmax&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

!-----------------------------------------------
! Local Variables
!-----------------------------------------------

! Indices
      INTEGER :: I, J, K
      INTEGER :: PSV

! terms of bm expression
      real(c_real) pSource
      real(c_real) vol

      vol = dx*dy*dz

!-----------------------------------------------
      PS_LP: do PSV = 1, DIMENSION_PS

         if(.not.ps_defined(psv)) cycle ps_lp
         if(ps_massflow_g(psv) < small_number) cycle ps_lp

         do k = PS_K_B(PSV), PS_K_T(PSV)
            do j = PS_J_S(PSV), PS_J_N(PSV)
               do i = PS_I_W(PSV), PS_I_E(PSV)

                  pSource = PS_MASSFLOW_G(PSV) * (VOL/PS_VOLUME(PSV))

                  b_m(I,J,K) = b_m(I,J,K) - pSource
                  B_MMAX(I,J,K) = max(abs(B_MMAX(I,J,K)), abs(b_m(I,J,K)))

               enddo
            enddo
         enddo

      enddo PS_LP

      END SUBROUTINE POINT_SOURCE_PP_G
end module solve_pp_module
