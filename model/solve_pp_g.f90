module solve_pp_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOLVE_Pp_g                                              !
!  Purpose: Solve fluid pressure correction equation                   !
!                                                                      !
!  Author: M. Syamlal                                 Date: 19-JUN-96  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine solve_pp_g(slo, shi, lo, hi, &
      u_g, v_g, w_g, p_g, ep_g, rop_g, rop_go, &
      ro_g, rop_ge, rop_gn, rop_gt, d_e,d_n, d_t, a_m, b_m, &
      flag, dt, normg, resg, dx, dy, dz)&
      bind(C, name="solve_pp_g")

! Module procedures ..................................................//
      use matrix  , only: init_ab_m
      use conv_pp_g_module, only: conv_pp_g
      use source_pp_module, only: source_pp_g
      use residual, only: calc_resid_pp

! Global data .........................................................//
! Fluid array bounds
! Flag for existence of point sources
      use ps, only: point_source
! Global data arrays for residuals
      use residual, only: resid_p, resid, max_resid
      use residual, only: num_resid, den_resid
      use residual, only: i_resid, j_resid, k_resid
! parameters, 0.0 and 1.0
      use param1, only: zero, one

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      real(c_real), intent(in   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
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
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: d_e&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: d_n&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: d_t&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(  out) :: a_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Normalization factor for gas pressure correction residual.
! At start of the iterate loop normg will either be 1 (i.e. not
! normalized) or a user defined value given by norm_g.  If norm_g
! was set to zero then the normalization is based on dominate
! term in the equation
      real(c_real), intent(in) :: normg
! gas pressure correction residual
      real(c_real), intent(out) :: resg

! Local parameters ...................................................//
! Parameter to make tolerance for residual scaled with max value
! compatible with residual scaled with first iteration residual.
! Increase it to tighten convergence.
      real(c_real), PARAMETER :: DEN = 1.0D1   !5.0D2
! Normalization factor for gas pressure correction residual
      real(c_real) :: NORMGloc
! dominate term in correction equation max(am,bm)
      real(c_real), allocatable :: B_MMAX(:,:,:)
!.....................................................................//

      ALLOCATE( B_MMAX(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )

! initializing
      call init_ab_m (slo, shi, a_m, b_m)

! Forming the sparse matrix equation.
      call conv_pp_g (slo, shi, lo, hi, a_m, rop_ge, rop_gn, rop_gt, flag, dx, dy, dz)

      call source_pp_g(slo, shi, lo, hi, a_m, b_m, b_mmax, dt, u_g, v_g, w_g, p_g, ep_g,&
         rop_g, rop_go, ro_g, d_e, d_n, d_t, flag, dx, dy, dz)

      if(point_source) call point_source_pp_g (slo, shi, lo, hi, b_m, b_mmax, flag, dx, dy, dz)

! Find average residual, maximum residual and location
      normgloc = normg
      if(abs(normg) < epsilon(zero)) then
! calculating the residual based on dominate term in correction equation
! and use this to form normalization factor
        call calc_resid_pp (slo, shi, lo, hi, &
         b_mmax, one, num_resid(resid_p), &
         den_resid(resid_p), resid(resid_p), max_resid(resid_p), &
         i_resid(resid_p),j_resid(resid_p),k_resid(resid_p), flag)
         normgloc = resid(resid_p)/den
      endif

      call calc_resid_pp (slo, shi, lo, hi, &
         b_m, normgloc, num_resid(resid_p),  &
         den_resid(resid_p), resid(resid_p), max_resid(resid_p), &
         i_resid(resid_p),j_resid(resid_p),k_resid(resid_p), flag)
      resg = resid(resid_p)


      RETURN
      END SUBROUTINE SOLVE_PP_G


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
      SUBROUTINE POINT_SOURCE_PP_G(slo, shi, lo, hi, b_m, B_mmax,flag, dx, dy, dz)

      use param1  , only: small_number
      use ps, only: dimension_ps, ps_defined, ps_massflow_g, ps_volume,&
         ps_i_w, ps_j_s, ps_k_b, ps_i_e, ps_j_n, ps_k_t

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Vector b_m
      real(c_real), INTENT(INOUT) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! maximum term in b_m expression
      real(c_real), INTENT(INOUT) :: B_mmax&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      INTEGER, DIMENSION(:,:,:,:), INTENT(IN) :: FLAG

      real(c_real), INTENT(IN   ) :: dx,dy,dz

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

         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         if(PS_MASSFLOW_G(PSV) < small_number) cycle PS_LP

         do k = PS_K_B(PSV), PS_K_T(PSV)
         do j = PS_J_S(PSV), PS_J_N(PSV)
         do i = PS_I_W(PSV), PS_I_E(PSV)

            if(1.eq.flag(i,j,k,1)) then
               pSource = PS_MASSFLOW_G(PSV) * (VOL/PS_VOLUME(PSV))

               b_m(I,J,K) = b_m(I,J,K) - pSource
               B_MMAX(I,J,K) = max(abs(B_MMAX(I,J,K)), abs(b_m(I,J,K)))
            endif

         enddo
         enddo
         enddo

      enddo PS_LP

      END SUBROUTINE POINT_SOURCE_PP_G
end module solve_pp_module
