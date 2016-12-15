module solve_pp_module
   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_Pp_g
!  Purpose: Solve fluid pressure correction equation                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOLVE_PP_G(u_g, v_g, w_g, p_g, ep_g, rop_g, rop_go, &
         ro_g, pp_g, rop_ge, rop_gn, rop_gt, d_e,d_n, d_t, A_m, b_m, &
         flag, NORMG, RESG, IER)

      USE compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE conv_pp_g_module, only: conv_pp_g
      USE leqsol  , only: leq_method, leq_it, leq_sweep, leq_tol, leq_pc
      USE matrix  , only: init_ab_m, lock_ambm, unlock_ambm
      USE param1  , only: zero, one
      USE ps, only: point_source
      USE residual, only: i_resid, j_resid,k_resid, den_resid, max_resid
      USE residual, only: resid_p, resid, num_resid
      USE solve_lin_eq_module, only: solve_lin_eq
      USE source_pp_module, only: source_pp_g
      use residual, only: CALC_RESID_PP

      IMPLICIT NONE
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! Parameter to make tolerance for residual scaled with max value
! compatible with residual scaled with first iteration residual.
! Increase it to tighten convergence.
      DOUBLE PRECISION, PARAMETER :: DEN = 1.0D1   !5.0D2
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------

      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: pp_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: rop_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: rop_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: rop_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: rop_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: d_e&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN  ) :: d_n&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN  ) :: d_t&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(inout) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,-3:3)
      DOUBLE PRECISION, INTENT(inout) :: b_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)

! Normalization factor for gas pressure correction residual.
! At start of the iterate loop normg will either be 1 (i.e. not
! normalized) or a user defined value given by norm_g.  If norm_g
! was set to zero then the normalization is based on dominate
! term in the equation
      DOUBLE PRECISION, INTENT(IN) :: NORMg
! gas pressure correction residual
      DOUBLE PRECISION, INTENT(OUT) :: RESg
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Normalization factor for gas pressure correction residual
      DOUBLE PRECISION :: NORMGloc
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI

! temporary use of global arrays:
! arraym1 (locally b_mmax)
! vector b_m based on dominate term in correction equation
!-----------------------------------------------
! Vector b_m
      DOUBLE PRECISION, allocatable :: B_MMAX(:,:,:)

      ALLOCATE(B_MMAX(istart3:iend3, jstart3:jend3, kstart3:kend3))

      call lock_ambm

! initializing
      PP_G = 0.0d0
      CALL INIT_AB_M (A_m, b_m)

! Forming the sparse matrix equation.
      CALL CONV_PP_G (A_m, rop_ge, rop_gn, rop_gt, flag)

      call source_pp_g(A_m, b_m, B_MMAX, u_g, v_g, w_g, p_g, ep_g,&
         rop_g, rop_go, ro_g, d_e, d_n, d_t, flag)

      IF(POINT_SOURCE) CALL POINT_SOURCE_PP_G (b_m, B_MMAX, flag)

! Find average residual, maximum residual and location
      NORMGloc = NORMG
      IF(NORMG == ZERO) THEN
! calculating the residual based on dominate term in correction equation
! and use this to form normalization factor
        CALL CALC_RESID_PP (B_MMAX, ONE, NUM_RESID(RESID_P), &
         DEN_RESID(RESID_P), RESID(RESID_P), MAX_RESID(RESID_P), &
         i_resid(resid_p),j_resid(resid_p),k_resid(resid_p), flag)
         NORMGloc = RESID(RESID_P)/DEN
      ENDIF
      CALL CALC_RESID_PP (b_m, NORMGloc, NUM_RESID(RESID_P),  &
         DEN_RESID(RESID_P), RESID(RESID_P), MAX_RESID(RESID_P), &
         i_resid(resid_p),j_resid(resid_p),k_resid(resid_p), flag)
      RESG = RESID(RESID_P)

! Solve P_g_prime equation
       LEQI = LEQ_IT(1)
       LEQM = LEQ_METHOD(1)

      CALL SOLVE_LIN_EQ ('Pp_g', 1, PP_G, A_m, b_m, 0, LEQI, LEQM, &
                         LEQ_SWEEP(1), LEQ_TOL(1), LEQ_PC(1), IER)

!      call out_array(Pp_g, 'Pp_g')

      call unlock_ambm

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
      SUBROUTINE POINT_SOURCE_PP_G(b_m, B_mmax,flag)

      use compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use geometry, only: vol
      use param1  , only: small_number
      use ps, only: dimension_ps, ps_defined, ps_massflow_g, ps_volume, ps_i_w, ps_j_s, ps_k_b, ps_i_e, ps_j_n, ps_k_t

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(INOUT) :: B_mmax&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      INTEGER, DIMENSION(:,:,:,:), INTENT(IN) :: FLAG

!-----------------------------------------------
! Local Variables
!-----------------------------------------------

! Indices
      INTEGER :: I, J, K
      INTEGER :: PSV

! terms of bm expression
      DOUBLE PRECISION pSource

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
