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

      SUBROUTINE SOLVE_PP_G(u_g, v_g, w_g, p_g, ep_g, rop_g, rop_go, ro_g, pp_g, &
                            rop_ge, rop_gn, rop_gt, d_e,d_n, d_t, NORMG, RESG, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE param1  , only: zero, one
      USE residual, only: resid_p, resid, num_resid
      USE residual, only: ijk_resid, den_resid, max_resid
      USE leqsol  , only: leq_method, leq_it, leq_sweep, leq_tol, leq_pc
      USE run
      use matrix  , only: a_m, b_m, init_ab_m, lock_ambm, unlock_ambm
      use ps

      use source_pp_module

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

      DOUBLE PRECISION, INTENT(INOUT) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: pp_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: d_e&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: d_n&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: d_t&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

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
! vector B_M based on dominate term in correction equation
!-----------------------------------------------
! Vector b_m
      DOUBLE PRECISION, allocatable :: B_MMAX(:,:,:)

      ALLOCATE(B_MMAX(istart3:iend3, jstart3:jend3, kstart3:kend3))

      call lock_ambm

! initializing
      PP_G = 0.0d0
      CALL INIT_AB_M (A_M, B_M)

! Forming the sparse matrix equation.
      CALL CONV_PP_G (A_M, B_M, rop_ge, rop_gn, rop_gt)

      call source_pp_g(A_M, B_M, B_MMAX, u_g, v_g, w_g, p_g, ep_g,&
                       rop_g, rop_go, ro_g, d_e, d_n, d_t)

      IF(POINT_SOURCE) CALL POINT_SOURCE_PP_G (B_M, B_MMAX)

! Find average residual, maximum residual and location
      NORMGloc = NORMG
      IF(NORMG == ZERO) THEN
! calculating the residual based on dominate term in correction equation
! and use this to form normalization factor
        CALL CALC_RESID_PP (B_MMAX, ONE, NUM_RESID(RESID_P), &
         DEN_RESID(RESID_P), RESID(RESID_P), MAX_RESID(RESID_P), &
         IJK_RESID(RESID_P))
         NORMGloc = RESID(RESID_P)/DEN
      ENDIF
      CALL CALC_RESID_PP (B_M, NORMGloc, NUM_RESID(RESID_P),  &
         DEN_RESID(RESID_P), RESID(RESID_P), MAX_RESID(RESID_P), &
         IJK_RESID(RESID_P))
      RESG = RESID(RESID_P)

!      write(*,*) resid(resid_p), max_resid(resid_p), &
!         ijk_resid(resid_p)

! Solve P_g_prime equation
       LEQI = LEQ_IT(1)
       LEQM = LEQ_METHOD(1)

      CALL SOLVE_LIN_EQ ('Pp_g', 1, PP_G, A_M, B_M, 0, LEQI, LEQM, &
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
      SUBROUTINE POINT_SOURCE_PP_G(B_M, B_mmax)

      use compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use geometry, only: vol
      use param1  , only: small_number
      use ps
      use run
      use functions, only: fluid_at

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

            if(fluid_at(i,j,k)) then
               pSource = PS_MASSFLOW_G(PSV) * (VOL/PS_VOLUME(PSV))

               B_M(I,J,K) = B_M(I,J,K) - pSource
               B_MMAX(I,J,K) = max(abs(B_MMAX(I,J,K)), abs(B_M(I,J,K)))
            endif

         enddo
         enddo
         enddo

      enddo PS_LP

      END SUBROUTINE POINT_SOURCE_PP_G
end module solve_pp_module
