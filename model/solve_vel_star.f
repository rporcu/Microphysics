!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_VEL_STAR                                          C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!                                                                      C
!  Purpose: Solve starred velocity components                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOLVE_VEL_STAR(IER)

      USE adjust_a
      USE compar
      USE discretelement
      USE drag
      USE fldvar
      USE geometry
      USE indices
      USE leqsol
      USE output
      USE param
      USE param1
      USE physprop
      USE ps
      USE residual
      USE run
      USE toleranc
      USE ur_facs
      use matrix

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! solids phase index
      INTEGER :: M
! fluid cell index
      INTEGER :: IJK
! temporary velocity arrays
      DOUBLE PRECISION, allocatable :: U_gtmp(:)
      DOUBLE PRECISION, allocatable :: V_gtmp(:)
      DOUBLE PRECISION, allocatable :: W_gtmp(:)
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI

!-----------------------------------------------

      allocate(U_gtmp(DIMENSION_3))
      allocate(V_gtmp(DIMENSION_3))
      allocate(W_gtmp(DIMENSION_3))

      call lock_ambm        ! locks arrys a_m and b_m

! Store the velocities so that the order of solving the momentum
! equations does not matter
      DO IJK = ijkstart3, ijkend3
         U_gtmp(IJK) = U_g(IJK)
         V_gtmp(IJK) = V_g(IJK)
         W_gtmp(IJK) = W_g(IJK)
      ENDDO


! Calculate U_m_star and residuals
! ---------------------------------------------------------------->>>
      CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0)

! calculate the convection-diffusion terms
      CALL CONV_DIF_U_G (A_M, B_M)

! calculate the source terms for the gas and solids phase u-momentum
! equations
      CALL SOURCE_U_G (A_M, B_M)
      IF(POINT_SOURCE) CALL POINT_SOURCE_U_G (A_M, B_M)

! evaluate local variable vxf_gs
      CALL VF_GS_X

! calculate coefficients for the pressure correction equation
      CALL CALC_D_E (A_M, IER)

! handle special case where center coefficient is zero
      IF (MOMENTUM_X_EQ(0)) CALL ADJUST_A_G ('U', A_M, B_M)

! calculate modifications to the A matrix center coefficient and B
! source vector for treating DEM drag terms
      IF(DES_CONTINUUM_COUPLED) THEN
         CALL GAS_DRAG_U(A_M, B_M, IER)
      ENDIF


      IF (MOMENTUM_X_EQ(0)) THEN
         CALL CALC_RESID_U (U_G, V_G, W_G, A_M, B_M, 0, &
            NUM_RESID(RESID_U,0), DEN_RESID(RESID_U,0), &
            RESID(RESID_U,0), MAX_RESID(RESID_U,0), &
            IJK_RESID(RESID_U,0))
         CALL UNDER_RELAX_U (U_G, A_M, B_M, 0, UR_FAC(3))
      ENDIF

      IF (MOMENTUM_X_EQ(0)) THEN
         CALL ADJUST_LEQ (RESID(RESID_U,0), LEQ_IT(3), LEQ_METHOD(3), &
            LEQI, LEQM)
         CALL SOLVE_LIN_EQ ('U_g', 3, U_Gtmp, A_M, B_M, 0, LEQI, LEQM, &
            LEQ_SWEEP(3), LEQ_TOL(3),  LEQ_PC(3), IER)
      ENDIF

! End U_m_star and residuals
! ----------------------------------------------------------------<<<



! Calculate V_m_star and residuals
! ---------------------------------------------------------------->>>
      CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0)

      CALL CONV_DIF_V_G (A_M, B_M, IER)

      CALL SOURCE_V_G (A_M, B_M)
      IF(POINT_SOURCE) CALL POINT_SOURCE_V_G (A_M, B_M)

      CALL VF_GS_Y

! calculate coefficients for the pressure correction equation
      CALL CALC_D_N (A_M, IER)

      IF (MOMENTUM_Y_EQ(0)) CALL ADJUST_A_G('V',A_M, B_M)

      IF(DES_CONTINUUM_COUPLED) THEN
         CALL GAS_DRAG_V(A_M, B_M, IER)
      ENDIF


      IF (MOMENTUM_Y_EQ(0)) THEN
         CALL CALC_RESID_V (U_G, V_G, W_G, A_M, B_M, 0, &
            NUM_RESID(RESID_V,0), DEN_RESID(RESID_V,0), &
            RESID(RESID_V,0), MAX_RESID(RESID_V,0), &
            IJK_RESID(RESID_V,0))
         CALL UNDER_RELAX_V (V_G, A_M, B_M, 0, UR_FAC(4))
      ENDIF

      IF (MOMENTUM_Y_EQ(0)) THEN
         CALL ADJUST_LEQ (RESID(RESID_V,0), LEQ_IT(4), LEQ_METHOD(4), &
            LEQI, LEQM)
         CALL SOLVE_LIN_EQ ('V_g', 4, V_Gtmp, A_M, B_M, 0, LEQI, LEQM, &
            LEQ_SWEEP(4), LEQ_TOL(4),  LEQ_PC(4), IER)
      ENDIF

! End V_m_star and residuals
! ----------------------------------------------------------------<<<


! Calculate W_m_star and residuals
! ---------------------------------------------------------------->>>
      IF (DO_K)THEN
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0)

         CALL CONV_DIF_W_G (A_M, B_M)

         CALL SOURCE_W_G (A_M, B_M)
         IF(POINT_SOURCE) CALL POINT_SOURCE_W_G (A_M, B_M)

         CALL VF_GS_Z

! calculate coefficients for the pressure correction equation
         CALL CALC_D_T (A_M, IER)

         IF (MOMENTUM_Z_EQ(0)) CALL ADJUST_A_G('W',A_M, B_M)

         IF(DES_CONTINUUM_COUPLED) THEN
            CALL GAS_DRAG_W(A_M, B_M, IER)
         ENDIF

         IF (MOMENTUM_Z_EQ(0)) THEN
            CALL CALC_RESID_W (U_G, V_G, W_G, A_M, B_M, 0, &
               NUM_RESID(RESID_W,0), DEN_RESID(RESID_W,0), &
               RESID(RESID_W,0), MAX_RESID(RESID_W,0), &
               IJK_RESID(RESID_W,0))
            CALL UNDER_RELAX_W (W_G, A_M, B_M, 0, UR_FAC(5))
         ENDIF


         IF (MOMENTUM_Z_EQ(0)) THEN
!            call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), &
!               1, DO_K, ier)
            CALL ADJUST_LEQ (RESID(RESID_W,0), LEQ_IT(5), &
               LEQ_METHOD(5), LEQI, LEQM)
            CALL SOLVE_LIN_EQ ('W_g', 5, W_Gtmp, A_M, B_M, 0, LEQI,&
               LEQM, LEQ_SWEEP(5), LEQ_TOL(5), LEQ_PC(5), IER)
!            call out_array(w_g, 'w_g')
         ENDIF

      ENDIF   ! end if (do_k)
! End W_m_star and residuals
! ----------------------------------------------------------------<<<

! Now update all velocity components
      DO IJK = ijkstart3, ijkend3
         U_g(IJK) = U_gtmp(IJK)
         V_g(IJK) = V_gtmp(IJK)
         W_g(IJK) = W_gtmp(IJK)
      ENDDO


      call unlock_ambm

      deallocate(U_gtmp)
      deallocate(V_gtmp)
      deallocate(W_gtmp)

      RETURN
      END SUBROUTINE SOLVE_VEL_STAR
