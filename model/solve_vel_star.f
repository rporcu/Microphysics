!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_VEL_STAR                                          C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!                                                                      C
!  Purpose: Solve starred velocity components                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOLVE_VEL_STAR(IER)

      USE adjust_a  , only: adjust_a_g
      USE calc_d_mod, only: calc_d
      USE compar    , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE discretelement, only: des_continuum_coupled
      USE fldvar  , only: d_e, d_n, d_t, u_g, v_g, w_g
      USE leqsol  , only: leq_it, leq_sweep, leq_method, leq_tol, leq_pc 
      USE fun_avg , only: do_k
      USE matrix  , only: a_m, b_m, init_ab_m, lock_ambm, unlock_ambm 
      USE ur_facs , only: under_relax
      USE ps      , only: point_source
      USE residual, only: resid, den_resid, max_resid, ijk_resid, &
                          resid_u, resid_v, resid_w, num_resid
      USE run    , only: momentum_x_eq, momentum_y_eq, momentum_z_eq

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! temporary velocity arrays
      DOUBLE PRECISION, allocatable :: U_gtmp(:,:,:)
      DOUBLE PRECISION, allocatable :: V_gtmp(:,:,:)
      DOUBLE PRECISION, allocatable :: W_gtmp(:,:,:)
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI

!-----------------------------------------------

      allocate(U_gtmp (istart3:iend3, jstart3:jend3, kstart3:kend3))
      allocate(V_gtmp (istart3:iend3, jstart3:jend3, kstart3:kend3))
      allocate(W_gtmp (istart3:iend3, jstart3:jend3, kstart3:kend3))

      call lock_ambm        ! locks arrys a_m and b_m

! Store the velocities so that the order of solving the momentum
! equations does not matter
     U_gtmp = U_g
     V_gtmp = V_g
     W_gtmp = W_g

! Calculate U_m_star and residuals
! ---------------------------------------------------------------->>>
      CALL INIT_AB_M (A_M, B_M)

! calculate the convection-diffusion terms
      CALL CONV_DIF_U_G (A_M, B_M)

! calculate the source terms for the gas and solids phase u-momentum
! equations
      CALL SOURCE_U_G (A_M, B_M)
      IF(POINT_SOURCE) CALL POINT_SOURCE_U_G (A_M, B_M)

! evaluate local variable vxf_gs
      CALL VF_GS_X

! calculate coefficients for the pressure correction equation
      IF (MOMENTUM_X_EQ(0)) THEN
         CALL CALC_D(D_E, "X", A_M)
      ENDIF

! handle special case where center coefficient is zero
      IF (MOMENTUM_X_EQ(0)) CALL ADJUST_A_G ('U', A_M, B_M)

! calculate modifications to the A matrix center coefficient and B
! source vector for treating DEM drag terms
      IF(DES_CONTINUUM_COUPLED) THEN
         CALL GAS_DRAG_U(A_M, B_M, IER)
      ENDIF


      IF (MOMENTUM_X_EQ(0)) THEN
         CALL CALC_RESID_U (U_G, V_G, W_G, A_M, B_M, 0, &
            NUM_RESID(RESID_U), DEN_RESID(RESID_U), &
            RESID(RESID_U), MAX_RESID(RESID_U), &
            IJK_RESID(RESID_U))
         CALL UNDER_RELAX (U_G, A_M, B_M, 'U', 3)
      ENDIF

      IF (MOMENTUM_X_EQ(0)) THEN
         CALL ADJUST_LEQ (RESID(RESID_U), LEQ_IT(3), LEQ_METHOD(3), &
            LEQI, LEQM)
         CALL SOLVE_LIN_EQ ('U_g', 3, U_Gtmp, A_M, B_M, 0, LEQI, LEQM, &
            LEQ_SWEEP(3), LEQ_TOL(3),  LEQ_PC(3), IER)
      ENDIF

! End U_m_star and residuals
! ----------------------------------------------------------------<<<



! Calculate V_m_star and residuals
! ---------------------------------------------------------------->>>
      CALL INIT_AB_M (A_M, B_M)

      CALL CONV_DIF_V_G (A_M, B_M, IER)

      CALL SOURCE_V_G (A_M, B_M)
      IF(POINT_SOURCE) CALL POINT_SOURCE_V_G (A_M, B_M)

      CALL VF_GS_Y

! calculate coefficients for the pressure correction equation
      IF (MOMENTUM_Y_EQ(0)) THEN
         CALL CALC_D(D_N, "Y", A_M)
      ENDIF

      IF (MOMENTUM_Y_EQ(0)) CALL ADJUST_A_G('V',A_M, B_M)

      IF(DES_CONTINUUM_COUPLED) THEN
         CALL GAS_DRAG_V(A_M, B_M, IER)
      ENDIF


      IF (MOMENTUM_Y_EQ(0)) THEN
         CALL CALC_RESID_V (U_G, V_G, W_G, A_M, B_M, 0, &
            NUM_RESID(RESID_V), DEN_RESID(RESID_V), &
            RESID(RESID_V), MAX_RESID(RESID_V), &
            IJK_RESID(RESID_V))
         CALL UNDER_RELAX (V_G, A_M, B_M, 'V', 4)
      ENDIF

      IF (MOMENTUM_Y_EQ(0)) THEN
         CALL ADJUST_LEQ (RESID(RESID_V), LEQ_IT(4), LEQ_METHOD(4), &
            LEQI, LEQM)
         CALL SOLVE_LIN_EQ ('V_g', 4, V_Gtmp, A_M, B_M, 0, LEQI, LEQM, &
            LEQ_SWEEP(4), LEQ_TOL(4),  LEQ_PC(4), IER)
      ENDIF

! End V_m_star and residuals
! ----------------------------------------------------------------<<<


! Calculate W_m_star and residuals
! ---------------------------------------------------------------->>>
      IF (DO_K)THEN
         CALL INIT_AB_M (A_M, B_M)

         CALL CONV_DIF_W_G (A_M, B_M)

         CALL SOURCE_W_G (A_M, B_M)
         IF(POINT_SOURCE) CALL POINT_SOURCE_W_G (A_M, B_M)

         CALL VF_GS_Z

! calculate coefficients for the pressure correction equation
         IF (MOMENTUM_Z_EQ(0)) THEN
            CALL CALC_D(D_T, "Z", A_M)
         ENDIF

         IF (MOMENTUM_Z_EQ(0)) CALL ADJUST_A_G('W',A_M, B_M)

         IF(DES_CONTINUUM_COUPLED) THEN
            CALL GAS_DRAG_W(A_M, B_M, IER)
         ENDIF

         IF (MOMENTUM_Z_EQ(0)) THEN
            CALL CALC_RESID_W (U_G, V_G, W_G, A_M, B_M, 0, &
               NUM_RESID(RESID_W), DEN_RESID(RESID_W), &
               RESID(RESID_W), MAX_RESID(RESID_W), &
               IJK_RESID(RESID_W))
            CALL UNDER_RELAX (W_G, A_M, B_M, 'W', 5)
         ENDIF


         IF (MOMENTUM_Z_EQ(0)) THEN
            CALL ADJUST_LEQ (RESID(RESID_W), LEQ_IT(5), &
               LEQ_METHOD(5), LEQI, LEQM)
            CALL SOLVE_LIN_EQ ('W_g', 5, W_Gtmp, A_M, B_M, 0, LEQI,&
               LEQM, LEQ_SWEEP(5), LEQ_TOL(5), LEQ_PC(5), IER)
!            call out_array(w_g, 'w_g')
         ENDIF

      ENDIF   ! end if (do_k)
! End W_m_star and residuals
! ----------------------------------------------------------------<<<

! Now update all velocity components
      U_g = U_gtmp
      V_g = V_gtmp
      W_g = W_gtmp

      call unlock_ambm

      deallocate(U_gtmp)
      deallocate(V_gtmp)
      deallocate(W_gtmp)

      RETURN
      END SUBROUTINE SOLVE_VEL_STAR
