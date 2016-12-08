module solve_vel_star_module
   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_VEL_STAR                                          C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!                                                                      C
!  Purpose: Solve starred velocity components                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOLVE_VEL_STAR(u_g, v_g, w_g, u_go, v_go, w_go,    &
         p_g, ro_g, rop_g, rop_go, ep_g, tau_u_g, tau_v_g, tau_w_g, &
         d_e, d_n, d_t, flux_ge, flux_gn, flux_gt ,mu_g, &
         f_gds, drag_am, drag_bm, flag, IER)

      USE u_g_conv_dif, only: conv_dif_u_g
      USE v_g_conv_dif, only: conv_dif_v_g
      USE w_g_conv_dif, only: conv_dif_w_g

      USE adjust_a  , only: adjust_a_g
      USE calc_d_mod, only: calc_d
      USE compar    , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE discretelement, only: des_continuum_coupled
      USE leqsol  , only: leq_it, leq_sweep, leq_method, leq_tol, leq_pc
      USE matrix  , only: a_m, b_m, init_ab_m, lock_ambm, unlock_ambm
      USE ur_facs , only: under_relax
      USE ps      , only: point_source
      USE run    , only: momentum_x_eq, momentum_y_eq, momentum_z_eq

      USE source_u_g_module, only: source_u_g, point_source_u_g
      USE source_v_g_module, only: source_v_g, point_source_v_g
      USE source_w_g_module, only: source_w_g, point_source_w_g
      use residual, only: CALC_RESID_VEL
      USE residual, only: resid, den_resid, max_resid, i_resid, j_resid, k_resid,&
                          resid_u, resid_v, resid_w, num_resid
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
      DOUBLE PRECISION, INTENT(INOUT) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: rop_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: tau_u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: tau_v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: tau_w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: d_e&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: d_n&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: d_t&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: flux_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: flux_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: flux_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: drag_am&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
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
      CALL CONV_DIF_U_G (A_M, MU_G, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt, flag)

      ! Calculate the source terms for the gas and solids phase u-momentum eqs
      CALL SOURCE_U_G(A_M, B_M, p_g, ep_g, ro_g, rop_g, rop_go, &
                      u_g, u_go, tau_u_g, flag)

      IF(POINT_SOURCE) CALL POINT_SOURCE_U_G (A_M, B_M)

! evaluate local variable vxf_gs
! calculate coefficients for the pressure correction equation
      IF (MOMENTUM_X_EQ(0)) THEN
         CALL CALC_D(D_E, "X", A_M, ep_g, f_gds, flag)
      ENDIF

! handle special case where center coefficient is zero
      IF (MOMENTUM_X_EQ(0)) CALL adjust_a_g ('U', A_M, B_M, ROP_G)

! calculate modifications to the A matrix center coefficient and B
! source vector for treating DEM drag terms
      IF(DES_CONTINUUM_COUPLED) &
         CALL GAS_DRAG_U(A_M, B_M, f_gds, drag_am, drag_bm, IER)

      IF (MOMENTUM_X_EQ(0)) THEN
         CALL CALC_RESID_VEL (U_G, V_G, W_G, A_M, B_M, 0, &
            NUM_RESID(RESID_U), DEN_RESID(RESID_U), &
            RESID(RESID_U), MAX_RESID(RESID_U), &
            i_resid(RESID_U),j_resid(RESID_U),k_resid(RESID_U), flag)
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

      CALL CONV_DIF_V_G (A_M, MU_G, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt, flag)

      CALL SOURCE_V_G(A_M, B_M, p_g, ep_g, ro_g, rop_g, rop_go, &
                      v_g, v_go, tau_v_g, flag)
      IF(POINT_SOURCE) CALL POINT_SOURCE_V_G (A_M, B_M)

! calculate coefficients for the pressure correction equation
      IF (MOMENTUM_Y_EQ(0)) THEN
         CALL CALC_D(D_N, "Y", A_M, ep_g, f_gds, flag)
      ENDIF

      IF (MOMENTUM_Y_EQ(0)) CALL adjust_a_g('V',A_M, B_M, ROP_G)

      IF(DES_CONTINUUM_COUPLED) &
         CALL GAS_DRAG_V(A_M, B_M, f_gds, drag_am, drag_bm, IER)

      IF (MOMENTUM_Y_EQ(0)) THEN
         ! Note we pass V first since that is the primary velocity component
         !   The order of the other two components doesn't matter
         CALL CALC_RESID_VEL (V_G, W_G, U_G, A_M, B_M, 0, &
            NUM_RESID(RESID_V), DEN_RESID(RESID_V), &
            RESID(RESID_V), MAX_RESID(RESID_V), &
            i_resid(RESID_V),j_resid(RESID_V),k_resid(RESID_V), flag)
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
      CALL INIT_AB_M (A_M, B_M)

      CALL CONV_DIF_W_G (A_M, MU_G, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt, flag)

      CALL SOURCE_W_G(A_M, B_M, p_g, ep_g, ro_g, rop_g, rop_go, &
                         w_g, w_go, tau_w_g, flag)
      IF(POINT_SOURCE) CALL POINT_SOURCE_W_G (A_M, B_M)

! calculate coefficients for the pressure correction equation
      IF (MOMENTUM_Z_EQ(0)) THEN
         CALL CALC_D(D_T, "Z", A_M, ep_g, f_gds, flag)
      ENDIF

      IF (MOMENTUM_Z_EQ(0)) CALL adjust_a_g('W',A_M, B_M, ROP_G)

      IF(DES_CONTINUUM_COUPLED) &
         CALL GAS_DRAG_W(A_M, B_M, f_gds, drag_am, drag_bm, IER)

      IF (MOMENTUM_Z_EQ(0)) THEN
            ! Note we pass W first since that is the primary velocity component
         CALL CALC_RESID_VEL (W_G, U_G, V_G, A_M, B_M, 0, &
            NUM_RESID(RESID_W), DEN_RESID(RESID_W), &
            RESID(RESID_W), MAX_RESID(RESID_W), &
            i_resid(RESID_W),j_resid(RESID_W),k_resid(RESID_W),flag)
         CALL UNDER_RELAX (W_G, A_M, B_M, 'W', 5)
      ENDIF


      IF (MOMENTUM_Z_EQ(0)) THEN
         CALL ADJUST_LEQ (RESID(RESID_W), LEQ_IT(5), &
            LEQ_METHOD(5), LEQI, LEQM)
         CALL SOLVE_LIN_EQ ('W_g', 5, W_Gtmp, A_M, B_M, 0, LEQI,&
            LEQM, LEQ_SWEEP(5), LEQ_TOL(5), LEQ_PC(5), IER)
!         call out_array(w_g, 'w_g')
      ENDIF

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

      END SUBROUTINE SOLVE_VEL_STAR
end module solve_vel_star_module
