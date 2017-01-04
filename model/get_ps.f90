!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: GET_PS                                                  !
!  Author: J.Musser                                   Date: 19-MAR-14  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for PS's               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GET_PS(PSV)

      USE calc_cell_module, only: calc_loc, calc_cell
      USE geometry, only: dx, dy, dz
      USE geometry, only: imax, jmax, kmax
      USE geometry, only: imax2, jmax2, kmax2
      USE location_check_module, only: location_check
      USE param1, only: undefined_i, is_undefined, is_defined
      USE ps, only: ps_x_w, ps_x_e, ps_i_w, ps_i_e
      USE ps, only: ps_y_n, ps_y_s, ps_j_s, ps_j_n
      USE ps, only: ps_z_t, ps_z_b, ps_k_b, ps_k_t

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Loop/variable indices
      INTEGER, INTENT(in) :: PSV

! Local Variables:
!---------------------------------------------------------------------//
! Error flag.
      INTEGER :: IER
! Calculated indices of the wall boundary
      INTEGER :: I_w , I_e , J_s , J_n , K_b , K_t
! Surface indictors
      LOGICAL :: X_CONSTANT, Y_CONSTANT, Z_CONSTANT
!......................................................................!

      CALL INIT_ERR_MSG('GET_PS')

      X_CONSTANT = .TRUE.
      Y_CONSTANT = .TRUE.
      Z_CONSTANT = .TRUE.

      IF(IS_DEFINED(PS_X_W(PSV)) .AND. IS_DEFINED(PS_X_E(PSV))) THEN
         CALL CALC_CELL(PS_X_W(PSV), DX, IMAX, I_W)
         CALL CALC_CELL(PS_X_E(PSV), DX, IMAX, I_E)
         IF (PS_X_W(PSV) /= PS_X_E(PSV)) THEN
            X_CONSTANT = .FALSE.
            I_W = I_W + 1
            IF(PS_I_W(PSV)/=UNDEFINED_I .OR.                           &
               PS_I_E(PSV)/=UNDEFINED_I) THEN
               CALL LOCATION_CHECK(PS_I_W(PSV), I_W, PSV, 'PS - west')
               CALL LOCATION_CHECK(PS_I_E(PSV), I_E, PSV, 'PS - east')
            ENDIF
         ENDIF
         PS_I_W(PSV) = I_W
         PS_I_E(PSV) = I_E
      ELSE
         IF(PS_I_W(PSV) /= UNDEFINED_I)                                &
            CALL CALC_LOC(DX, PS_I_W(PSV), PS_X_W(PSV))
         IF(PS_I_E(PSV) /= UNDEFINED_I)                                &
            CALL CALC_LOC(DX, PS_I_E(PSV), PS_X_E(PSV))
         IF (PS_X_W(PSV) /= PS_X_E(PSV)) X_CONSTANT = .FALSE.
      ENDIF

!  If there is no variation in the I direction set indices to 1
      IF (IS_DEFINED(PS_Y_S(PSV)).AND. IS_DEFINED(PS_Y_N(PSV))) THEN
         CALL CALC_CELL(PS_Y_S(PSV), DY, JMAX, J_S)
         CALL CALC_CELL(PS_Y_N(PSV), DY, JMAX, J_N)
         IF (PS_Y_S(PSV) /= PS_Y_N(PSV)) THEN
            Y_CONSTANT = .FALSE.
            J_S = J_S + 1
            IF(PS_J_S(PSV)/=UNDEFINED_I .OR.                           &
               PS_J_N(PSV)/=UNDEFINED_I) THEN
               CALL LOCATION_CHECK(PS_J_S(PSV), J_S, PSV, 'PS - south')
               CALL LOCATION_CHECK(PS_J_N(PSV), J_N, PSV, 'PS - north')
            ENDIF
         ENDIF
         PS_J_S(PSV) = J_S
         PS_J_N(PSV) = J_N
      ELSE
         IF(PS_J_S(PSV) /= UNDEFINED_I)                                &
            CALL CALC_LOC(DY, PS_J_S(PSV), PS_Y_S(PSV))
         IF(PS_J_N(PSV) /= UNDEFINED_I)                                &
            CALL CALC_LOC(DY, PS_J_N(PSV), PS_Y_N(PSV))
         IF (PS_Y_S(PSV) /= PS_Y_N(PSV)) Y_CONSTANT = .FALSE.
      ENDIF

      IF (IS_DEFINED(PS_Z_B(PSV)) .AND. IS_DEFINED(PS_Z_T(PSV))) THEN
         CALL CALC_CELL(PS_Z_B(PSV), DZ, KMAX, K_B)
         CALL CALC_CELL(PS_Z_T(PSV), DZ, KMAX, K_T)
         IF (PS_Z_B(PSV) /= PS_Z_T(PSV)) THEN
            Z_CONSTANT = .FALSE.
            K_B = K_B + 1
            IF (PS_K_B(PSV)/=UNDEFINED_I .OR.                          &
               PS_K_T(PSV)/=UNDEFINED_I) THEN
               CALL LOCATION_CHECK(PS_K_B(PSV), K_B, PSV, 'PS - bottom')
               CALL LOCATION_CHECK(PS_K_T(PSV), K_T, PSV, 'PS - top')
            ENDIF
         ENDIF
         PS_K_B(PSV) = K_B
         PS_K_T(PSV) = K_T
      ELSE
         IF(PS_K_B(PSV) /= UNDEFINED_I)                                &
            CALL CALC_LOC(DZ, PS_K_B(PSV), PS_Z_B(PSV))
         IF(PS_K_T(PSV) /= UNDEFINED_I)                                &
            CALL CALC_LOC(DZ, PS_K_T(PSV), PS_Z_T(PSV))
         IF (PS_Z_B(PSV) /= PS_Z_T(PSV)) Z_CONSTANT = .FALSE.
      ENDIF

! CHECK FOR VALID VALUES
      IER = 0
      IF(PS_I_W(PSV)<1 .OR. PS_I_W(PSV)>IMAX2) IER = 1
      IF(PS_I_E(PSV)<1 .OR. PS_I_E(PSV)>IMAX2) IER = 1
      IF(PS_J_S(PSV)<1 .OR. PS_J_S(PSV)>JMAX2) IER = 1
      IF(PS_J_N(PSV)<1 .OR. PS_J_N(PSV)>JMAX2) IER = 1
      IF(PS_K_B(PSV)<1 .OR. PS_K_B(PSV)>KMAX2) IER = 1
      IF(PS_K_T(PSV)<1 .OR. PS_K_T(PSV)>KMAX2) IER = 1
      IF(PS_K_B(PSV) > PS_K_T(PSV)) IER = 1
      IF(PS_J_S(PSV) > PS_J_N(PSV)) IER = 1
      IF(PS_I_W(PSV) > PS_I_E(PSV)) IER = 1

      IF(IER /= 0)THEN
         WRITE(ERR_MSG,1101) PSV,                                      &
            'X', PS_X_W(PSV), PS_X_E(PSV),'I',PS_I_W(PSV),PS_I_E(PSV), &
            'Y', PS_Y_S(PSV), PS_Y_N(PSV),'J',PS_J_S(PSV),PS_J_N(PSV), &
            'Z', PS_Z_B(PSV), PS_Z_T(PSV),'K',PS_K_B(PSV),PS_K_T(PSV)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1101 FORMAT('Error 1101: Invalid location specified for PS ',I3,'.',  &
         3(/3x,A1,': ',g12.5,',',g12.5,8x,A1,': ',I8,',',I8),/         &
         'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE GET_PS
