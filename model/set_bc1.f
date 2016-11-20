!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_BC1                                                 C
!  Purpose: Set transient flow boundary conditions                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC1

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_defined, bc_type
      USE param, only: dimension_bc
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! index for boundary condition
      INTEGER :: L


! Set the boundary conditions
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

            SELECT CASE(TRIM(BC_TYPE(L)))
            CASE ('P_OUTFLOW')
               CALL SET_OUTFLOW(L)
               CALL SET_BC1_REPORT_OUTFLOW(L)
            CASE ('MASS_OUTFLOW')
               CALL SET_OUTFLOW(L)
               CALL SET_BC1_ADJUST_OUTFLOW(L)
            CASE ('MASS_INFLOW')
               CALL SET_BC1_JET(L)
            CASE ('P_INFLOW')
               CALL SET_OUTFLOW(L)
            CASE ('OUTFLOW')
               CALL SET_OUTFLOW(L)
               CALL SET_BC1_REPORT_OUTFLOW(L)
            END SELECT
         ENDIF   ! end if (bc_defined(l))
      ENDDO    ! end do loop (l=1,dimension_bc)

      RETURN
      END SUBROUTINE SET_BC1



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc1_jet                                             C
!  Purpose: update transient jet conditions                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC1_JET(BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_plane
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_i_w, bc_i_e
      use bc, only: bc_jet_g
      use bc, only: bc_time
      use bc, only: bc_jet_gh, bc_dt_h
      use bc, only: bc_jet_gl, bc_dt_l

      use fldvar, only: u_g, v_g, w_g

      use run, only: time, dt

      use param1, only: undefined

      use functions, only: is_on_mype_plus2layers
      use functions, only: funijk, iminus, jminus, kminus
      use compar, only: dead_cell_at

      IMPLICIT NONE
! Dummy arguments
!---------------------------------------------------------------------//
! index for boundary condition
      INTEGER, INTENT(IN) :: BCV

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I, J, K, IJK
! IJK index for setting velocity bc
      INTEGER :: IJK2
!---------------------------------------------------------------------//

      IF (TIME + 0.1d0*DT>=BC_TIME(BCV) .AND. &
          BC_JET_G(BCV)/=UNDEFINED) THEN

         IF (BC_JET_G(BCV) == BC_JET_GH(BCV)) THEN
            BC_JET_G(BCV) = BC_JET_GL(BCV)
            BC_TIME(BCV) = TIME + BC_DT_L(BCV)
         ELSEIF (BC_JET_G(BCV) == BC_JET_GL(BCV)) THEN
            BC_JET_G(BCV) = BC_JET_GH(BCV)
            BC_TIME(BCV) = TIME + BC_DT_H(BCV)
         ELSE
            BC_JET_G(BCV) = BC_JET_GH(BCV)
            BC_TIME(BCV) = TIME + BC_DT_H(BCV)
         ENDIF

         DO K = BC_K_B(BCV), BC_K_T(BCV)
         DO J = BC_J_S(BCV), BC_J_N(BCV)
         DO I = BC_I_W(BCV), BC_I_E(BCV)
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I,J,K)
! Why is the velocity of the boundary cell not always set (in case of
! w, s, or b plane then velocity of the adjacent fluid cell is set)?
! It should not really matter for MI...
            SELECT CASE (TRIM(BC_PLANE(BCV)))
            CASE ('W')
               IJK2 = funijk(iminus(i,j,k),j,k)
               U_G(IJK2) = BC_JET_G(BCV)
            CASE ('E')
               U_G(IJK) = BC_JET_G(BCV)
            CASE ('S')
               IJK2 =  funijk(i,jminus(i,j,k),k)
               V_G(IJK2) = BC_JET_G(BCV)
            CASE ('N')
               V_G(IJK) = BC_JET_G(BCV)
            CASE ('B')
               IJK2 =  funijk(i,j,kminus(i,j,k))
               W_G(IJK2) = BC_JET_G(BCV)
            CASE ('T')
               W_G(IJK) = BC_JET_G(BCV)
            END SELECT
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SET_BC1_JET


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc1_report_outflow                                  C
!  Purpose: print out outflow conditions                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC1_REPORT_OUTFLOW(BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_dt_0, bc_time
      use bc, only: bc_out_n
      use bc, only: bc_mout_g
      use bc, only: bc_vout_g
      use funits, only: dmp_log, unit_log
      use param1, only: undefined, zero
      use run, only: time, dt, tstop
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! index for boundary condition
      INTEGER, INTENT(IN) :: BCV

      IF (BC_DT_0(BCV) == UNDEFINED) RETURN

      CALL CALC_OUTFLOW(BCV)

! Calculate and accumulate the actual mass and volume outflow
      IF (TIME + 0.1d0*DT>=BC_TIME(BCV) .OR. &
          TIME+0.1d0*DT>=TSTOP) THEN
         BC_TIME(BCV) = TIME + BC_DT_0(BCV)

! Average and print out the flow rates
         BC_MOUT_G(BCV) = ABS(BC_MOUT_G(BCV))/BC_OUT_N(BCV)
         BC_VOUT_G(BCV) = ABS(BC_VOUT_G(BCV))/BC_OUT_N(BCV)
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) BCV, TIME
         IF(DMP_LOG)WRITE (UNIT_LOG, 1100) BC_MOUT_G(BCV), &
            BC_VOUT_G(BCV)
         BC_MOUT_G(BCV) = ZERO
         BC_VOUT_G(BCV) = ZERO
         BC_OUT_N(BCV) = 0
      ENDIF

 1000 FORMAT(/,1X,'Average outflow rates at BC No. ',I2,'  At Time = ',G12.5)
 1100 FORMAT(3X,'Gas : Mass flow = ',G12.5,'     Volumetric flow = ',G12.5)

      RETURN
      END SUBROUTINE SET_BC1_REPORT_OUTFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc1_adjust_outflow                                  C
!  Purpose: Adjust velocities to get specified mass or volumetric      C
!  flow rate based on average outflow rate                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC1_ADJUST_OUTFLOW(BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_dt_0, bc_time
      use bc, only: bc_i_w, bc_i_e
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_massflow_g
      use bc, only: bc_mout_g
      use bc, only: bc_out_n
      use bc, only: bc_plane
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use bc, only: bc_volflow_g
      use bc, only: bc_vout_g
      use compar, only: dead_cell_at
      use fldvar, only: u_g, v_g, w_g
      use functions, only: funijk, iminus, jminus, kminus
      use functions, only: is_on_mype_plus2layers
      use funits, only: dmp_log, unit_log
      use param1, only: undefined, zero, small_number
      use run, only: time, dt, tstop

      IMPLICIT NONE
! Dummy arguments
!---------------------------------------------------------------------//
! index for boundary condition
      INTEGER, INTENT(IN) :: BCV

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I, J, K, IJK
! IJK index for setting velocity bc
      INTEGER :: IJK2
!---------------------------------------------------------------------//

      CALL CALC_OUTFLOW(BCV)

! Calculate and accumulate the actual mass and volume outflow
      IF (TIME + 0.1d0*DT>=BC_TIME(BCV) .OR. &
          TIME+0.1d0*DT>=TSTOP) THEN
         BC_TIME(BCV) = TIME + BC_DT_0(BCV)

! Average and print out the flow rates
         BC_MOUT_G(BCV) = ABS(BC_MOUT_G(BCV))/BC_OUT_N(BCV)
         BC_VOUT_G(BCV) = ABS(BC_VOUT_G(BCV))/BC_OUT_N(BCV)
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) BCV, TIME
         IF(DMP_LOG)WRITE (UNIT_LOG, 1100) BC_MOUT_G(BCV), &
            BC_VOUT_G(BCV)
         BC_OUT_N(BCV) = 0

! Now that we know the mass and volume outflow update the bc velocities
! (gas phase)
         IF (BC_MASSFLOW_G(BCV) /= UNDEFINED) THEN
            IF (BC_MOUT_G(BCV) > SMALL_NUMBER) THEN
               SELECT CASE (TRIM(BC_PLANE(BCV)))
               CASE ('W', 'E')
                  BC_U_G(BCV) = BC_U_G(BCV)*BC_MASSFLOW_G(BCV)/&
                     BC_MOUT_G(BCV)
               CASE ('S', 'N')
                  BC_V_G(BCV) = BC_V_G(BCV)*BC_MASSFLOW_G(BCV)/&
                     BC_MOUT_G(BCV)
               CASE ('B', 'T')
                  BC_W_G(BCV) = BC_W_G(BCV)*BC_MASSFLOW_G(BCV)/&
                     BC_MOUT_G(BCV)
               END SELECT
            ENDIF
         ELSEIF (BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN
            IF (BC_VOUT_G(BCV) > SMALL_NUMBER) THEN
               SELECT CASE (TRIM(BC_PLANE(BCV)))
               CASE ('W', 'E')
                  BC_U_G(BCV) = BC_U_G(BCV)*BC_VOLFLOW_G(BCV)/&
                     BC_VOUT_G(BCV)
               CASE ('S', 'N')
                  BC_V_G(BCV) = BC_V_G(BCV)*BC_VOLFLOW_G(BCV)/&
                     BC_VOUT_G(BCV)
               CASE ('B', 'T')
                  BC_W_G(BCV) = BC_W_G(BCV)*BC_VOLFLOW_G(BCV)/&
                     BC_VOUT_G(BCV)
               END SELECT
            ENDIF
         ENDIF
! zero out counter for new cycle
         BC_MOUT_G(BCV) = zero
         BC_VOUT_G(BCV) = zero

! Apply updated boundary velocities - Define the field variables at the
! boundaries according to user specifications with modifications from
! the above calculations.
! If the boundary plane is W, S, or B (i.e., the fluid cell is on the
! west, south or bottom of the boundary cell) then define the velocity
! of the adjacent fluid cell according to the boundary velocity rather
! than the velocity of the boundary cell.
! Why not set the velocity in the boundary cell itself?  Based on the
! momentum bc routine it should not really matter for MO.
         DO K = BC_K_B(BCV), BC_K_T(BCV)
         DO J = BC_J_S(BCV), BC_J_N(BCV)
         DO I = BC_I_W(BCV), BC_I_E(BCV)
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I,J,K)
            SELECT CASE (TRIM(BC_PLANE(BCV)))
            CASE ('W')
               IJK2 = funijk(iminus(i,j,k),j,k)
               U_G(IJK2) = BC_U_G(BCV)
            CASE ('E')
               U_G(IJK) = BC_U_G(BCV)
            CASE ('S')
               IJK2 = funijk(i,jminus(i,j,k),k)
               V_G(IJK2) = BC_V_G(BCV)
            CASE ('N')
               V_G(IJK) = BC_V_G(BCV)
            CASE ('B')
               IJK2 = funijk(i,j,kminus(i,j,k))
               W_G(IJK2) = BC_W_G(BCV)
            CASE ('T')
               W_G(IJK) = BC_W_G(BCV)
            END SELECT

         ENDDO
         ENDDO
         ENDDO
      ENDIF   ! if time to update outflow condition

 1000 FORMAT(/,1X,'Average outflow rates at BC No. ',I2,'  At Time = ',G12.5)
 1100 FORMAT(3X,'Gas : Mass flow = ',G12.5,'     Volumetric flow = ',G12.5)

      RETURN
      END SUBROUTINE SET_BC1_ADJUST_OUTFLOW
