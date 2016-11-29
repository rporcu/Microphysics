!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: SET_ICBC_FLAG                                            !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE SET_ICBC_FLAG


      CALL INIT_ICBC_FLAG

      CALL SET_IC_FLAGS

      CALL SET_BC_FLAGS_WALL

      CALL SET_BC_FLAGS_FLOW

! Verify that ICBC flags are set for all fluid cells.
      CALL CHECK_ICBC_FLAG

      END SUBROUTINE SET_ICBC_FLAG



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: INIT_ICBC_FLAG                                           !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE INIT_ICBC_FLAG

      use param1, only: zero
      use ic
      use functions
      use run, only: run_type

      implicit none
      INTEGER :: I, J, K, IJK

! Initialize the icbc_flag array.
      DO k = kstart3, kend3
      DO j = jstart3, jend3
      DO i = istart3, iend3

         IJK = FUNIJK(I,J,K)

         ! Initialize the ICBC Flag
         if (run_type == 'NEW') then
            icbc_flag(i,j,k) = icbc_undef
         else
            icbc_flag(i,j,k) = icbc_fluid
         endif

! If at domain boundaries then set default values (wall or, if
! specified, cyclic)
         IF (DO_K) THEN
            IF(K==KMIN3 .OR. K==KMIN2 .OR. K==KMAX2 .OR. K==KMAX3)THEN
               IF (CYCLIC_Z_PD) THEN
                  ICBC_FLAG(i,j,k) = 1000 + icbc_cyclp
               ELSEIF (CYCLIC_Z) THEN
                  ICBC_FLAG(i,j,k) = 1000 + icbc_cycl
               ELSE
                  ICBC_FLAG(i,j,k) = 1000 + icbc_no_s
               ENDIF
            ENDIF
         ENDIF

         IF(DO_J)THEN
            IF(J==JMIN3 .OR. J==JMIN2 .OR. J==JMAX2 .OR. J==JMAX3)THEN
               IF (CYCLIC_Y_PD) THEN
                  ICBC_FLAG(i,j,k) = 1000 + icbc_cyclp
               ELSEIF (CYCLIC_Y) THEN
                  ICBC_FLAG(i,j,k) = 1000 + icbc_cycl
               ELSE
                 ICBC_FLAG(i,j,k) = 1000 + icbc_no_s
               ENDIF
            ENDIF
         ENDIF

         IF(DO_I)THEN
            IF(I==IMIN3 .OR. I==IMIN2 .OR. I==IMAX2 .OR. I==IMAX3)THEN
               IF (CYCLIC_X_PD) THEN
                  ICBC_FLAG(i,j,k) = 1000 + icbc_cyclp
               ELSEIF (CYCLIC_X) THEN
                  ICBC_FLAG(i,j,k) = 1000 + icbc_cycl
               ELSE
                  ICBC_FLAG(i,j,k) = 1000 + icbc_no_s
               ENDIF
            ENDIF
         ENDIF
! corner cells are wall cells
         IF ((I==IMIN3 .OR. I==IMIN2 .OR. I==IMAX2 .OR. I==IMAX3) .AND. &
             (J==JMIN3 .OR. J==JMIN2 .OR. J==JMAX2 .OR. J==JMIN3) .AND. &
             (K==KMIN3 .OR. K==KMIN2 .OR. K==KMAX2 .OR. K==KMAX3)) THEN
            IF (ICBC_FLAG(i,j,k) /= icbc_free) ICBC_FLAG(i,j,k) = icbc_no_s
         ENDIF

      ENDDO ! end do loop (i=istart3, iend3)
      ENDDO ! end do loop (j=jstart3, jend3)
      ENDDO ! end do loop (k=kstart3, kend3)

      RETURN

      END SUBROUTINE INIT_ICBC_FLAG



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_ICBC_FLAG                                         !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Verify that data was not given for undefined BC regions.   !
!  Note that the error message may be incomplete
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_ICBC_FLAG

      use run, only: RUN_TYPE
      use error_manager
      use ic, only: icbc_undef
      use functions

      IMPLICIT NONE

      LOGICAL :: ERROR = .FALSE.

      INTEGER :: I, J ,K, IER

      IF(RUN_TYPE(1:3) /= 'NEW') RETURN

      CALL INIT_ERR_MSG("CHECK_ICBC_FLAG")

      ! First check for any errors.
      DO K = kStart2, kEnd2
      DO J = jStart2, jEnd2
      DO I = iStart2, iEnd2
         IF (ICBC_FLAG(i,j,k) == icbc_undef) ERROR = .TRUE.
      ENDDO
      ENDDO
      ENDDO

! Sync up the error flag across all processes.
      ! CALL GLOBAL_ALL_OR(ERROR)

! If an error is detected, have each rank open a log file and write
! it's own message. Otherwise, we need to send all the data back to
! PE_IO and that's too much work!
      IF(ERROR) THEN

         CALL OPEN_PE_LOG(IER)

         WRITE(ERR_MSG, 1100) trim(iVal(myPE))
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

         DO K = kStart2, kEnd2
         DO J = jStart2, jEnd2
         DO I = iStart2, iEnd2
            IF (ICBC_FLAG(i,j,k) == icbc_undef) then
               WRITE(ERR_MSG,1101) I, J, K
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF

         ENDDO
         ENDDO
         ENDDO

         WRITE(ERR_MSG, 1102)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)

      ELSE
! If no erros, sync up the ghost cell layers.
         ! CALL SEND_RECV(ICBC_FLAG,1)
      ENDIF

! Clean up and return.
      CALL FINL_ERR_MSG

      RETURN

 1100 FORMAT('Error 1100 (PE ',A,') : No initial or boundary ',        &
         'condtions specified in','the following cells:',/             &
         '    I       J       K')

 1101 FORMAT(I5,3X,I5,3X,I5)

 1102 FORMAT('Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_ICBC_FLAG



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IC_FLAGS                                            !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Set the IC portions of the ICBC_Flag array.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_IC_FLAGS

      use ic, only: IC_DEFINED
      use ic, only: IC_TYPE

      use ic, only: IC_I_W, IC_I_E
      use ic, only: IC_J_S, IC_J_N
      use ic, only: IC_K_B, IC_K_T
      use ic, only: icbc_fluid

      use param, only: dimension_ic

      use error_manager
      use functions

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: ICV
      INTEGER :: I, J, K

      CALL INIT_ERR_MSG("SET_IC_FLAGS")


      IC_LP: DO ICV=1, DIMENSION_IC

         IF(.NOT.IC_DEFINED(ICV)) CYCLE IC_LP

! Skip checks for PATCH restarts.
         IF (IC_TYPE(ICV) == 'PATCH') CYCLE IC_LP

!  Set ICBC flag
         DO K = IC_K_B(ICV), IC_K_T(ICV)
         DO J = IC_J_S(ICV), IC_J_N(ICV)
         DO I = IC_I_W(ICV), IC_I_E(ICV)
            if  (.not.is_on_myPE_plus2layers(I,J,K)) CYCLE
            if (dead_cell_at(I,J,K)) CYCLE
            icbc_flag(i,j,k) = 1000*ICV + icbc_fluid
         ENDDO
         ENDDO
         ENDDO


      ENDDO IC_LP

! Update the ICBC flag on ghost cells.
      ! CALL SEND_RECV(ICBC_FLAG, 1)


! Clean up and return.
      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE SET_IC_FLAGS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_FLAGS_WALL                                       !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for walls BC's         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_FLAGS_WALL

      USE param
      USE param1
      USE geometry
      USE fldvar
      USE physprop
      USE bc
      USE ic
      USE funits
      USE compar
      USE functions

      use error_manager

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop/variable indices
      INTEGER :: I , J , K
! loop index
      INTEGER :: BCV

! Total number of valid BC types
      INTEGER, PARAMETER :: DIM_BCTYPE = 21

!-----------------------------------------------

      CALL INIT_ERR_MSG("SET_BC_FLAGS_WALL")

! Set the wall flags.
      DO BCV=1, DIMENSION_BC
         IF(.NOT.BC_DEFINED(BCV)) CYCLE

         IF(BC_TYPE(BCV)=='FREE_SLIP_WALL' .OR. &
            BC_TYPE(BCV)=='NO_SLIP_WALL'   .OR. &
            BC_TYPE(BCV)=='PAR_SLIP_WALL') THEN

            DO K = BC_K_B(BCV), BC_K_T(BCV)
            DO J = BC_J_S(BCV), BC_J_N(BCV)
            DO I = BC_I_W(BCV), BC_I_E(BCV)

               if (.not. is_on_myPE_plus2layers(i,j,k)) CYCLE
               if (dead_cell_at(I,J,K)) CYCLE

               SELECT CASE (TRIM(BC_TYPE(BCV)))
                  CASE('FREE_SLIP_WALL'); ICBC_FLAG(i,j,k) = 1000*BCV + icbc_free
                  CASE('NO_SLIP_WALL');   ICBC_FLAG(i,j,k) = 1000*BCV + icbc_no_s
                  CASE('PAR_SLIP_WALL');  ICBC_FLAG(i,j,k) = 1000*BCV + icbc_pslip
               END SELECT
            ENDDO
            ENDDO
            ENDDO

         ENDIF
      ENDDO

      ! CALL SEND_RECV(ICBC_FLAG,1)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_BC_FLAGS_WALL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_FLAGS_FLOW                                       !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for flow BC's. Also    !
!           set value of bc_plane for flow BC's.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_FLAGS_FLOW

      USE param
      USE param1
      USE geometry
      USE fldvar
      USE physprop
      USE bc
      USE ic
      USE funits
      USE compar
      use functions

      use error_manager

      IMPLICIT NONE

! loop/variable indices
      INTEGER :: BCV, I, J, K, IJK

      INTEGER :: IER

! error indicator
      LOGICAL :: ERROR
! surface indictors
      LOGICAL :: X_CONSTANT, Y_CONSTANT, Z_CONSTANT

      CALL INIT_ERR_MSG("SET_BC_FLAGS_FLOW")


! FIND THE FLOW SURFACES
      ERROR = .FALSE.

      DO BCV = 1, DIMENSION_BC

         IF(.NOT.BC_DEFINED(BCV)) CYCLE

         IF(BC_TYPE(BCV)=='MASS_INFLOW'  .OR. &
            BC_TYPE(BCV)=='MASS_OUTFLOW' .OR. &
            BC_TYPE(BCV)=='P_INFLOW'     .OR. &
            BC_TYPE(BCV)=='P_OUTFLOW'    .OR. &
            BC_TYPE(BCV)=='OUTFLOW') THEN

            X_CONSTANT = (BC_X_W(BCV) == BC_X_E(BCV))
            Y_CONSTANT = (BC_Y_S(BCV) == BC_Y_N(BCV))
            Z_CONSTANT = (BC_Z_B(BCV) == BC_Z_T(BCV))

            IF(X_CONSTANT .AND. BC_X_W(BCV)/=UNDEFINED)                &
               CALL MOD_BC_I(BCV)

            IF(Y_CONSTANT .AND. BC_Y_S(BCV)/=UNDEFINED)                &
               CALL MOD_BC_J(BCV)

            IF(Z_CONSTANT .AND. BC_Z_B(BCV)/=UNDEFINED)                &
               CALL MOD_BC_K(BCV)

! Extend the boundaries for cyclic implementation
            IF(BC_I_W(BCV) == 2 .AND. BC_I_E(BCV) == (IMAX2 - 1) .AND. &
               CYCLIC_X .AND. NODESI > 1) THEN
                   BC_I_W(BCV) = 1
                   BC_I_E(BCV) = IMAX2
            ENDIF
            IF(BC_J_S(BCV) == 2 .AND. BC_J_N(BCV) == (JMAX2 - 1) .AND. &
               CYCLIC_Y .AND. NODESJ > 1) THEN
               BC_J_S(BCV) = 1
               BC_J_N(BCV) = JMAX2
            ENDIF
            IF(BC_K_B(BCV) == 2 .AND. BC_K_T(BCV) == (KMAX2 - 1) .AND. &
               CYCLIC_Z .AND. NODESK > 1) THEN
               BC_K_B(BCV) = 1
               BC_K_T(BCV) = KMAX2
            ENDIF

! Set add the BC to the ICBC_FLAG. If a "non-wall" BC is found, then flag
! this as an error. The next triple-loop will take care of reporting the
! error.
            ERROR = .FALSE.
            DO K = BC_K_B(BCV), BC_K_T(BCV)
            DO J = BC_J_S(BCV), BC_J_N(BCV)
            DO I = BC_I_W(BCV), BC_I_E(BCV)

               IF(.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IF(DEAD_CELL_AT(I,J,K)) CYCLE

               ! Verify that the FLOW BC is overwriting a wall.
               IF(WALL_ICBC_FLAG(i,j,k)) THEN

                  SELECT CASE (TRIM(BC_TYPE(BCV)))
                     CASE ('P_OUTFLOW');    ICBC_FLAG(i,j,k) = 1000*BCV + icbc_p_out
                     CASE ('MASS_INFLOW');  ICBC_FLAG(i,j,k) = 1000*BCV + icbc_m_inf
                     CASE ('MASS_OUTFLOW'); ICBC_FLAG(i,j,k) = 1000*BCV + icbc_m_out
                     CASE ('OUTFLOW');      ICBC_FLAG(i,j,k) = 1000*BCV + icbc_outfl
                     CASE ('P_INFLOW');     ICBC_FLAG(i,j,k) = 1000*BCV + icbc_p_inf
                  END SELECT

               ELSE
                  ERROR = .TRUE.
               ENDIF

            ENDDO
            ENDDO
            ENDDO

! Sync the error flag over all ranks.
            ! CALL GLOBAL_ALL_OR(ERROR)

! Report errors and exit.
            IF(ERROR)THEN

               CALL OPEN_PE_LOG(IER)

               WRITE(ERR_MSG, 1200) BCV
               CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 1200 FORMAT('Error 1200: Boundary condition ',I3,' overlaps with ',&
         'another BC.',2/7x,'I',7x,'J',7x,'K',3x,'ICBC')

               DO K = BC_K_B(BCV), BC_K_T(BCV)
               DO J = BC_J_S(BCV), BC_J_N(BCV)
               DO I = BC_I_W(BCV), BC_I_E(BCV)

                  IF(.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                  IF(DEAD_CELL_AT(I,J,K)) CYCLE

                  IJK = FUNIJK(I,J,K)

! Verify that the FLOW BC is overwriting a wall.
                  IF(.NOT.WALL_ICBC_FLAG(i,j,k)) THEN
                     WRITE(ERR_MSG, 1201) I,J,K, ICBC_FLAG(i,j,k)
                     CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
                  ENDIF

 1201 FORMAT(1x,3(2x,I6),3x,A3)

               ENDDO
               ENDDO
               ENDDO

               WRITE(ERR_MSG,"('Please correct the mfix.dat file.')")
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)

            ENDIF ! IF(ERROR)
         ENDIF ! IF(not a wall BC)
      ENDDO ! BC Loop

! Sync the ICBC flag across ghost layers
      ! CALL SEND_RECV(ICBC_FLAG,1)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_BC_FLAGS_FLOW
