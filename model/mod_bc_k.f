!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: MOD_BC_K (BC, I_w, J_s, K_b, K_t, PLANE)               !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: modify the "K" values for the b.c. plane                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MOD_BC_K(BCV)

      use bc, only: BC_I_W, BC_I_E
      use bc, only: BC_J_S, BC_J_N
      use bc, only: BC_K_B, BC_K_T
      use bc, only: BC_PLANE

      use geometry, only: FLAG
      use ic, only: FLUID_
      use ic, only: NSW_, FSW_, PSW_

      use compar, only: myPE

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      use open_files_mod, only: open_pe_log

      IMPLICIT NONE

! boundary condition index
      INTEGER, INTENT(in) :: BCV

! calculated cell indices in I,J,K directions
      INTEGER :: OWNER

      INTEGER :: I, J, K

      INTEGER :: IER
      LOGICAL :: ERROR
      INTEGER :: K_FLUID
      INTEGER :: K_WALL


!-----------------------------------------------

      CALL INIT_ERR_MSG("MOD_BC_K")

      K = BC_K_B(BCV)
      J = BC_J_S(BCV)
      I = BC_I_W(BCV)


! Establish the OWNER of the BC
      OWNER = myPE
      ! CALL GLOBAL_ALL_SUM(OWNER)

      IF(myPE == OWNER) THEN

         if(FLAG(i,j,k+1,1) ==FLUID_ .and. (&
            flag(i,j,k,1) == NSW_ .or. &
            flag(i,j,k,1) == FSW_ .or. &
            flag(i,j,k,1) == PSW_)) then
            BC_PLANE(BCV) = 'T'

         ELSEIF(FLAG(i,j,k,1) == FLUID_ .and. (&
            flag(i,j,k+1,1) == NSW_ .or. &
            flag(i,j,k+1,1) == FSW_ .or. &
            flag(i,j,k+1,1) == PSW_)) then
            BC_K_b(BCV) = BC_K_b(BCV) + 1
            BC_K_t(BCV) = BC_K_t(BCV) + 1
            BC_PLANE(BCV) = 'B'
         ELSE
            BC_PLANE(BCV) = '.'
         ENDIF
      ENDIF

! The owner distributes the new Iw/Ie coordinates to the other ranks.
      !CALL BCAST(K_B,OWNER)
      !CALL BCAST(K_T,OWNER)
      !CALL BCAST(BC_PLANE(BCV),OWNER)

! If there is an error, send i,j,k to all ranks. Report and exit.
      IF(BC_PLANE(BCV) == '.') THEN

         WRITE(ERR_MSG, 1100) BCV, BC_K_B(BCV), BC_K_T(BCV), &
            BC_I_W(BCV), BC_J_S(BCV)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Cannot locate flow plane for boundary ',     &
         'condition ',I3,'.',2/3x,'K Bottom =  ',I6,' K Top    = ',I6,/&
         3x,'I West   =  ',I6,' J South  = ',I6)

! Set up the I-indices for checking the entire BC region.
      K_WALL = BC_K_B(BCV)
      K_FLUID = merge(K_WALL-1, K_WALL+1, BC_PLANE(BCV)=='B')


      ERROR = .FALSE.
      DO J = BC_J_S(BCV), BC_J_N(BCV)
      DO I = BC_I_W(BCV), BC_I_E(BCV)

! Only check cells that you own and contain fluid.
         IF(FLAG(i,j,k_fluid,1) /= FLUID_ .and. (&
            flag(i,j,k_wall,1) /= NSW_ .or. &
            flag(i,j,k_wall,1) /= FSW_ .or. &
            flag(i,j,k_wall,1) /= PSW_)) ERROR = .TRUE.

      ENDDO
      ENDDO

! Sync up the error flag across all processes.
      ! CALL GLOBAL_ALL_OR(ERROR)
! If an error is detected, have each rank open a log file and write
! it's own message. Otherwise, we need to send all the data back to
! PE_IO and that's too much work!
      IF(ERROR) THEN

         CALL OPEN_PE_LOG(IER)

         WRITE(ERR_MSG, 1200) BCV
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 1200 FORMAT('Error 1200: Illegal geometry for boundary condition:',I3)

         DO J = BC_J_S(BCV), BC_J_N(BCV)
         DO I = BC_I_W(BCV), BC_I_E(BCV)

! Only check cells that you own and contain fluid.
            IF(FLAG(i,j,k_fluid,1) /= FLUID_ .and. (&
               flag(i,j,k_wall,1) /= NSW_ .or. &
               flag(i,j,k_wall,1) /= FSW_ .or. &
               flag(i,j,k_wall,1) /= PSW_)) then

               WRITE(ERR_MSG, 1201) I, J, K_WALL,  FLAG(i,j,k_wall,1),  &
                  I, J, K_FLUID, FLAG(i,j,k_fluid,1)
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF

 1201 FORMAT(' ',/14X,'I',7X,'J',7X,'K',7X,'FLAG',/3x,        &
         'WALL ',3(2x,I6),3x,I3,/3x,'FLUID',3(2x,I6),3x,I3)

         ENDDO
         ENDDO

         WRITE(ERR_MSG,"('Please correct the mfix.dat file.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)

      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MOD_BC_K
