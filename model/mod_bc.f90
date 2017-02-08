MODULE MOD_BC
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: MOD_BC_I                                                C
!  Author: P. Nicoletti                               Date: 10-DEC-91  C
!                                                                      C
!  Purpose: modify the "I" values for the b.c. plane                   C
!     This is a yz plane                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE MOD_BC_I(BCV, flag, slo, shi)

      use bc, only: BC_I_W, BC_I_E
      use bc, only: BC_J_S, BC_J_N
      use bc, only: BC_K_B, BC_K_T
      use bc, only: BC_PLANE

      use ic, only: FLUID_
      use ic, only: NSW_, FSW_, PSW_

      use compar       , only: myPE
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      use open_files_mod, only: open_pe_log

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      ! boundary condition index
      integer, intent(in   ) :: BCV

      integer, intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      ! i cell indices defining location of yz plane
      INTEGER :: i,j,k

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: OWNER

      INTEGER :: IER
      LOGICAL :: ERROR
      INTEGER :: I_FLUID
      INTEGER :: I_WALL

      CALL INIT_ERR_MSG("MOD_BC_I")

      K = BC_K_B(BCV)
      J = BC_J_S(BCV)
      I = BC_I_W(BCV)

! Establish the OWNER of the BC
      OWNER = myPE
      ! CALL GLOBAL_ALL_SUM(OWNER)

      IF(myPE == OWNER) THEN

! Flow on west boundary (fluid cell on east).
         if(flag(i+1,j,k,1) == FLUID_ .and. (&
            flag(i,j,k,1) == NSW_ .or. &
            flag(i,j,k,1) == FSW_ .or. &
            flag(i,j,k,1) == PSW_)) then

            BC_PLANE(BCV) = 'E'

! Flow on east boundary (fluid cell on west).
         elseif(flag(i,j,k,1) == FLUID_ .and. (&
            flag(i+1,j,k,1) == NSW_ .or. &
            flag(i+1,j,k,1) == FSW_ .or. &
            flag(i+1,j,k,1) == PSW_)) then

            BC_I_W(BCV) = BC_I_W(BCV) + 1
            BC_I_E(BCV) = BC_I_E(BCV) + 1
            BC_PLANE(BCV) = 'W'

! Set the plane of a value we know to be wrong so we can detect the error.
         ELSE
            BC_PLANE(BCV) = '.'
         ENDIF
      ENDIF

      !CALL BCAST(I_W,   OWNER)
      !CALL BCAST(I_E,   OWNER)
      !CALL BCAST(BC_PLANE(BCV), OWNER)

! If there is an error, send I,J,K to all ranks. Report and exit.
      IF(BC_PLANE(BCV) == '.') THEN
         WRITE(ERR_MSG, 1100) BCV, BC_I_W(BCV), BC_I_E(BCV), &
            BC_J_S(BCV), BC_K_B(BCV)
!           FLAG(i,j,k,1), FLAG(i+1,j,k,1)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Cannot locate flow plane for boundary ',     &
         'condition ',I3,'.',2/3x,'I West   =  ',I6,' I East   = ',I6,/&
         3x,'J South  =  ',I6,' K Bottom = ',I6)

! Set up the I-indices for checking the entire BC region.
      I_WALL = BC_I_W(BCV)
      I_FLUID = merge(I_WALL-1, I_WALL+1, BC_PLANE(BCV)=='W')


! First pass through all of the BC region and verify that you have
! inflow/outflow specified against a wall. Flag any errors.
      ERROR = .FALSE.
      DO K = BC_K_B(BCV), BC_K_T(BCV)
         DO J = BC_J_S(BCV), BC_J_N(BCV)

            ! Verify that the the fluid and wall cells match the FLAG.
            ! Only check cells that you own and contain fluid.
            IF(FLAG(i_fluid,j,k,1) /= FLUID_ .and. (&
               flag(i_wall,j,k,1) /= NSW_ .or. &
               flag(i_wall,j,k,1) /= FSW_ .or. &
               flag(i_wall,j,k,1) /= PSW_)) ERROR = .TRUE.

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

         DO K = BC_K_B(BCV), BC_K_T(BCV)
            DO J = BC_J_S(BCV), BC_J_N(BCV)

               IF(FLAG(i_fluid,j,k,1) /= FLUID_ .and. (&
                  flag(i_wall,j,k,1) /= NSW_ .or. &
                  flag(i_wall,j,k,1) /= FSW_ .or. &
                  flag(i_wall,j,k,1) /= PSW_)) then

                  WRITE(ERR_MSG, 1201) I_WALL, J, K, FLAG(i_wall,j,k,1),  &
                     I_FLUID, J, K, FLAG(i_fluid,j,k,1)
                  CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
               ENDIF
            ENDDO
         ENDDO

         WRITE(ERR_MSG,"('Please correct the mfix.dat file.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)

      ENDIF


 1200 FORMAT('Error 1200: Illegal geometry for boundary condition:',I3)

 1201 FORMAT(' ',/14X,'I',7X,'J',7X,'K',7X,'IJK',4x,'FLAG',/3x,        &
         'WALL ',3(2x,I6),3x,I3,/3x,'FLUID',3(2x,I6),3x,I3)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MOD_BC_I
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: MOD_BC_J(BC, I_w, J_s, J_n, K_b, PLANE)                !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: modify the "J" values for the b.c. plane                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MOD_BC_J(BCV, flag, slo, shi)

      use bc, only: BC_I_W, BC_I_E
      use bc, only: BC_J_S, BC_J_N
      use bc, only: BC_K_B, BC_K_T
      use bc, only: BC_PLANE

      use ic, only: FLUID_
      use ic, only: NSW_, FSW_, PSW_

      use compar, only: myPE

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      use open_files_mod, only: open_pe_log

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      ! boundary condition index
      INTEGER, intent(in) :: BCV

      integer, intent(in) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      ! Calculated cell indices in I,J,K directions
      INTEGER :: OWNER
      INTEGER :: I, J, K

      INTEGER :: IER
      LOGICAL :: ERROR
      INTEGER :: J_FLUID
      INTEGER :: J_WALL


!-----------------------------------------------

      CALL INIT_ERR_MSG("MOD_BC_J")

      K = BC_K_B(BCV)
      J = BC_J_S(BCV)
      I = BC_I_W(BCV)

! Establish the OWNER of the BC
      OWNER = myPE
      ! CALL GLOBAL_ALL_SUM(OWNER)

      IF(myPE == OWNER)THEN

         if(flag(i,j+1,k,1) == FLUID_ .and. (&
            flag(i,j  ,k,1) == NSW_ .or. &
            flag(i,j  ,k,1) == FSW_ .or. &
            flag(i,j  ,k,1) == PSW_)) then

            BC_PLANE(BCV) = 'N'

         elseif(flag(i,j,k,1) == FLUID_ .and. (&
            flag(i,j+1,k,1) == NSW_ .or. &
            flag(i,j+1,k,1) == FSW_ .or. &
            flag(i,j+1,k,1) == PSW_)) then

            BC_J_S(BCV) = BC_J_S(BCV) + 1
            BC_J_N(BCV) = BC_J_N(BCV) + 1
            BC_PLANE(BCV) = 'S'

         ELSE
            BC_PLANE(BCV) = '.'
         ENDIF
      ENDIF

      !CALL BCAST(J_S,OWNER)
      !CALL BCAST(J_N,OWNER)
      !CALL BCAST(BC_PLANE(BCV),OWNER)

! If there is an error, send i,j,k to all ranks. Report and exit.
      IF(BC_PLANE(BCV) == '.') THEN
         WRITE(ERR_MSG, 1100) BCV, BC_J_S(BCV), BC_J_N(BCV), &
            BC_I_W(BCV), BC_K_B(BCV)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


 1100 FORMAT('Error 1100: Cannot locate flow plane for boundary ',     &
         'condition ',I3,'.',2/3x,'J South   =  ',I6,' J North   = ',I6,/&
         3x,'I West  =  ',I6,' K Bottom = ',I6)


      J_WALL = BC_J_S(BCV)
      J_FLUID = merge(J_WALL-1, J_WALL+1, BC_PLANE(BCV)=='S')


! First pass through all of the BC region and verify that you have
! inflow/outflow specified against a wall. Flag any errors.
      ERROR = .FALSE.
      DO K = BC_K_B(BCV), BC_K_T(BCV)
         DO I = BC_I_W(BCV), BC_I_E(BCV)
            IF(flag(i,j_fluid,k,1) /= FLUID_ .and. (&
               flag(i,j_wall,k,1) /= NSW_ .or. &
               flag(i,j_wall,k,1) /= FSW_ .or. &
               flag(i,j_wall,k,1) /= PSW_)) ERROR = .TRUE.
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

         DO K = BC_K_B(BCV), BC_K_T(BCV)
            DO I = BC_I_W(BCV), BC_I_E(BCV)

               IF(FLAG(i,j_fluid,k,1) /= FLUID_ .and. (&
                  flag(i,j_wall,k,1) /= NSW_ .or. &
                  flag(i,j_wall,k,1) /= FSW_ .or. &
                  flag(i,j_wall,k,1) /= PSW_)) then

                  WRITE(ERR_MSG, 1201) &
                     I, J_WALL,  K, FLAG(i,j_wall, k,1), &
                     I, J_FLUID, K, FLAG(i,j_fluid,k,1)
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
      END SUBROUTINE MOD_BC_J
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: MOD_BC_K (BC, I_w, J_s, K_b, K_t, PLANE)               !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: modify the "K" values for the b.c. plane                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MOD_BC_K(BCV,flag,slo,shi)

      use bc, only: BC_I_W, BC_I_E
      use bc, only: BC_J_S, BC_J_N
      use bc, only: BC_K_B, BC_K_T
      use bc, only: BC_PLANE

      use ic, only: FLUID_
      use ic, only: NSW_, FSW_, PSW_

      use compar, only: myPE

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      use open_files_mod, only: open_pe_log

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

! boundary condition index
      INTEGER, intent(in) :: BCV
      integer, intent(in) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

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
END MODULE MOD_BC
