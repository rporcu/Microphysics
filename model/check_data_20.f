!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_DATA_20                                           C
!  Purpose:                                                            C
!     - check whether field variables are initialized in all cells     C
!     - check whether the sum of void and volume fractions is 1.0      C
!       in all fluid and mass inflow cells                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_DATA_20

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE toleranc
      USE fldvar
      USE run
      USE geometry
      USE constant
      USE physprop
      USE indices
      USE funits
      USE compar
      USE sendrecv
      USE discretelement
      USE mfix_pic
      USE functions

      use mpi_utility
      use error_manager

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IMJK, IJMK, IJKM
! Solids phase
      INTEGER :: M
! Species index
      INTEGER :: N
! Logical variable to set, if there is an error
      LOGICAL :: ABORT
! 1.0 - sum of all volume fractions
      DOUBLE PRECISION DIF
!-----------------------------------------------

      CALL INIT_ERR_MSG("CHECK_DATA_20")

      call send_recv(p_g,2)
      call send_recv(ep_g,2)
      call send_recv(w_g,2)
      call send_recv(u_g,2)
      call send_recv(v_g,2)
      call send_recv( ROP_G, 2 )
      call send_recv( RO_G, 2 )

      ABORT = .FALSE.

! Check whether all field variables are initialized in all fluid cells
! and flow boundary cells
! ---------------------------------------------------------------->>>
      DO K = kstart2, kend2
      DO J = jstart2, jend2
      DO I = istart2, iend2
         IJK = FUNIJK(I,J,K)
         IF (.NOT.WALL_AT(IJK)) THEN

            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

! check gas phase fields
            IF(EP_G(IJK) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K, 'EP_G')
            IF(P_G(IJK) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K, 'P_G')
            IF(RO_G(IJK) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K, 'RO_G')
            IF(ROP_G(IJK) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K, 'ROP_G')

            IF(U_G(IJK) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K, 'U_G')


            IF(V_G(IJK) == UNDEFINED) then
               write(*,*) 'here'
               CALL REPORT_ERROR(ABORT, I, J, K, 'V_G')
            endif

            IF(W_G(IJK) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K, 'W_G')

            IF(U_G(IMJK) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I-1, J, K, 'U_G')


            IF(V_G(IJMK) == UNDEFINED) then
               write(*,*) 'or here'
               CALL REPORT_ERROR(ABORT, I, J-1, K, 'V_G')
            endif


            IF(W_G(IJKM) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K-1, 'W_G')

         ENDIF  ! IF (.NOT.WALL_AT(IJK)) THEN
      ENDDO  ! end do I = istart2, iend2
      ENDDO  ! end do J = jstart2, jend2
      ENDDO  ! end do K = kstart2, kend2


      CALL GLOBAL_ALL_OR(ABORT)
      IF(ABORT) THEN
         WRITE(ERR_MSG, 2000)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)
         CALL FINL_ERR_MSG
         RETURN
      ENDIF


! Additional check for fluid or mass inflow cells
! --------------------------------------------------------------------//

      DO K = kstart2, kend2
      DO J = jstart2, jend2
      DO I = istart2, iend2
         IJK = FUNIJK(I,J,K)

         IF (FLAG(IJK)==1 .OR. FLAG(IJK)==20) THEN

! Ep_g must have a value > 0 and < 1
            IF(EP_G(IJK) < SMALL_NUMBER .OR. EP_G(IJK) > ONE) &
               CALL REPORT_UNPHYSICAL(ABORT, I, J, K, 'EP_G', EP_G(IJK))

         ENDIF   ! IF (FLAG(IJK)==1 .OR. FLAG(IJK)==20) THEN
      ENDDO   ! I = istart2, iend2
      ENDDO   ! J = jstart2, jend2
      ENDDO   ! K = kstart2, kend2


      CALL GLOBAL_ALL_OR(ABORT)
      IF(ABORT) THEN
         WRITE(ERR_MSG, 2000)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)
         CALL FINL_ERR_MSG
         RETURN
      ENDIF


      CALL FINL_ERR_MSG
      RETURN

 2000 FORMAT('Please correct the mfix.dat file.')



      CONTAINS


      SUBROUTINE REPORT_ERROR(ABORT, pI, pJ, pK, VAR, LC1, LC2)

      LOGICAL, INTENT(INOUT) :: ABORT
      INTEGER, INTENT(IN) :: pI, pJ, pK
      CHARACTER(LEN=*), INTENT(IN) :: VAR
      INTEGER, INTENT(IN), OPTIONAL :: LC1, LC2
      CHARACTER(LEN=32) :: VAR_FULL

      VAR_FULL=''
      IF(PRESENT(LC2)) THEN
         VAR_FULL = iVAR(VAR,LC1,LC2)
      ELSEIF(PRESENT(LC1)) THEN
         VAR_FULL = iVAR(VAR,LC1)
      ELSE
         VAR_FULL = VAR
      ENDIF

      IF(.NOT.ABORT) THEN
         WRITE(ERR_MSG,1000)
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
         ABORT = .TRUE.
      ENDIF

 1000 FORMAT('Error 1000: The following field variables are undefined')


      WRITE(ERR_MSG, 1010) I, J, K, trim(VAR_FULL)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1010 FORMAT(1X,'I = ',I6,' J = ',I6,' K = ',I6,5X,A)



      END SUBROUTINE REPORT_ERROR

      SUBROUTINE REPORT_UNPHYSICAL(ABORT, pI, pJ, pK, VAR, VALUE)

      LOGICAL, INTENT(INOUT) :: ABORT
      INTEGER, INTENT(IN) :: pI, pJ, pK
      CHARACTER(LEN=*), INTENT(IN) :: VAR
      DOUBLE PRECISION, INTENT(IN) :: VALUE

      IF(.NOT.ABORT) THEN
         WRITE(ERR_MSG,1100)
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
         ABORT = .TRUE.
      ENDIF

 1100 FORMAT('Error 1100: The following field variables are ',&
         'out of range')

      WRITE(ERR_MSG, 1110) I, J, K, trim(VAR), VALUE
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1110 FORMAT(1X,'I = ',I6,' J = ',I6,' K = ',I6,2X,A,'Value:',g11.4)


      END SUBROUTINE REPORT_UNPHYSICAL


      END SUBROUTINE CHECK_DATA_20
