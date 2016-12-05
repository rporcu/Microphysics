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

      SUBROUTINE CHECK_DATA_20(ep_g,p_g,ro_g,rop_g,u_g,v_g,w_g, flag)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1   , only: one, undefined, small_number
      USE compar   , only: istart2, iend2, jstart2, jend2, kstart2, kend2
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE functions, only: iminus, jminus, kminus

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,0:4)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K
! Logical variable to set, if there is an error
      LOGICAL :: ABORT
!-----------------------------------------------

      CALL INIT_ERR_MSG("CHECK_DATA_20")

      ! call send_recv(p_g,2)
      ! call send_recv(ep_g,2)
      ! call send_recv(w_g,2)
      ! call send_recv(u_g,2)
      ! call send_recv(v_g,2)
      ! call send_recv( ROP_G, 2 )
      ! call send_recv( RO_G, 2 )

      ABORT = .FALSE.

! Check whether all field variables are initialized in all fluid cells
! and flow boundary cells
! ---------------------------------------------------------------->>>
      DO K = kstart2, kend2
      DO J = jstart2, jend2
      DO I = istart2, iend2
         IF (flag(i,j,k,1)<100) THEN

! check gas phase fields
            IF(EP_G(I,J,K) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K, 'EP_G')
            IF(P_G(I,J,K) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K, 'P_G')
            IF(RO_G(I,J,K) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K, 'RO_G')
            IF(ROP_G(I,J,K) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K, 'ROP_G')

            IF(U_G(I,J,K) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K, 'U_G')


            IF(V_G(I,J,K) == UNDEFINED) then
               write(*,*) 'here'
               CALL REPORT_ERROR(ABORT, I, J, K, 'V_G')
            endif

            IF(W_G(I,J,K) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K, 'W_G')

            IF(U_G(iminus(i,j,k),j,k) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I-1, J, K, 'U_G')


            IF(V_G(i,jminus(i,j,k),k) == UNDEFINED) then
               write(*,*) 'or here'
               CALL REPORT_ERROR(ABORT, I, J-1, K, 'V_G')
            endif


            IF(W_G(i,j,kminus(i,j,k)) == UNDEFINED) &
               CALL REPORT_ERROR(ABORT, I, J, K-1, 'W_G')

         ENDIF  ! IF (flag(i,j,k,1)<100) THEN
      ENDDO  ! end do I = istart2, iend2
      ENDDO  ! end do J = jstart2, jend2
      ENDDO  ! end do K = kstart2, kend2


      ! CALL GLOBAL_ALL_OR(ABORT)
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

         IF (FLAG(i,j,k,1)==1 .OR. FLAG(i,j,k,1)==20) THEN

! Ep_g must have a value > 0 and < 1
            IF(EP_G(I,J,K) < SMALL_NUMBER .OR. EP_G(I,J,K) > ONE) &
               CALL REPORT_UNPHYSICAL(ABORT, I, J, K, 'EP_G', EP_G(I,J,K))

         ENDIF   ! IF (FLAG(i,j,k)==1 .OR. FLAG(i,j,k)==20) THEN
      ENDDO   ! I = istart2, iend2
      ENDDO   ! J = jstart2, jend2
      ENDDO   ! K = kstart2, kend2


      ! CALL GLOBAL_ALL_OR(ABORT)
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
      INTEGER, INTENT(IN   ) :: pI, pJ, pK
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
      INTEGER, INTENT(IN   ) :: pI, pJ, pK
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
