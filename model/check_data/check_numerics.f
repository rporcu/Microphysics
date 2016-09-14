!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_NUMERICS                                          !
!  Purpose: Check the numerics control namelist section                !
!                                                                      !
!  Author: P. Nicoletti                               Date: 27-NOV-91  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_NUMERICS


! Global Variables:
!---------------------------------------------------------------------//
! Discretization scheme for various equations
      USE run, only: DISCRETIZE

      use param, only: dim_eqs

! Global Parameters:
!---------------------------------------------------------------------//
! NONE

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: L


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_NUMERICS")

      DO L = 1,DIM_EQS
         IF(DISCRETIZE(L) > 9 .OR. DISCRETIZE(L) < 0) THEN
            WRITE(ERR_MSG,2002) trim(ivar('DISCRETIZE',L)),&
               trim(ival(DISCRETIZE(L)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO
 2002 FORMAT('Error 2002: Invalid option ', A,' = ', A, '.',/  &
         'Please correct the mfix.dat file.')


! Finalize the error msg.
      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_NUMERICS
