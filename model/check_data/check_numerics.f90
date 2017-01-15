module check_numerics_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

contains

!  Subroutine: CHECK_NUMERICS                                          !
!  Purpose: Check the numerics control namelist section                !
!                                                                      !
!  Author: P. Nicoletti                               Date: 27-NOV-91  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine CHECK_NUMERICS


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
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival

      IMPLICIT NONE
      INTEGER :: L

      DO L = 1,DIM_EQS
         IF(DISCRETIZE(L) > 9 .OR. DISCRETIZE(L) < 0) THEN
            WRITE(ERR_MSG,2002) trim(ivar('DISCRETIZE',L)),&
               trim(ival(DISCRETIZE(L)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

 2002 FORMAT('Error 2002: Invalid option ', A,' = ', A, '.',/  &
         'Please correct the mfix.dat file.')

      end subroutine check_numerics
end module check_numerics_module
