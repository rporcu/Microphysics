!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: ambm                                                        !
!  Purpose:                                                            !
!     IMPORTANT:  For using these arrays in a subroutine               !
!     -lock the module in the beginning of the subroutine              !
!      call lock_ambm                                                  !
!     -and unlock the module at the end of the subroutine              !
!      call unlock_ambm                                                !
! Contains the following subroutines:                                  !
!      lock_ambm, unlock_ambm                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE ambm

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE funits
!-----------------------------------------------

! linear equation matrix and vector
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: A_m
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: B_m

      LOGICAL :: ambm_locked = .false.


! Definitions for sparse matrix
      INTEGER, PARAMETER :: e = 1
      INTEGER, PARAMETER :: w =-1
      INTEGER, PARAMETER :: n = 2
      INTEGER, PARAMETER :: s =-2
      INTEGER, PARAMETER :: t = 3
      INTEGER, PARAMETER :: b =-3


      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE lock_ambm
      IF(ambm_locked) THEN
         IF (DMP_LOG) WRITE(*,*) &
            'Error:  Multiple use of ambm (ambm_mod.f)'
         CALL MFIX_EXIT(myPE)
      ELSE
         ambm_locked = .true.
      ENDIF
      END SUBROUTINE lock_ambm


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE unlock_ambm
      ambm_locked = .false.
      END SUBROUTINE unlock_ambm

      END MODULE ambm


