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
      MODULE matrix

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


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Init_Ab_m(A_m, b_m, IJKMAX2, M, IER)                   C                     C
!  Author: M. Syamlal                                 Date: 16-MAY-96  C
!                                                                      C
!  Purpose:Initialiize the sparse matrix coefficients and the          C
!           source vector.                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE INIT_AB_M(A_M, B_M, IJKMAX2A, M)

      USE param
      USE param1
      USE compar

      IMPLICIT NONE
!
!                      Phase index
      INTEGER          M
!
!                      Maximum dimension
      INTEGER          IJKMAX2A
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Source vector
      DOUBLE PRECISION b_m(DIMENSION_3, 0:DIMENSION_M)
!
!-----------------------------------------------
!
      A_M(:,B,M) = ZERO
      A_M(:,S,M) = ZERO
      A_M(:,W,M) = ZERO
      A_M(:,0,M) = -ONE
      A_M(:,E,M) = ZERO
      A_M(:,N,M) = ZERO
      A_M(:,T,M) = ZERO
      B_M(:,M) = ZERO

      RETURN
      END SUBROUTINE INIT_AB_M


      END MODULE matrix
