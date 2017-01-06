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

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      USE compar, only: iend3, jend3, kend3
      USE compar, only: mype
      USE funits, only: dmp_log
      USE param1, only: zero, one

! linear equation matrix and vector
      real(c_real), DIMENSION(:,:,:, :), ALLOCATABLE :: A_m
      real(c_real), DIMENSION(:,:,:), ALLOCATABLE :: B_m

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
         use exit_mod, only: mfix_exit
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
!  Module name: Init_Ab_m(A_m, b_m)                                    C
!  Author: M. Syamlal                                 Date: 16-MAY-96  C
!                                                                      C
!  Purpose:Initialiize the sparse matrix coefficients and the          C
!           source vector.                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE INIT_AB_M(slo, shi, A_M, B_M)

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      ! Septadiagonal matrix A_m
      real(c_real) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      ! Vector b_m
      real(c_real) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

!-----------------------------------------------
!
      A_M(:,:,:,B) = ZERO
      A_M(:,:,:,S) = ZERO
      A_M(:,:,:,W) = ZERO
      A_M(:,:,:,0) = -ONE
      A_M(:,:,:,E) = ZERO
      A_M(:,:,:,N) = ZERO
      A_M(:,:,:,T) = ZERO
      B_M = ZERO

      RETURN
      END SUBROUTINE INIT_AB_M


      END MODULE matrix
