MODULE ur_facs

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   use param, only: DIM_EQS

! Under relaxation factors for coefficient update:
!  [0]  every time step (explicit)
!  [1]  every iteration (implicit)
! (0,1) under-relaxed
   real(c_real) :: UR_FAC(DIM_EQS)

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: UNDER_RELAX                                             !
!  Author: M. Syamlal                                 Date: 24-MAY-96  !
!                                                                      !
!  Purpose: Under-relax equation.                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE UNDER_RELAX(slo, shi, VAR, A_M, B_M, AXIS, flag, EQ)

   use param1, only: one

   implicit none

   integer     , intent(in   ) :: slo(3),shi(3)

! Dummy arguments:
!---------------------------------------------------------------------//
! Variable
      real(c_real) :: Var(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
! Septadiagonal matrix
      real(c_real) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
!   Vector b_m
      real(c_real) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer, intent(in   ) ::  flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
! Equation ID
      INTEGER :: EQ
! Axis ID: U, V, W, S (scalar)
      CHARACTER :: AXIS

! Local variables:
!---------------------------------------------------------------------//
! Loop index
      INTEGER :: i, j, k
! Functions of under-relaxation factor
      real(c_real) :: f1, f2
! Center coefficient
      real(c_real) :: Ap

!.......................................................................!
      F1 = ONE/UR_FAC(EQ)
      F2 = F1 - ONE

      if (axis.eq.'S') then

         do k = slo(3),shi(3)
            do j = slo(2),shi(2)
               do i = slo(1),shi(1)
                  IF(flag(i,j,k,1) == 1) THEN
                     AP = A_M(I,J,K,0)
                     IF (AP /= (-ONE)) THEN
                        A_M(I,J,K,0) = AP*F1
                        B_M(I,J,K) = B_M(I,J,K) + AP*VAR(i,j,k)*F2
                     ENDIF
                  ENDIF
               end do
            end do
         end do

      else if (axis.eq.'U') then
         do k = slo(3),shi(3)
            do j = slo(2),shi(2)
               do i = slo(1),shi(1)
                  IF(flag(i,j,k,2) >= 2000 .and. &
                     flag(i,j,k,2) <= 2011) THEN
                     AP = A_M(I,J,K,0)
                     IF (AP /= (-ONE)) THEN
                        A_M(I,J,K,0) = AP*F1
                        B_M(I,J,K) = B_M(I,J,K) + AP*VAR(i,j,k)*F2
                     ENDIF
                  ENDIF
               end do
            end do
         end do

      else if (axis.eq.'V') then
         do k = slo(3),shi(3)
            do j = slo(2),shi(2)
               do i = slo(1),shi(1)
                  IF(flag(i,j,k,3) >= 2000 .and. &
                     flag(i,j,k,3) <= 2011) THEN
                     AP = A_M(I,J,K,0)
                     IF (AP /= (-ONE)) THEN
                        A_M(I,J,K,0) = AP*F1
                        B_M(I,J,K) = B_M(I,J,K) + AP*VAR(i,j,k)*F2
                     ENDIF
                  ENDIF
               end do
            end do
         end do

      else if (axis.eq.'W') then
         do k = slo(3),shi(3)
            do j = slo(2),shi(2)
               do i = slo(1),shi(1)
                  IF(flag(i,j,k,4) >= 2000 .and. &
                     flag(i,j,k,4) <= 2011) THEN
                     AP = A_M(I,J,K,0)
                     IF (AP /= (-ONE)) THEN
                        A_M(I,J,K,0) = AP*F1
                        B_M(I,J,K) = B_M(I,J,K) + AP*VAR(i,j,k)*F2
                     ENDIF
                  ENDIF
               end do
            end do
         end do
      endif

      RETURN

   contains

      include 'functions.inc'

   END SUBROUTINE UNDER_RELAX

END MODULE ur_facs
