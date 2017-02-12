module ur_facs

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
   subroutine under_relax(var, varlo, varhi, flag, slo, shi, A_m, b_m, alo, ahi, AXIS, EQ)

    use param1, only: one, equal

    implicit none

    integer     , intent(in   ) :: varlo(3),varhi(3)
    integer     , intent(in   ) ::   slo(3),  shi(3)
    integer     , intent(in   ) ::   alo(3),  ahi(3)

   ! Variable
   real(c_real) :: var(varlo(1):varhi(1),varlo(2):varhi(2),varlo(3):varhi(3))

   ! Septadiagonal matrix
   real(c_real) :: A_m&
      (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

   !   Vector b_m
   real(c_real) :: B_m&
      (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

   integer, intent(in   ) ::  flag&
      (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

   ! Equation ID
   INTEGER :: EQ

   ! Axis ID: U, V, W, S (scalar)
   CHARACTER :: AXIS

   ! Loop index
   integer :: i, j, k
   ! Functions of under-relaxation factor
   real(c_real) :: f1, f2

   ! Center coefficient
   real(c_real) :: Ap

   F1 = ONE/UR_FAC(EQ)
   F2 = F1 - ONE

      if (axis.eq.'S') then

         do k = slo(3),shi(3)
            do j = slo(2),shi(2)
               do i = slo(1),shi(1)
                  IF(flag(i,j,k,1) == 1) THEN
                     AP = A_m(I,J,K,0)
                     IF (.NOT.EQUAL(AP, (-ONE))) THEN
                        A_m(I,J,K,0) = AP*F1
                        b_m(I,J,K) = b_m(I,J,K) + AP*VAR(i,j,k)*F2
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
                     AP = A_m(I,J,K,0)
                     IF (.NOT.EQUAL(AP, (-ONE))) THEN
                        A_m(I,J,K,0) = AP*F1
                        b_m(I,J,K) = b_m(I,J,K) + AP*VAR(i,j,k)*F2
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
                     AP = A_m(I,J,K,0)
                     IF (.NOT.EQUAL(AP, (-ONE))) THEN
                        A_m(I,J,K,0) = AP*F1
                        b_m(I,J,K) = b_m(I,J,K) + AP*VAR(i,j,k)*F2
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
                     AP = A_m(I,J,K,0)
                     IF (.NOT.EQUAL(AP, (-ONE))) THEN
                        A_m(I,J,K,0) = AP*F1
                        b_m(I,J,K) = b_m(I,J,K) + AP*VAR(i,j,k)*F2
                     ENDIF
                  ENDIF
               end do
            end do
         end do
      endif

      RETURN

   contains

      include 'functions.inc'

   end subroutine under_relax

end module ur_facs
