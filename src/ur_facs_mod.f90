module ur_facs

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   use param, only: DIM_eqS

! Under relaxation factors for coefficient update:
!  [0]  every time step (explicit)
!  [1]  every iteration (implicit)
! (0,1) under-relaxed
   real(c_real) :: UR_FAC(DIM_eqS)

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: UNDER_RELAX                                             !
!  Author: M. Syamlal                                 Date: 24-MAY-96  !
!                                                                      !
!  Purpose: Under-relax equation.                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine under_relax(var, varlo, varhi, A_m, b_m, alo, ahi, eq)

   use param, only: one, equal

   implicit none

   integer(c_int), intent(in   ) :: varlo(3),varhi(3)
   integer(c_int), intent(in   ) ::   alo(3),  ahi(3)
   real(c_real)  , intent(in   ) :: var(varlo(1):varhi(1),varlo(2):varhi(2),varlo(3):varhi(3))
   real(c_real)  , intent(inout) :: A_m(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
   real(c_real)  , intent(inout) :: b_m(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
   integer(c_int), intent(in   ) :: eq

   ! Loop index
   integer :: i, j, k

   ! Functions of under-relaxation factor
   real(c_real) :: f1, f2

   ! Center coefficient
   real(c_real) :: Ap

   F1 = ONE/UR_FAC(eq)
   F2 = F1 - ONE

   do k = alo(3),ahi(3)
      do j = alo(2),ahi(2)
         do i = alo(1),ahi(1)
            Ap = A_m(I,J,K,0)
            A_m(I,J,K,0) = Ap*F1
            b_m(I,J,K) = b_m(I,J,K) + Ap*VAR(i,j,k)*F2
         end do
      end do
   end do

   end subroutine under_relax

end module ur_facs
