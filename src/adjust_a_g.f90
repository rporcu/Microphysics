   module adjust_a

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: adjust_a(A_m, B_m, IER)                            C
!  Author: M. Syamlal                                 Date:  2-AUG-96  C
!                                                                      C
!  Purpose: Handle the special case of the center coefficient in       C
!  U_g momentum eq. becoming zero.                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

contains

   subroutine adjust_a_g(axis, slo, shi, alo, ahi, lo, hi, A_m, b_m, rop_g, dx, dy, dz)

      use functions, only: avg
      use param, only: ONE, ZERO, small_number

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3),alo(3),ahi(3),lo(3),hi(3)
      character     , intent(in   ) :: axis

      ! Septadiagonal matrix A_m
      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3)

      ! Vector b_m
      real(c_real), intent(inout) :: B_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      real(c_real), intent(in   ) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: dx, dy, dz
!---------------------------------------------------------------------//

      integer      :: i, j, k, is, js, ks
      real(c_real) :: denominator, xxxm, xxxp
      real(c_real) :: axy, ayz, axz

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

      is = lo(1)
      js = lo(2)
      ks = lo(3)

!     if (axis.eq.'x') then
!         is = lo(1)-1
!     else if (axis.eq.'y') then
!         js = lo(2)-1
!     else if (axis.eq.'z') then
!         ks = lo(3)-1
!     endif

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

         if (abs(A_m(i,j,k,0)) < SMALL_NUMBER) THEN

            A_m(i,j,k,:) = ZERO
            A_m(i,j,k,0) = -ONE

            if (b_m(i,j,k) < ZERO) THEN

               if (axis .eq. 'U') then
                  denominator = rop_g(i+1,j,k)*AYZ
               else if (axis .eq. 'V') then
                  denominator = rop_g(i,j+1,k)*AXZ
               else if (axis .eq. 'W') then
                  denominator = rop_g(i,j,k+1)*AXY
               end if

               xxxm = ONE
               xxxp = ZERO

            else if (b_m(i,j,k) > ZERO) THEN

               if (axis .eq. 'U') then
                  denominator = rop_g(i,j,k)*AYZ
               else if (axis .eq. 'V') then
                  denominator = rop_g(i,j,k)*AXZ
               else if (axis .eq. 'W') then
                  denominator = rop_g(i,j,k)*AXY
               end if

               xxxm = ZERO
               xxxp = ONE

            else
               denominator = ZERO
            end if

            if (denominator > SMALL_NUMBER) THEN
               b_m(i,j,k) = SQRT(ABS(b_m(i,j,k))/(denominator*AVG(xxxm,xxxp)))
            else
               b_m(i,j,k) = ZERO
            end if
         end if
      end do
      end do
      end do

      end subroutine adjust_a_g

   end module adjust_a
