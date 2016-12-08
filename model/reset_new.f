!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RESET_NEW                                              C
!  Purpose: Reset the new variables with the stored previous-time-step C
!           values of field variables.                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine reset_new(new,old)

! Modules
!-----------------------------------------------
      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3

      implicit none

      double precision, intent(inout) :: new(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) :: old(istart3:iend3,jstart3:jend3,kstart3:kend3)

! Local Variables
!-----------------------------------------------
      integer :: i,j,k

      do k=kstart3, kend3
         do j=jstart3,jend3
            do i=istart3, iend3
               new(i,j,k) = old(i,j,k)
            enddo
         enddo
      enddo

      end subroutine reset_new
