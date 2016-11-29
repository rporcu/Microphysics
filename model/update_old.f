!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: UPDATE_OLD                                             C
!  Purpose: Update the stored previous-time-step values of certain     C
!           field variables                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE UPDATE_OLD(old,new)

! Modules
!-----------------------------------------------
      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3

      implicit none 

      double precision, intent(in   ) :: new(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: old(istart3:iend3,jstart3:jend3,kstart3:kend3)

! Local Variables
!-----------------------------------------------
      integer :: i, j, k

      do k=kstart3, kend3
         do j=jstart3,jend3
            do i=istart3, iend3
               old(i,j,k) = new(i,j,k)
            enddo
         enddo
      enddo
      
      END SUBROUTINE UPDATE_OLD
