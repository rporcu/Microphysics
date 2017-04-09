module check_plane_module

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: check_plane                                            C
!  Purpose: make sure the flow boundary condition or internal surface  C
!           is a plane                                                 C
!                                                                      C
!  Author: P. Nicoletti                               Date: 10-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      subroutine check_plane(X_CONSTANT, Y_CONSTANT, Z_CONSTANT, BC, NAME)

      use compar, only: mype
      use exit_mod, only: mfix_exit

      IMPLICIT NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------

      ! surface indicators
      logical, intent(IN) :: X_CONSTANT,Y_CONSTANT,Z_CONSTANT

      ! boundary condition or internal surface index
      integer, intent(IN) ::  BC

      ! BC or IS
      CHARACTER(LEN=2), intent(IN) :: NAME
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! number of directions that are not constant (must equal 2)
      integer :: N
!-----------------------------------------------


! number of directions that are not constant (must equal 2)
      N = 3
      IF (X_CONSTANT) N = N - 1
      IF (Y_CONSTANT) N = N - 1
      IF (Z_CONSTANT) N = N - 1

      IF (N /= 2) THEN
         WRITE (*, 1000) NAME, BC
         call mfix_exit(myPE)
      ENDIF

      RETURN

 1000 FORMAT(/70('*')//' From: check_plane',/'Message: ',A,' No ',I3,&
         ' is not a plane',/70('*')/)
      end subroutine check_plane

end module check_plane_module
