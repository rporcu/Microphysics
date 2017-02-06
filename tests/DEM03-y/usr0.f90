!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR0                                                   !
!  Author: J.Musser                                   Date: dd-mmm-yy  !
!  Purpose: This routine is called before the time loop starts and is  !
!           user-definable.  The user may insert code in this routine  !
!           or call appropriate user defined subroutines.  This        !
!           can be used for setting constants and checking errors in   !
!           data.  This routine is not called from an IJK loop, hence  !
!           all indices are undefined.                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine usr0

      use constant, only: gravity
      use usr, only: f1b, f2b, gx1, gx2, gy1, gy2

      IMPLICIT NONE

      gy1 = 0.00045d0
      gy2 = 0.00135d0
      gx1 = 0.0d0
      gx2 = 0.0d0

      ! Body forces. (Gravity)
      F1b = gravity(2)
      F2b = gravity(2)

      end subroutine usr0
