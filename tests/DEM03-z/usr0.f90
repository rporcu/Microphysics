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

      use usr, only: f1b, f2b, gy1, gy2, gz1, gz2

      IMPLICIT NONE

      gz1 = 0.00045d0
      gz2 = 0.00135d0

      gy1 = 0.0d0
      gy2 = 0.0d0

      ! Body forces. (Gravity)
      F1b = gravity(3)
      F2b = gravity(3)

      end subroutine usr0
