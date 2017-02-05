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

      use discretelement, only: particles
      use usr, only: f1b, f2b, gx1, gx2, gy1, gy2

      IMPLICIT NONE

      gY1 = 0.00045d0
      gY2 = 0.00135d0

      gX1 = 0.0d0
      gX2 = 0.0d0

      ! Body forces -- zero in this case.
      F1b = 0.d0
      F2b = 0.d0

      end subroutine usr0
