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
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR0

      use usr

      IMPLICIT NONE

! Initialize:
! Previoius position and velocity.
      yPOSO = 0.5d0
      yVELO = 0.0d0
! Bounce counter
      BOUNCE_COUNT = 0
! Max bounce height
      MAX_HEIGHT(0) = yPOSO
! Drop height
      h0 = yPOSO

      return
      END SUBROUTINE USR0
