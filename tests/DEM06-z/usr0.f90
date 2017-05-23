!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: usr0
!                                                                      C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
!           all indices are undefined.                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine usr0

      use usr, only: rk4_pos, rk4_vel

      implicit none

      ! Store the initial position and velocity
      RK4_POS(:) = (/0.0050d0, 0.0050d0, 0.0925d0/)
      RK4_VEL(:) = (/0.0000d0, 0.0000d0, 0.0000d0/)

      end subroutine usr0
