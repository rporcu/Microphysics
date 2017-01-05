MODULE CHECK_GEOMETRY_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   USE check_axis_module, only: check_axis

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_GEOMETRY                                          !
!  Purpose: Check the distributed parallel namelist variables.         !
!                                                                      !
!  Author: P. Nicoletti                               Date: 14-DEC-99  !
!  Reviewer: J.Musser                                 Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_GEOMETRY(dx,dy,dz)


! Global Variables:
!---------------------------------------------------------------------//
! Domain partitions in various directions.
      use geometry, only: XLENGTH
      use geometry, only: YLENGTH
      use geometry, only: ZLENGTH

      use geometry, only: IMAX, IMAX3
      use geometry, only: JMAX, JMAX3
      use geometry, only: KMAX, KMAX3

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival

      implicit none

      real(c_real), intent(inout) :: dx,dy,dz

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_GEOMETRY")

      CALL CHECK_AXIS(IMAX, IMAX3, XLENGTH, DX, 'X', 'I')
      CALL CHECK_AXIS(JMAX, JMAX3, YLENGTH, DY, 'Y', 'J')
      CALL CHECK_AXIS(KMAX, KMAX3, ZLENGTH, DZ, 'Z', 'K')

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_GEOMETRY

END MODULE CHECK_GEOMETRY_MODULE
