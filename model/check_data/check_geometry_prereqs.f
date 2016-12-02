!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_GEOMETRY_PREREQS                                  !
!  Purpose: Check the distributed parallel namelist variables.         !
!                                                                      !
!  Author: P. Nicoletti                               Date: 14-DEC-99  !
!  Reviewer: J.Musser                                 Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_GEOMETRY_PREREQS



! Global Variables:
!---------------------------------------------------------------------//
! Domain partitions in various directions.
      use geometry, only: IMAX
      use geometry, only: JMAX
      use geometry, only: KMAX, DZ, ZLENGTH

! Runtime flag specifying 2D simulations

      use geometry, only: COORDINATES

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ONE, ZERO, UNDEFINED_I, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival

      implicit none

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_GEOMETRY_PREREQS")

! Verify that the domain decomposition was specified.
      IF(IMAX == UNDEFINED_I .OR. JMAX == UNDEFINED_I .OR.             &
         KMAX == UNDEFINED_I ) THEN
         WRITE(ERR_MSG,1000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1000 FORMAT('Error 1000: IMAX or JMAX or KMAX not specified in ',     &
          'mfix.dat')


      SELECT CASE(trim(COORDINATES))

      CASE ('CARTESIAN')

      CASE DEFAULT
         WRITE(ERR_MSG, 1103)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1103 FORMAT('Error 1103: Unknown COORDINATES specified. Please ',     &
         'correct the ',/'mfix.dat file.')

      END SELECT

      CALL FINL_ERR_MSG

      RETURN



      END SUBROUTINE CHECK_GEOMETRY_PREREQS
