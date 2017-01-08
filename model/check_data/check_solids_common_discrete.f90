MODULE CHECK_SOLIDS_COMMON_DISCRETE_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_SOLIDS_COMMON_DISCRETE                            !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE


! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Store DES_*_OLD arrays.
      USE discretelement, only: DO_OLD
! Number of solids phases.
      USE constant, only: MMAX
! TFM solids phase diameters and densities. (DEM default)
      USE constant, only: D_p0

! User specified integration method.
      USE discretelement, only: DES_INTG_METHOD
      USE discretelement, only: INTG_ADAMS_BASHFORTH
      USE discretelement, only: INTG_EULER
! Max/Min particle radii
      USE discretelement, only: MAX_RADIUS, MIN_RADIUS

! Subroutine access.
      use constant, only: MMAX

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: undefined

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival, err_msg
      use geometry, only: cyclic_x, cyclic_y, cyclic_z
      use geometry, only: cyclic_x_pd, cyclic_y_pd, cyclic_z_pd

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Loop index.
      INTEGER :: lM  ! Solids phase Index

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE")

      MAX_RADIUS = -UNDEFINED
      MIN_RADIUS =  UNDEFINED

! Determine the maximum particle size in the system (MAX_RADIUS), which
! in turn is used for various tasks
      DO lM=1, MMAX
         MAX_RADIUS = MAX(MAX_RADIUS, 0.5d0*D_P0(lM))
         MIN_RADIUS = MIN(MIN_RADIUS, 0.5d0*D_P0(lM))
      ENDDO

! Check for valid integration method
      SELECT CASE(trim(DES_INTG_METHOD))
      CASE ('EULER')
         INTG_EULER = .TRUE.
         INTG_ADAMS_BASHFORTH = .FALSE.
         !DES_INTG_METHOD_ENUM = 1
      CASE ('ADAMS_BASHFORTH')
         INTG_EULER = .FALSE.
         INTG_ADAMS_BASHFORTH = .TRUE.
         !DES_INTG_METHOD_ENUM = 2
      CASE DEFAULT
         WRITE(ERR_MSG,2020) trim(DES_INTG_METHOD)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2020 FORMAT('Error 2020:Invalid DES_INGT_METHOD: ',A,/'Please ',      &
         'correct the mfix.dat file.')

      END SELECT

      DO_OLD = INTG_ADAMS_BASHFORTH

! Check geometry constrains.
      CALL CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY

      CALL FINL_ERR_MSG


      RETURN

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY                   !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Purpose: Check user input data                                      !
!                                                                      !
!  Comments: Geometry checks were moved here from CHECK_DES_DATA.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE geometry, only: ZLENGTH
! Flag: Use DES E-L model
      USE discretelement, only: DES_CONTINUUM_COUPLED
      USE discretelement, only: MAX_RADIUS

      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival, err_msg

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      real(c_real) :: MIN_DEPTH

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY")


      IF(DES_CONTINUUM_COUPLED)THEN
! Check that the depth of the simulation exceeds the largest particle
! to ensure correct calculation of volume fraction. This is important
! for coupled simulations.
         MIN_DEPTH = 2.0d0*MAX_RADIUS
         IF(ZLENGTH < MIN_DEPTH)THEN
            WRITE(ERR_MSG, 1300)
            CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
         ENDIF
      ENDIF

 1300 FORMAT('Error 1300: The maximum particle diameter exceeds the ', &
         'simulation',/'depth (ZLENGTH). Please correct the mfix.dat ',&
         'file.')

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY

END MODULE CHECK_SOLIDS_COMMON_DISCRETE_MODULE
