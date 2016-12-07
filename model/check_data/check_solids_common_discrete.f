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
! Runtime Flag: Generate initial particle configuration.
      USE discretelement, only: GENER_PART_CONFIG
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
! Runtime Flag: Periodic boundaries
      USE discretelement, only: DES_PERIODIC_WALLS
      USE discretelement, only: DES_PERIODIC_WALLS_X
      USE discretelement, only: DES_PERIODIC_WALLS_Y
      USE discretelement, only: DES_PERIODIC_WALLS_Z

! Subroutine access.
      use constant, only: MMAX

      USE run, only: MOMENTUM_X_EQ
      USE run, only: MOMENTUM_Y_EQ
      USE run, only: MOMENTUM_Z_EQ

      use run, only: RUN_TYPE
      use discretelement, only: GENER_PART_CONFIG

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: undefined
      use param, only: dim_m

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


! Turn off the 'continuum' equations for discrete solids if the user
! specified them.  We could make use of these flags.
      MOMENTUM_X_EQ((MMAX+1):DIM_M) = .FALSE.
      MOMENTUM_Y_EQ((MMAX+1):DIM_M) = .FALSE.
      MOMENTUM_Z_EQ((MMAX+1):DIM_M) = .FALSE.

! Derive periodicity from cyclic boundary flags.
      DES_PERIODIC_WALLS_X = CYCLIC_X .OR. CYCLIC_X_PD
      DES_PERIODIC_WALLS_Y = CYCLIC_Y .OR. CYCLIC_Y_PD
      DES_PERIODIC_WALLS_Z = CYCLIC_Z .OR. CYCLIC_Z_PD

      DES_PERIODIC_WALLS = (DES_PERIODIC_WALLS_X .OR.                  &
        DES_PERIODIC_WALLS_Y .OR. DES_PERIODIC_WALLS_Z)


! Overwrite for restart cases.
      IF(TRIM(RUN_TYPE) .NE. 'NEW') GENER_PART_CONFIG = .FALSE.

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

! Check interpolation input.
      CALL CHECK_SOLIDS_COMMON_DISCRETE_INTERP

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
      DOUBLE PRECISION :: MIN_DEPTH

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

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_DISCRETE_INTERP                     !
!  Author: J.Musser                                   Date: 25-Nov-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_INTERP

! Runtime Flag: Invoke gas/solids coupled simulation.
      use discretelement, only: DES_CONTINUUM_COUPLED
! User input for DES interpolation scheme.
      use particle_filter, only: DES_INTERP_SCHEME
! Enumerated interpolation scheme for faster access
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE
      use particle_filter, only: DES_INTERP_GARG
      use particle_filter, only: DES_INTERP_DPVM
      use particle_filter, only: DES_INTERP_GAUSS
! User specified filter width
      use particle_filter, only: DES_INTERP_WIDTH
! Flag: Interpolate continuum fields
      use particle_filter, only: DES_INTERP_MEAN_FIELDS
! Flag: Interplate variables for drag calculation.
      use particle_filter, only: DES_INTERP_ON

      use param1, only: UNDEFINED

      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival, err_msg

      IMPLICIT NONE

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_INTERP")

! Set the interpolation ENUM value.
      SELECT CASE(trim(adjustl(DES_INTERP_SCHEME)))
      CASE ('NONE')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_NONE
! Cannot use interpolation when no scheme is selected.
         IF(DES_INTERP_ON)THEN
            WRITE(ERR_MSG,2001) 'DES_INTERP_ON'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(DES_INTERP_MEAN_FIELDS)THEN
            WRITE(ERR_MSG,2001) 'DES_INTERP_MEAN_FIELDS'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF(DES_CONTINUUM_COUPLED) THEN
         ENDIF

      CASE ('GARG_2012')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_GARG

      CASE ('SQUARE_DPVM')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_DPVM

      CASE ('GAUSS_DPVM')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_GAUSS

      CASE DEFAULT
         WRITE(ERR_MSG,2000) trim(adjustl(DES_INTERP_SCHEME))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      END SELECT

 2000 FORMAT('Error 2000: Invalid DES_INTERP_SCHEME: ',A,/'Please ',   &
         'correct the mfix.dat file.')

 2001 FORMAT('Error 2001: No interpolation scheme specified when ',A,/ &
         'is enabled. Please correct the mfix.dat file.')

      SELECT CASE(DES_INTERP_SCHEME_ENUM)

      CASE(DES_INTERP_NONE)

         IF(DES_INTERP_WIDTH /= UNDEFINED) THEN
            WRITE(ERR_MSG,2100) trim(adjustl(DES_INTERP_SCHEME))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2100 FORMAT('Error 2100: The selected interpolation scheme (',A,') ', &
         'does',/'not support an adjustable interpolation width.',/    &
         'Please correct the input file.')


      CASE(DES_INTERP_GARG)
         DES_INTERP_MEAN_FIELDS= .TRUE.

         IF(DES_INTERP_WIDTH /= UNDEFINED) THEN
            WRITE(ERR_MSG,2100) trim(adjustl(DES_INTERP_SCHEME))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF


      END SELECT

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_INTERP
