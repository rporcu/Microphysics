!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: SET_GEOMETRY                                           !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Calculate X, X_E,  oX, oX_E                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_GEOMETRY(flag)

! Global Variables:
!---------------------------------------------------------------------//
! Domain decomposition and dimensions
      use geometry, only: DX, oDX
      use geometry, only: DY, oDZ
      use geometry, only: DZ, oDY
! Cyclic domain flags.
      use geometry, only: CYCLIC
      use geometry, only: CYCLIC_X, CYCLIC_X_PD, CYCLIC_X_MF
      use geometry, only: CYCLIC_Y, CYCLIC_Y_PD, CYCLIC_Y_MF
      use geometry, only: CYCLIC_Z, CYCLIC_Z_PD, CYCLIC_Z_MF
! Flag for specificed constant mass flux.
      use bc, only: Flux_g

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: HALF, ONE, UNDEFINED

! Module procedures
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar
      use allocate_mod, only: allocate_arrays_geometry

      IMPLICIT NONE

      INTEGER, DIMENSION(:,:,:,:), allocatable :: FLAG

! Local Variables:
!---------------------------------------------------------------------//
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("SET_GEOMETRY")

! Allocate geometry arrays.
      CALL ALLOCATE_ARRAYS_GEOMETRY(flag)

!  Determine the cyclic direction with a specified mass flux
      CYCLIC_X_MF = (FLUX_G /= UNDEFINED .AND. CYCLIC_X_PD)
      CYCLIC_Y_MF = (FLUX_G /= UNDEFINED .AND. CYCLIC_Y_PD)
      CYCLIC_Z_MF = (FLUX_G /= UNDEFINED .AND. CYCLIC_Z_PD)

! Force the cyclic flag if cyclic with pressure drop.
      IF (CYCLIC_X_PD) CYCLIC_X = .TRUE.
      IF (CYCLIC_Y_PD) CYCLIC_Y = .TRUE.
      IF (CYCLIC_Z_PD) CYCLIC_Z = .TRUE.
      CYCLIC = CYCLIC_X .OR. CYCLIC_Y .OR. CYCLIC_Z

      ODX = ONE/DX
      ODY = ONE/DY
      ODZ = ONE/DZ

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_GEOMETRY
