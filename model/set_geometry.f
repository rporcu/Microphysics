!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: SET_GEOMETRY                                           !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Calculate X, X_E,  oX, oX_E                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_GEOMETRY

! Global Variables:
!---------------------------------------------------------------------//
! Domain decomposition and dimensions
      use geometry, only: DX, oDX
      use geometry, only: DY, oDZ
      use geometry, only: DZ, oDY
! Domain indices.
      use geometry, only: DO_I, IMIN1, IMAX1, IMAX2, IMAX3, IMIN3
      use geometry, only: DO_J, JMIN1, JMAX1, JMAX2, JMAX3, JMIN3
      use geometry, only: DO_K, KMIN1, KMAX1, KMAX2, KMAX3, KMIN3
! Cyclic domain flags.
      use geometry, only: CYCLIC
      use geometry, only: CYCLIC_X, CYCLIC_X_PD, CYCLIC_X_MF
      use geometry, only: CYCLIC_Y, CYCLIC_Y_PD, CYCLIC_Y_MF
      use geometry, only: CYCLIC_Z, CYCLIC_Z_PD, CYCLIC_Z_MF
! MPI-Domain decompoint and rank flags.
      use compar, only: NODESI, NODESJ, NODESK
! Flag for specificed constant mass flux.
      use bc, only: Flux_g

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ZERO, HALF, ONE, UNDEFINED

! Module procedures
!---------------------------------------------------------------------//
      use error_manager
      use toleranc

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("SET_GEOMETRY")

! Allocate geometry arrays.
      CALL ALLOCATE_ARRAYS_GEOMETRY

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
