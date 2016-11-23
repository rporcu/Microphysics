!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_GEOMETRY1                                          C
!  Author: M. Syamlal                                 Date: 1-MAY-96   C
!                                                                      C
!  Purpose: Calculate cell volumes and face areas                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_GEOMETRY1

! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE run
      USE geometry
      USE compar
      USE functions

      IMPLICIT NONE


! Local Variables
!-----------------------------------------------
      INTEGER :: I, J, K

      k = kstart3
      j = jstart3
      i = istart3

      VOL = DX*DY*DZ
      AYZ = DY*DZ
      AXY = DX*DY
      AXZ = DX*DZ

      RETURN

      END SUBROUTINE SET_GEOMETRY1
