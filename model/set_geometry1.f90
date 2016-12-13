!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_GEOMETRY1                                          C
!  Author: M. Syamlal                                 Date: 1-MAY-96   C
!                                                                      C
!  Purpose: Calculate cell volumes and face areas                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_GEOMETRY1

      USE compar, only: istart3, jstart3, kstart3
      USE geometry, only: dx, dy, dz, ayz, axy, axz, vol

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
