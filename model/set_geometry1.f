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

      VOL = DX(I)*DY(J)*(X(I)*DZ(K))
      AYZ = DY(J)*(X_E(I)*DZ(K))
      AXY = DX(I)*DY(J)
      AXZ = DX(I)*(X(I)*DZ(K))

      RETURN

      END SUBROUTINE SET_GEOMETRY1
