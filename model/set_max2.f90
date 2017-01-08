MODULE SET_MAX2_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: SET_MAX2                                               !
!  Purpose: calculate domain bounds.                                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_MAX2

! X-Axix partition specifications.
      USE compar, only: NODESI
      USE geometry, only: IMAX
      USE geometry, only: IMIN1, IMIN2, IMIN3, IMIN4
      USE geometry, only: IMAX1, IMAX2, IMAX3, IMAX4

! Y-Axix partition specifications.
      USE compar, only: NODESJ
      USE geometry, only: JMAX
      USE geometry, only: JMIN1, JMIN2, JMIN3, JMIN4
      USE geometry, only: JMAX1, JMAX2, JMAX3, JMAX4

! Z-Axix partition specifications.
      USE compar, only: NODESK
      USE geometry, only: KMAX
      USE geometry, only: KMIN1, KMIN2, KMIN3, KMIN4
      USE geometry, only: KMAX1, KMAX2, KMAX3, KMAX4

      IMPLICIT NONE


! Initialize I's
      IMIN1=1;  IMIN2=1;  IMIN3=1;  IMIN4=1
      IMAX1=1;  IMAX2=1;  IMAX3=1;  IMAX4=1

! Set the domain specific values.
      IMIN1 = 2
      IMAX1 = IMAX + 1
      IMAX2 = IMAX + 2
      IMIN2 = 1
      IF(NODESI.NE.1) THEN
         IMIN3 = 0
         IMAX3 = IMAX + 3
         IMIN4 = -1
         IMAX4 = IMAX + 4
      ELSE
         IMIN3 = IMIN2
         IMAX3 = IMAX2
         IMIN4 = IMIN3
         IMAX4 = IMAX3
      ENDIF

! Initialize J's
      JMIN1=1;  JMIN2=1;  JMIN3=1;  JMIN4=1
      JMAX1=1;  JMAX2=1;  JMAX3=1;  JMAX4=1

! Set the domain specific values.
      JMIN1 = 2
      JMAX1 = JMAX + 1
      JMAX2 = JMAX + 2
      JMIN2 = 1
      IF(NODESJ.NE.1) THEN
         JMIN3 = 0
         JMAX3 = JMAX + 3
         JMIN4 = -1
         JMAX4 = JMAX + 4
      ELSE
         JMIN3 = JMIN2
         JMAX3 = JMAX2
         JMIN4 = JMIN3
         JMAX4 = JMAX3
      ENDIF


! Initialize J's
      KMIN1=1;  KMIN2=1;  KMIN3=1;
      KMAX1=1;  KMAX2=1;  KMAX3=1;

! Set the domain specific values.
      KMIN1 = 2
      KMAX1 = KMAX + 1
      KMAX2 = KMAX + 2
      KMIN2 = 1
      IF(NODESK.NE.1) THEN
         KMIN3 = 0
         KMAX3 = KMAX + 3
      ELSE
         KMIN3 = KMIN2
         KMAX3 = KMAX2
      ENDIF

      RETURN
      END SUBROUTINE SET_MAX2
END MODULE SET_MAX2_MODULE
