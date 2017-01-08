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
      USE geometry, only: IMAX
      USE geometry, only: IMIN1, IMIN2
      USE geometry, only: IMAX1, IMAX2

! Y-Axix partition specifications.
      USE geometry, only: JMAX
      USE geometry, only: JMIN1, JMIN2
      USE geometry, only: JMAX1, JMAX2

! Z-Axix partition specifications.
      USE geometry, only: KMAX
      USE geometry, only: KMIN1, KMIN2
      USE geometry, only: KMAX1, KMAX2

      IMPLICIT NONE

! Set the domain specific values.
      IMIN2 = 1
      IMIN1 = 2
      IMAX1 = IMAX + 1
      IMAX2 = IMAX + 2

! Set the domain specific values.
      JMIN2 = 1
      JMIN1 = 2
      JMAX1 = JMAX + 1
      JMAX2 = JMAX + 2

! Set the domain specific values.
      KMIN2 = 1
      KMIN1 = 2
      KMAX1 = KMAX + 1
      KMAX2 = KMAX + 2

      RETURN
      END SUBROUTINE SET_MAX2
END MODULE SET_MAX2_MODULE
