MODULE SET_MAX2_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: SET_MAX2                                               !
!  Purpose: calculate domain bounds.                                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_MAX2

      use geometry, only: imax, jmax, kmax
      use geometry, only: domlo, domhi

      IMPLICIT NONE

      domlo(1) = 0 + 2
      domlo(2) = 0 + 2
      domlo(3) = 0 + 2

      domhi(1) = imax-1 + 2
      domhi(2) = jmax-1 + 2
      domhi(3) = kmax-1 + 2

      RETURN
      END SUBROUTINE SET_MAX2
END MODULE SET_MAX2_MODULE
