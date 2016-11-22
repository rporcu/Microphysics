! This file is included by fun_avg_mod.f

! These function are are defined in this file so it can be included by
! individual subprograms in order to inline the (short) functions
! defined in this file.

! All functions that may benefits from inlining should be defined in this file.

! For code that is not performance-critical, do:

! SUBROUTINE subname
! USE fun_avg
! ...
! END SUBROUTINE subname

! For code that is performance-critical, do:

! SUBROUTINE subname
! ...
! CONTAINS
! INCLUDE 'fun_avg.inc'
! END SUBROUTINE subname

! Arithmetic average
      DOUBLE PRECISION FUNCTION AVG(XXXm, XXXp)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) ::XXXp, XXXm
      AVG = 0.5d0 *(XXXm + XXXp)
      END FUNCTION AVG

! Harmonic average
      DOUBLE PRECISION FUNCTION AVG_H(XXXm, XXXp)
      USE param1, only: SMALL_NUMBER
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      AVG_H = XXXm * XXXp / MAX(SMALL_NUMBER, 0.5d0*(XXXm + XXXp) )
      END FUNCTION AVG_H
