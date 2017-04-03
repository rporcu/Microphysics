! This file is included by functions_mod.f
! These function are are defined in this file so it can be included by
! individual subprograms in order to inline the (short) functions
! defined in this file.

! All functions that may benefits from inlining should be defined in this file.

! For code that is not performance-critical, do:

! SUBROUTINE subname
! use functions
! ...
! END SUBROUTINE subname

! For code that is performance-critical, do:

! SUBROUTINE subname
! ...
! CONTAINS
! INCLUDE 'functions.inc'
! END SUBROUTINE subname

!---------------------------------------------------------------------//
! Misc
!---------------------------------------------------------------------//

! Arithmetic average
      DOUBLE PRECISION FUNCTION AVG(XXXm, XXXp)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) ::XXXp, XXXm
      AVG = 0.5d0 *(XXXm + XXXp)
      END FUNCTION AVG

! Harmonic average
      DOUBLE PRECISION FUNCTION AVG_H(XXXm, XXXp)
      use param1, only: SMALL_NUMBER
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      AVG_H = XXXm * XXXp / MAX(SMALL_NUMBER, 0.5d0*(XXXm + XXXp) )
      END FUNCTION AVG_H
