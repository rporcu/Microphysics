!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: constant                                                    C
!  Author: M. Syamlal                                 Date: 5-FEB-92   C
!                                                                      C
!  Purpose: Common block containing physical constants and constants   C
!           used in the numerical technique                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE constant


! Modules
!---------------------------------------------------------------------//
      Use param, only: dim_m, dimension_c
!---------------------------------------------------------------------//

! Coefficients for calibrating Syamlal-O'Brien drag correlation with
! Umf data
      DOUBLE PRECISION :: drag_c1, drag_d1

! Gravitational acceleration
      DOUBLE PRECISION :: GRAVITY(3)

! Universal gas constant; (Pa.m3/kmol.K)
      DOUBLE PRECISION, parameter :: GAS_CONST = 8314.56D0

! Pi, the ubiquitous irrational number
      DOUBLE PRECISION, PARAMETER :: Pi = 4.D0*ATAN(1.D0)

! Maximum pressure correction allowed in one iteration
      DOUBLE PRECISION :: MAX_DELP

! User defined constants
      DOUBLE PRECISION :: C (DIMENSION_C)

! Names of user defined constants (for output file only)
      CHARACTER(LEN=20) :: C_NAME (DIMENSION_C)


! Solids phase constants
!-----------------------------------------------------------------------
! Number of solids phases
      integer :: mmax = 0
! Particle diameters
      double precision :: d_p0(dim_m)
! Constant solids phase densities.
      double precision :: ro_s0(dim_m)

      END MODULE constant
