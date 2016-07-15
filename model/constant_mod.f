!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: constant                                                    C
!  Purpose: Common block containing physical constants and constants   C
!           used in the numerical technique                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 5-FEB-92   C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!    Gera, D., Syamlal, M., and O'Brien, T. J., "Hydrodynamics of      C
!      particle segregation in fluidized beds", Int. J. of Multiphase  C
!      Flow, Vol 30, 2004, pp. 419-428.                                C
!    Johnson, P. C., and Jackson, R., "Frictional-collisional          C
!      constitutive relations for granluar materials, with application C
!      to plane shearing", JFM, Vol. 176, 1987, pp. 67-93.             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE constant


! Modules
!---------------------------------------------------------------------//
      Use param, only: dim_m, dimension_c
!---------------------------------------------------------------------//

! Packed bed (close packed) void fraction
      DOUBLE PRECISION :: EP_star

! Coefficients for calibrating Syamlal-O'Brien drag correlation with
! Umf data
      DOUBLE PRECISION :: drag_c1, drag_d1

! UNIT conversion factor for pressure (Barye to Pa if SI)
      DOUBLE PRECISION :: to_SI

! Gravitational acceleration
      DOUBLE PRECISION :: GRAVITY, GRAVITY_X, GRAVITY_Y, GRAVITY_Z

! Universal gas constant
      DOUBLE PRECISION :: GAS_CONST

! Universal gas constant in cal/mol.K
      DOUBLE PRECISION, PARAMETER :: GAS_CONST_cal = 1.987207D0

! Pi, the ubiquitous irrational number
      DOUBLE PRECISION, PARAMETER :: Pi = 4.D0*ATAN(1.D0)

! Square root of Pi
      DOUBLE PRECISION, PARAMETER :: SQRT_Pi = 2.D0*SQRT(ATAN(1.D0))

! Maximum pressure correction allowed in one iteration
      DOUBLE PRECISION :: MAX_DELP

! User defined constants
      DOUBLE PRECISION :: C (DIMENSION_C)

! Names of user defined constants (for output file only)
      CHARACTER(LEN=20) :: C_NAME (DIMENSION_C)


      END MODULE constant
