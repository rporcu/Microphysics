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

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

! Modules
!---------------------------------------------------------------------//
      Use param, only: dim_m, dimension_c
!---------------------------------------------------------------------//

! Coefficients for calibrating Syamlal-O'Brien drag correlation with
! Umf data
      real(c_real) :: drag_c1, drag_d1

! Gravitational acceleration
      real(c_real) :: GRAVITY(3)

! Universal gas constant; (Pa.m3/kmol.K)
      real(c_real), parameter :: GAS_CONST = 8314.56D0

! Pi, the ubiquitous irrational number
      real(c_real), PARAMETER :: Pi = 4.D0*ATAN(1.D0)

! Maximum pressure correction allowed in one iteration
      real(c_real) :: MAX_DELP

! User defined constants
      real(c_real) :: C (DIMENSION_C)

! Names of user defined constants (for output file only)
      CHARACTER(LEN=20) :: C_NAME (DIMENSION_C)


! Solids phase constants
!-----------------------------------------------------------------------
! Number of solids phases
      integer :: mmax = 0
! Particle diameters
      real(c_real) :: d_p0(dim_m)
! Constant solids phase densities.
      real(c_real) :: ro_s0(dim_m)

      END MODULE constant
