      MODULE param1

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

! Maximum number of cell classes
      INTEGER, PARAMETER :: MAX_CLASS = 1000000
! Maximum number of corner cells
      INTEGER, PARAMETER :: MAX_NCORN = 4000

! Dimension for the upper triangle of an MxM matrix
      INTEGER :: DIMENSION_LM
! Total number of species
      INTEGER :: DIMENSION_N_all

! Parameters for testing if user input was specifed.
      real(c_real), PARAMETER :: UNDEFINED = 9.87654321D31
      INTEGER, PARAMETER :: UNDEFINED_I = 987654321
      CHARACTER, PARAMETER :: UNDEFINED_C = ' '

! Cutoffs for large and small numbers
      real(c_real), PARAMETER :: LARGE_NUMBER = 1.0D32
      real(c_real), PARAMETER :: SMALL_NUMBER = 1.0D-15

! ZERO, HALF, ONE
      real(c_real), PARAMETER :: ZERO = 0.0d0
      real(c_real), PARAMETER :: HALF = 0.5d0
      real(c_real), PARAMETER :: ONE  = 1.0d0

      interface is_defined
         module procedure is_defined_db
         module procedure is_defined_i
      end interface is_defined

      interface is_undefined
         module procedure is_undefined_db
         module procedure is_undefined_i
      end interface is_undefined

   CONTAINS

      PURE LOGICAL FUNCTION IS_DEFINED_DB(x)
         real(c_real), INTENT(IN) :: x
         IS_DEFINED_DB = (x /= UNDEFINED)
      END FUNCTION IS_DEFINED_DB

      PURE LOGICAL FUNCTION IS_DEFINED_I(x)
         INTEGER, INTENT(IN) :: x
         IS_DEFINED_I = (x /= UNDEFINED_I)
      END FUNCTION IS_DEFINED_I

      PURE LOGICAL FUNCTION IS_UNDEFINED_DB(x)
         real(c_real), INTENT(IN) :: x
         IS_UNDEFINED_DB = (x == UNDEFINED)
      END FUNCTION IS_UNDEFINED_DB

      PURE LOGICAL FUNCTION IS_UNDEFINED_I(x)
         INTEGER, INTENT(IN) :: x
         IS_UNDEFINED_I = (x == UNDEFINED_I)
      END FUNCTION IS_UNDEFINED_I

      END MODULE param1
