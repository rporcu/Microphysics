      MODULE param1

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

! Maximum number of corner cells
      integer, parameter :: MAX_NCORN = 4000

! Total number of species
      integer :: DIMENSION_N_all

! Parameters for testing if user input was specifed.
      real(c_real), parameter :: UNDEFINED = 9.87654321D31
      integer, parameter :: UNDEFINED_I = 987654321
      CHARACTER, parameter :: UNDEFINED_C = ' '

! Cutoffs for large and small numbers
      real(c_real), parameter :: LARGE_NUMBER = 1.0D32
      real(c_real), parameter :: SMALL_NUMBER = 1.0D-15

! ZERO, HALF, ONE
      real(c_real), parameter :: ZERO = 0.0d0
      real(c_real), parameter :: HALF = 0.5d0
      real(c_real), parameter :: ONE  = 1.0d0

      interface is_defined
         module procedure is_defined_db
         module procedure is_defined_i
      end interface is_defined

      interface is_undefined
         module procedure is_undefined_db
         module procedure is_undefined_i
      end interface is_undefined

   contains

      pure logical function is_defined_db(x)
         real(c_real), intent(in) :: x
         is_defined_db = .not.equal(x, undefined)
      end function is_defined_db

      pure logical function is_defined_i(x)
         integer, intent(in) :: x
         is_defined_i = (x /= undefined_i)
      end function is_defined_i

      pure logical function is_undefined_db(x)
         real(c_real), intent(in) :: x
         is_undefined_db = equal(x, undefined)
      end function is_undefined_db

      pure logical function is_undefined_i(x)
         integer, intent(in) :: x
         is_undefined_i = (x == undefined_i)
      end function is_undefined_i

      pure logical function equal(x, y)
         real(c_real), intent(in) :: x, y
         equal = (abs(x-y) < epsilon(x))
      end function equal

      end MODULE param1
