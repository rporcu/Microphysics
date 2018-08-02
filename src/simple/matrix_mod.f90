      module matrix

      use iso_c_binding , only: c_int

      ! Definitions for sparse matrix
      integer(c_int), parameter :: e = 1
      integer(c_int), parameter :: w =-1
      integer(c_int), parameter :: n = 2
      integer(c_int), parameter :: s =-2
      integer(c_int), parameter :: t = 3
      integer(c_int), parameter :: b =-3

      end module matrix
