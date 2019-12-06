!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: random_nb_mod.f                                        C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      module random_number_module

      use amrex_fort_module, only : rt => amrex_real, amrex_random
      use iso_c_binding , only: c_int

      implicit none

      contains

      real(rt) function get_random() bind(C)
      implicit none
      get_random = amrex_random()
      end function get_random

      end module random_number_module
