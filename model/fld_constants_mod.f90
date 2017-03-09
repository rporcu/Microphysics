!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: fld_constants                                               C
!  Purpose: Common block containing field variable constants           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      module fld_const

         use amrex_fort_module, only : c_real => amrex_real
         use iso_c_binding , only: c_int

         ! Specified constant gas density
         real(c_real) :: ro_g0

         ! Specified constant gas viscosity
         real(c_real) :: mu_g0

         ! Average molecular weight of gas
         real(c_real) :: mw_avg

      end module fld_const
