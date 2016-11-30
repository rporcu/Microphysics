!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: fld_constants                                               C
!  Purpose: Common block containing field variable constants           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      module fld_const

         ! Specified constant gas density
         double precision :: ro_g0

         ! Specified constant gas viscosity
         double precision :: mu_g0

         ! Average molecular weight of gas
         double precision :: mw_avg

      end module fld_const
