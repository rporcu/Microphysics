module physical_prop_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PHYSICAL_PROP                                           !
!                                                                      !
!  Purpose: Calculate the indicated physical properties that vary      !
!           with time if directed to do so by the corresponding flag   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine physical_prop(slo, shi, level, ro_g, p_g, ep_g, rop_g, &
      flag) bind(c, name="physical_prop")

      use calc_ro_g_module, only: calc_ro_g
      use funits, only: unit_log
      use param1, only: is_undefined, zero

      USE compar, only: myPE, PE_IO, numPEs
      USE exit_mod, only: mfix_exit
      use fld_const, only: ro_g0

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: level

      real(c_real), intent(inout) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

!-----------------------------------------------------------------------

! Calculate density only. This is invoked several times within iterate,
! making it the most frequently called.
      if(level == 0) then
         if(is_undefined(ro_g0)) &
            call calc_ro_g(slo, shi, ro_g, rop_g, p_g, ep_g, flag)

! Calculate everything except density. This is called at the start of
! each iteration.
      elseif(level == 1) then

! Calculate everything. This is invoked via calc_coeff_all as part of
! the initialization (before starting the time march) and at the start
! of each step step thereafter.
      elseif(level == 2) then
         if(is_undefined(ro_g0)) &
            call calc_ro_g(slo, shi, ro_g, rop_g, p_g, ep_g, flag)
      endif


      return

   end subroutine physical_prop
end module physical_prop_module
