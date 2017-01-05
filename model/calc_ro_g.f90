module calc_ro_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  subroutine: set_ro_g                                                !
!                                                                      !
!  Purpose: Initialize the gas density.                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_ro_g(slo, shi, ro_g, rop_g, p_g, ep_g, flag)

      use eos      , only: eosg
      use fld_const, only: mw_avg, ro_g0

      IMPLICIT NONE

! Dummy arguments ....................................................//
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      real(c_real), intent(inout) ::  ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) ::   p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) ::  ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer, intent(in   ) ::  flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

! Local variables .....................................................//
      integer :: i,j,k


      do k = slo(3), shi(3)
         do j = slo(2), shi(2)
            do i = slo(1), shi(1)
               if (flag(i,j,k,1)<100) then
                  ro_g(i,j,k) = eosg(mw_avg,p_g(i,j,k),295.15d0)
                  rop_g(i,j,k) = ep_g(i,j,k)*ro_g(i,j,k)
               endif
            enddo
         enddo
      enddo

      return
   end subroutine calc_ro_g

end module calc_ro_g_module
