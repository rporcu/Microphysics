module calc_ro_g_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  subroutine: set_ro_g                                                !
!                                                                      !
!  Purpose: Initialize the gas density.                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_ro_g(slo, shi, lo, hi, ro_g, rop_g, p_g, ep_g)

      use eos      , only: eosg
      use fld_const, only: ro_g0, mw_avg
      use param,     only: is_undefined

      implicit none

! Dummy arguments ....................................................//
      integer(c_int), intent(in   ) :: slo(3), shi(3), lo(3), hi(3)

      real(c_real), intent(in   ) :: &
          p_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
         ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) ::  &
          ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
         rop_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Local variables .....................................................//
      integer :: i,j,k

      if ( is_undefined( ro_g0 ) ) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ro_g(i,j,k)  = eosg(mw_avg,p_g(i,j,k),295.15d0)
                  rop_g(i,j,k) = ep_g(i,j,k) * ro_g(i,j,k)
               enddo
            enddo
         enddo
         
      else

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ro_g(i,j,k)  = ro_g0
                  rop_g(i,j,k) = ep_g(i,j,k) * ro_g0
               enddo
            enddo
         enddo
         
      endif

   end subroutine calc_ro_g

end module calc_ro_g_module
