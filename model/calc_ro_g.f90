module calc_ro_g_module
contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  subroutine: set_ro_g                                                !
!                                                                      !
!  Purpose: Initialize the gas density.                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_ro_g(ro_g,rop_g,p_g,ep_g,flag)

      use compar   , only: istart3, jstart3, kstart3, iend3, jend3, kend3
      use eos      , only: eosg
      use fld_const, only: mw_avg, ro_g0

      IMPLICIT NONE

      double precision, intent(inout) ::  ro_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) ::   p_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) ::  ep_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      integer, intent(in   ) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer :: i,j,k
!-----------------------------------------------


      do k = kstart3, kend3
         do j = jstart3, jend3
            do i = istart3, iend3
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
