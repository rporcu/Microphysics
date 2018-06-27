!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_epg                                                 C
!  Purpose: This subroutine hard-wires the volume fraction in cases    C
!  where we want to use this for debugging.                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

   subroutine set_epg ( ep, slo, shi, domlo, domhi ) bind(C,name="set_epg")

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int
      use param         , only: half, one

      integer(c_int), intent(in) :: slo(3), shi(3), domlo(3), domhi(3)
      real(rt), intent(inout)    :: ep(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer :: ih, i

      ih = ( domhi(1) - domlo(1) ) / 2

      do i = slo(1), shi(1)
         if( i<ih ) then
            ep(i,:,:) = half
         else
            ep(i,:,:) = one
         endif
      end do

   end subroutine set_epg
