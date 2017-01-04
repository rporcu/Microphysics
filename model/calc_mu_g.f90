module calc_mu_g_module
contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_mu_g                                               !
!                                                                      !
!  Purpose: Calculate the molecular viscosity in (Pa.s) using the      !
!  Sutherland formulation with Sutherland constant (C) given by        !
!  Vogel's equation C=1.47*Tb.                                         !
!                                                                      !
! For air, Tb=74.82 so that mu = 1.71*10-4 poise at T = 273K           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_mu_g(lambda_g,mu_g,flag)

      use param1, only: zero, is_undefined
      use compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

      implicit none

! Dummy arguments .....................................................//
      double precision, intent(  out) :: lambda_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(  out) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      integer         , intent(in   ) :: flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

! Local variables .....................................................//
      integer :: i,j,k
      double precision, parameter :: f2o3 = 2.d0/3.d0

!-----------------------------------------------------------------------!

      do k = kstart3, kend3
         do j = jstart3, jend3
            do i = istart3, iend3

               if(flag(i,j,k,1) == 1) then
                  mu_g(i,j,k) = 1.7d-5 * (293.15d0/273.0d0)**1.5d0 *&
                     (383.d0/(293.15d0+110.d0))
                  lambda_g(i,j,k) = -f2o3*mu_g(i,j,k)
               endif
            enddo
         enddo
      enddo

      return
   end subroutine calc_mu_g
end module calc_mu_g_module
