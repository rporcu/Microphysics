module calc_mu_g_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

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
   subroutine calc_mu_g(slo,shi,lambda_g,mu_g)

      use param, only: is_undefined

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)

! Dummy arguments .....................................................//
      real(c_real), intent(  out) :: lambda_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Local variables .....................................................//
      integer :: i,j,k
      real(c_real), parameter :: f2o3 = 2.d0/3.d0

!-----------------------------------------------------------------------!

      do k = slo(3),shi(3)
         do j = slo(2),shi(2)
            do i = slo(1),shi(1)

                  mu_g(i,j,k) = 1.7d-5 * (293.15d0/273.0d0)**1.5d0 *&
                     (383.d0/(293.15d0+110.d0))
                  lambda_g(i,j,k) = -f2o3*mu_g(i,j,k)

            enddo
         enddo
      enddo

   end subroutine calc_mu_g
end module calc_mu_g_module
