module calc_grad_des_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   contains

      subroutine calc_grad_des(slo, shi, lo, hi, PHI, DEL_PHI, dx, dy, dz, domlo, domhi)

! Modules
!-----------------------------------------------

      use functions, only: avg
      use param, only: zero

      IMPLICIT NONE

      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) ::  lo(3), hi(3)
      integer, intent(in   ) :: domlo(3), domhi(3)

      real(c_real), intent(in   ) :: PHI&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: DEL_PHI&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)
      real(c_real), intent(in   ) :: dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! general i, j, k indices
      integer :: I, J, K
      real(c_real) :: odx, ody, odz
!......................................................................!

      odx = 1.d0 / dx
      ody = 1.d0 / dy
      odz = 1.d0 / dz

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               del_phi(i,j,k,:) = zero

               if((i > domlo(1)).and.(i < domhi(1))) then
                  del_phi(i,j,k,1) = odx*(avg(phi(i,j,k),phi(i+1,j,k)) -    &
                     avg(phi(i-1,j,k),phi(i,j,k)))
               elseif(i == domlo(1)) then
                  del_phi(i,j,k,1) = 2.0d0*odx *                            &
                     (avg(phi(i,j,k),phi(i+1,j,k)) -  phi(i,j,k))
               elseif(i == domhi(1)) then
                  del_phi(i,j,k,1) = 2.0d0*odx *                            &
                     (phi(i,j,k) - avg(phi(i-1,j,k), phi(i,j,k)))
               else
                  del_phi(i,j,k,1) = zero
               endif


               if((j > domlo(2)) .and. (j < domhi(2))) then
                  del_phi(i,j,k,2) = ody*(avg(phi(i,j,k),phi(i,j+1,k)) -    &
                     avg(phi(i,j-1,k),phi(i,j,k)))
               elseif(j == domlo(2)) then
                  del_phi(i,j,k,2) = 2.0d0*ody *                            &
                     (avg(phi(i,j,k),phi(i,j+1,k)) - phi(i,j,k))
               elseif(j == domhi(2)) then
                  del_phi(i,j,k,2) = 2.0d0*ody *                            &
                     (phi(i,j,k)- avg(phi(i,j-1,k),phi(i,j,k)))
               else
                  del_phi(i,j,k,2) = zero
               endif


               if((k > domlo(3)) .and. (k < domhi(3))) then
                  del_phi(i,j,k,3) = odz*(avg(phi(i,j,k),phi(i,j,k+1)) -    &
                     avg(phi(i,j,k-1),phi(i,j,k)))
               elseif(k == domlo(3)) then
                  del_phi(i,j,k,3) = 2.0d0*odz *                            &
                     (avg(phi(i,j,k),phi(i,j,k+1)) - phi(i,j,k))
               elseif(k == domhi(3)) then
                  del_phi(i,j,k,3) = 2.0d0*odz *                            &
                     (phi(i,j,k) - avg(phi(i,j,k-1),phi(i,j,k)))
               else
                  del_phi(i,j,k,3) = zero
               endif
            enddo
         enddo
      enddo

      end subroutine calc_grad_des
   end module calc_grad_des_module
