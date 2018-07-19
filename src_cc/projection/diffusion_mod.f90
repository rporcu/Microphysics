! 
!              
!  This module contains the subroutines to compute the three components
!  of the diffusion term div(tau) where
!
!      tau = mu ( grad(u) + grad(u)^T )
!  
!  Author: Michele Rosso
! 
!  Date: October 12, 2017
!
! 
module diffusion_mod
   
   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one
   use bc,                only: minf_, nsw_, fsw_, psw_
   
   implicit none
   private

   public compute_divtau
   
contains
   !
   ! Computes  d(txx)/dx + d(txy)/dy + d(txz)/dz 
   !
   !  txx = 2 * mu * du/dx 
   !  txy =  mu * ( du/dy + dv/dx ) 
   !  txz =  mu * ( du/dz + dw/dx )
   ! 
   subroutine compute_divtau ( lo, hi, divtau, dlo, dhi, & 
                               vel_in, vlo, vhi, &
                               mu, lambda, rop, slo, shi, &
                               domlo, domhi, &
                               bc_ilo_type, bc_ihi_type, &
                               bc_jlo_type, bc_jhi_type, &
                               bc_klo_type, bc_khi_type, dx, ng ) bind(C)


      ! Loops bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Number of ghost cells 
      integer(c_int),  intent(in   ) :: ng

      ! Array bounds
      integer(c_int),  intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: slo(3), shi(3)
      integer(c_int),  intent(in   ) :: dlo(3), dhi(3)
      integer(c_int),  intent(in   ) :: domlo(3), domhi(3)

      ! Grid
      real(ar),        intent(in   ) :: dx(3)

      ! Arrays
      real(ar),        intent(in   ) ::                        &
           & vel_in(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3), &
           &    rop(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           &     mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           &  lambda(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      
      real(ar),        intent(inout) ::                        &
           & divtau(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),3)

      ! BC types
      integer(c_int), intent(in   ) ::  &
           & bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           & bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

      ! Temporary array just to handle bc's
      real(ar) &
           & vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)

      integer(c_int)                 :: i, j, k, n
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: txx, tyy, tzz
      
      real(ar)                       :: mu_e, mu_w
      real(ar)                       :: mu_n, mu_s
      real(ar)                       :: mu_t, mu_b
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      do k = lo(3)-1, hi(3)+1
         do j = lo(2)-1, hi(2)+1
            do i = lo(1)-1, hi(1)+1
                vel(i,j,k,:) = vel_in(i,j,k,:)
            end do
         end do
      end do

      ! 
      ! Put values into ghost cells so we can use the easy derivatives
      ! 
      if ( lo(1) == domlo(1) ) then
         i = lo(1)
         do n = 1, 3
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)

               if ( ( bc_ilo_type(j,k,1) == MINF_ ) .or. &
                    ( bc_ilo_type(j,k,1) == NSW_ )  .or. &
                    ( bc_ilo_type(j,k,1) == FSW_ )  .or. &
                    ( bc_ilo_type(j,k,1) == PSW_ )  ) then

                  vel(lo(1)-1,j,k,n) = 2.d0*vel_in(lo(1)-1,j,k,n) - vel_in(lo(1),j,k,n)

               end if
            end do
         end do
         end do
      end if

      if ( hi(1) == domhi(1) ) then

         i = hi(1) 

         do n = 1, 3
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)


               if ( ( bc_ihi_type(j,k,1) == MINF_ ) .or. &
                    ( bc_ihi_type(j,k,1) == NSW_ )  .or. &
                    ( bc_ihi_type(j,k,1) == FSW_ )  .or. &
                    ( bc_ihi_type(j,k,1) == PSW_ )  ) then

                  vel(hi(1)+1,j,k,n) = 2.d0*vel_in(hi(1)+1,j,k,n) - vel_in(hi(1),j,k,n)

               end if
            end do
         end do
         end do
      end if

      if ( lo(2) == domlo(2) ) then
         j = lo(2)
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)

               if ( ( bc_jlo_type(i,k,1) == MINF_ ) .or. &
                    ( bc_jlo_type(i,k,1) == NSW_ )  .or. &
                    ( bc_jlo_type(i,k,1) == FSW_ )  .or. &
                    ( bc_jlo_type(i,k,1) == PSW_ )  ) then

                  vel(i,lo(2)-1,k,:) = 2.d0*vel_in(i,lo(2)-1,k,:) - vel_in(i,lo(2),k,:)

               end if
            end do
         end do
      end if

      if ( hi(2) == domhi(2) ) then

         j = hi(2) 

         do n = 1, 3
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)


               if ( ( bc_jhi_type(i,k,1) == MINF_ ) .or. &
                    ( bc_jhi_type(i,k,1) == NSW_ )  .or. &
                    ( bc_jhi_type(i,k,1) == FSW_ )  .or. &
                    ( bc_jhi_type(i,k,1) == PSW_ )  ) then

                  vel(i,hi(2)+1,k,n) = 2.d0*vel_in(i,hi(2)+1,k,n) - vel_in(i,hi(2),k,n)

               end if
            end do
         end do
         end do
      end if

      if ( lo(3) == domlo(3) ) then

         k = lo(3)

         do n = 1, 3
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               if ( ( bc_klo_type(i,j,1) == MINF_ ) .or. &
                    ( bc_klo_type(i,j,1) == NSW_ )  .or. &
                    ( bc_klo_type(i,j,1) == FSW_ )  .or. &
                    ( bc_klo_type(i,j,1) == PSW_ )  ) then

                  vel(i,j,lo(3)-1,n) = 2.d0*vel_in(i,j,lo(3)-1,n) - vel_in(i,j,lo(3),n)

               end if
            end do
         end do
         end do
      end if

      if ( hi(3) == domhi(3) ) then

         i = hi(1) 

         do n = 1, 3
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               if ( ( bc_khi_type(i,j,1) == MINF_ ) .or. &
                    ( bc_khi_type(i,j,1) == NSW_ )  .or. &
                    ( bc_khi_type(i,j,1) == FSW_ )  .or. &
                    ( bc_khi_type(i,j,1) == PSW_ )  ) then

                  vel(i,j,hi(3)+1,n) = 2.d0*vel_in(i,j,hi(3)+1,n) - vel_in(i,j,hi(3),n)

               end if
            end do
         end do
         end do
      end if

      do n = 1, 3
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               mu_w = half * (mu(i,j,k) + mu(i-1,j,k))
               mu_e = half * (mu(i,j,k) + mu(i+1,j,k))
               mu_s = half * (mu(i,j,k) + mu(i,j-1,k))
               mu_n = half * (mu(i,j,k) + mu(i,j+1,k))
               mu_b = half * (mu(i,j,k) + mu(i,j,k-1))
               mu_t = half * (mu(i,j,k) + mu(i,j,k+1))

               ! txx
               txx = ( mu_e * ( vel(i+1,j,k,n) - vel(i  ,j,k,n) ) &
                      -mu_w * ( vel(i  ,j,k,n) - vel(i-1,j,k,n) ) ) * idx * idx
               tyy = ( mu_n * ( vel(i,j+1,k,n) - vel(i,j  ,k,n) ) &
                      -mu_s * ( vel(i,j  ,k,n) - vel(i,j-1,k,n) ) ) * idy * idy
               tzz = ( mu_t * ( vel(i,j,k+1,n) - vel(i,j,k  ,n) ) &
                      -mu_b * ( vel(i,j,k  ,n) - vel(i,j,k-1,n) ) ) * idz * idz

               ! Assemble divtau
               divtau(i,j,k,n) = (txx + tyy + tzz) / rop(i,j,k)

            end do
         end do
      end do
      end do

   end subroutine compute_divtau
   
end module diffusion_mod
