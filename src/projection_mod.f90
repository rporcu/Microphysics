! 
!              
!  This module contains the subroutines to perform some of the steps of the
!  projection method.
!  
!  Author: Michele Rosso
! 
!  Date: October 13, 2017
!
! 
module projection_mod
   
   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one
   
   implicit none
   private

contains
   
   !
   ! Computes  u_i* = u_i + dt * RHS
   ! where RHS is the Right-Hand Side of the momentum equation
   !
   ! Here u_i is the i-th component of velocity.
   ! Previous pressure gradient is not included.
   !
   ! dir = 1,2,3 indicates x, y, z direction respectively.
   !  
   subroutine compute_intermediate_velocity ( lo, hi, u_i_star, ulo, uhi,     &
        & u_i, ugradu_i,  divtau_i, drag_i, f_gds_i, rop, slo, shi, dt, dir ) &
        & bind(C, name="compute_intermediate_velocity")

      use constant, only: gravity
      
      integer(c_int),  intent(in)    :: lo(3),  hi(3)   ! Tile indeces for mf associated to ug
      integer(c_int),  intent(in)    :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: slo(3), shi(3)
      integer(c_int),  intent(in   ) :: dir
      real(ar),        intent(in   ) :: dt
      
      
      real(ar),        intent(inout) ::                           &
           & u_i_star(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           & u_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),      &
           & divtau_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           & ugradu_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), & 
           & drag_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),   &
           & f_gds_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),  &
           & rop(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
           
      integer(c_int)                 :: i, j, k
      real(ar)                       :: rhs, irho, f 
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               irho = half / rop(i,j,k) + half / rop(i-1,j,k) 

               ! Remeber the minus sign between the two parts
               f    = f_gds_i(i,j,k) * u_i(i,j,k) - drag_i(i,j,k)
               
               rhs  = - ugradu_i(i,j,k) + gravity(dir) + &
                    & irho * ( divtau_i(i,j,k) + f )

               u_i_star(i,j,k) = u_i(i,j,k) + dt * rhs      
              
               
            end do
         end do
      end do

   end subroutine compute_intermediate_velocity

      
end module projection_mod
