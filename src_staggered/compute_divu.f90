   !
   ! Compute the divergence of the velocity field
   ! It is equivalent to cal_trd_g but the arguments are arranged
   ! in a way to make it easier to pass in amrex stuff 
   ! 
   subroutine compute_divu ( lo, hi, divu, slo, shi, u_g, ulo, uhi, v_g, vlo, vhi, &
        & w_g, wlo, whi, dx )  bind(C, name="compute_divu")

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

      real(rt),   intent(in   ) :: dx(3)
      real(rt),   intent(  out) :: &
           divu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      
      real(rt), intent(in   ) :: &
           u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      ! Local variables
      !-----------------------------------------------
      integer      :: i, j, k
      real(rt) :: odx, ody, odz

      odx = 1.d0 / dx(1)
      ody = 1.d0 / dx(2)
      odz = 1.d0 / dx(3)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               divu(i,j,k) = &
                    (u_g(i+1,j,k)-u_g(i,j,k))*odx + &
                    (v_g(i,j+1,k)-v_g(i,j,k))*ody  + &
                    (w_g(i,j,k+1)-w_g(i,j,k))*odz
            end do
         end do
      end do

   end subroutine compute_divu
