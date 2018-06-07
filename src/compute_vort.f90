module compute_vort_module

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   implicit none

contains

   !
   ! Compute the vorticity 
   ! 
   subroutine compute_vort ( lo, hi, vort, slo, shi, u_g, ulo, uhi, v_g, vlo, vhi, &
        & w_g, wlo, whi, dx )  bind(C, name="compute_vort")

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

      real(rt),   intent(in   ) :: dx(3)
      real(rt),   intent(  out) :: &
           vort(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      
      real(rt), intent(in   ) :: &
           u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      ! Local variables
      !-----------------------------------------------
      integer      :: i, j, k
      real(rt) :: odx, ody, odz
      real(rt) :: uy,uz,vx,vz,wx,wy

      odx = 1.d0 / dx(1)
      ody = 1.d0 / dx(2)
      odz = 1.d0 / dx(3)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               uy = 0.25d0 * ody * ( u_g(i+1,j+1,k) + u_g(i,j+1,k) &
                                    -u_g(i+1,j-1,k) - u_g(i,j-1,k))
               uz = 0.25d0 * odz * ( u_g(i+1,j,k+1) + u_g(i,j,k+1) &
                                    -u_g(i+1,j,k-1) - u_g(i,j,k-1))
               vx = 0.25d0 * odx * ( v_g(i+1,j+1,k) + v_g(i+1,j,k) &
                                    -v_g(i-1,j+1,k) - v_g(i-1,j,k))
               vz = 0.25d0 * odz * ( v_g(i,j+1,k+1) + v_g(i,j,k+1) &
                                    -v_g(i,j+1,k-1) - v_g(i,j,k-1))
               wx = 0.25d0 * odx * ( w_g(i+1,j,k+1) + w_g(i+1,j,k) &
                                    -w_g(i-1,j,k+1) - w_g(i-1,j,k))
               wy = 0.25d0 * ody * ( w_g(i,j+1,k+1) + w_g(i,j+1,k) &
                                    -w_g(i,j-1,k+1) - w_g(i,j-1,k))
               vort(i,j,k) = &
                   sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)

            end do
         end do
      end do

   end subroutine compute_vort

end module compute_vort_module
