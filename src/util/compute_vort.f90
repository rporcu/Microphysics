module compute_vort_module

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   implicit none

contains

   !-----------------------------------------------------------------------!

   !
   ! Compute the vorticity
   !
   subroutine compute_vort ( lo, hi, vort, slo, shi, vel, vlo, vhi, dx) &
      bind(C, name="compute_vort")

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: vlo(3),vhi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      real(rt),   intent(in   ) :: dx(3)
      real(rt),   intent(  out) :: &
         vort(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(rt), intent(in   ) :: &
         vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)

      ! Local variables
      !-----------------------------------------------
      integer      :: i, j, k
      real(rt) :: idx, idy, idz
      real(rt) :: uy, uz, vx, vz, wx, wy

      idx = 1.d0 / dx(1)
      idy = 1.d0 / dx(2)
      idz = 1.d0 / dx(3)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               uy = idy * ( vel(i,j+1,k,1) - vel(i,j-1,k,1))
               uz = idz * ( vel(i,j,k+1,1) - vel(i,j,k-1,1))

               vx = idz * ( vel(i+1,j,k,2) - vel(i-1,j,k,2))
               vz = idz * ( vel(i,j,k+1,2) - vel(i,j,k-1,2))

               wx = idz * ( vel(i+1,j,k,3) - vel(i-1,j,k,3))
               wy = idy * ( vel(i,j+1,k,3) - vel(i,j-1,k,3))

               ! The factor half is included here instead of in each of the above
               vort(i,j,k) = 0.d5 * sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)

            end do
         end do
      end do

   end subroutine compute_vort

   !-----------------------------------------------------------------------!

end module compute_vort_module
