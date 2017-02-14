MODULE CORRECT_0_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CORRECT_0                                               C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
!                                                                      C
!  Purpose: Correct the fluid pressure and gas velocities              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CORRECT_0(slo,shi,ulo,uhi,vlo,vhi,wlo,whi,lo,hi,&
           p_g,pp_g,u_g,v_g,w_g,d_e,d_n,d_t)&
           bind(C, name="correct_0")

      USE ur_facs  , only: ur_fac

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

      real(c_real), intent(inout   ) :: pp_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: d_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: d_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: d_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I,J,K
      integer :: funit
!-----------------------------------------------

! Underrelax pressure correction.  Velocity corrections should not be
! underrelaxed, so that the continuity eq. is satisfied.




      if(vhi(3) == 1) then ! FLD01-X
         k = 0
         write(2333,"(2/,'Pass=',i2)") 0
         do j = vhi(2)-1,vlo(2)+1,-1
            write(2333,"(i3,1x)",advance='no') j
            do i = vlo(1),vhi(1)-1
               write(2333,"(1x,es10.2)",advance='no') v_g(i,j,k)
            end do
            i=vhi(1)
            write(2333,"(1x,es10.2)",advance='yes') v_g(i,j,k)
         end do

         k = 0
         write(2334,"(2/,'Pass=',i4)") 2334
         do j = vhi(2)-1,vlo(2)+1,-1
            write(2334,"(i3,1x)",advance='no') j
            do i = vlo(1),vhi(1)-1
               write(2334,"(1x,es10.2)",advance='no') d_n(i,j,k)
            end do
            i=vhi(1)
            write(2334,"(1x,es10.2)",advance='yes') d_n(i,j,k)
         end do

         k = 0
         write(2335,"(2/,'Pass=',i4)") 2335
         do j = vhi(2)-1,vlo(2)+1,-1
            write(2335,"(i3,1x)",advance='no') j
            do i = vlo(1),vhi(1)-1
               write(2335,"(1x,es10.2)",advance='no') (pp_g(i,j+1,k)-pp_g(i,j,k))
            end do
            i=vhi(1)
            write(2335,"(1x,es10.2)",advance='yes') (pp_g(i,j+1,k)-pp_g(i,j,k))
         end do

         funit = 9100! + count
         k = 0
         write(funit,"(2/,'Pass=',i2)") k
         do i = slo(1),shi(1)-1
            write(funit,"(5x,i2,4x)",advance='no')i
         enddo
         write(funit,"(5x,i2,4x)",advance='yes')i
         i=shi(1)
         do j = shi(2),slo(2),-1
            do i = slo(1),shi(1)-1
               write(funit,"(1x,es10.2)",advance='no') pp_g(i,j,k)
            end do
            i=shi(1)
            write(funit,"(1x,es10.2)",advance='yes') pp_g(i,j,k)
         end do

         k = 0
         write(2330,"(2/,'Pass=',i4)") 2330
         do j = uhi(2)-1,ulo(2)+1,-1
            write(2330,"(i3,1x)",advance='no') j
            do i = ulo(1),uhi(1)-1
               write(2330,"(1x,es10.2)",advance='no') u_g(i,j,k)
            end do
            i=uhi(1)
            write(2330,"(1x,es10.2)",advance='yes') u_g(i,j,k)
         end do

         k = 0
         write(2331,"(2/,'Pass=',i4)") 2331
         do j = uhi(2)-1,ulo(2)+1,-1
            write(2331,"(i3,1x)",advance='no') j
            do i = ulo(1),uhi(1)-1
               write(2331,"(1x,es10.2)",advance='no') d_e(i,j,k)
            end do
            i=uhi(1)
            write(2331,"(1x,es10.2)",advance='yes') d_e(i,j,k)
         end do

      else if(vhi(1) == 1) then  ! FLD01-Y
         i = 0
         write(2333,"(2/,'Pass=',i2)") 0
         do k = vhi(3)-1,vlo(3)+1,-1
            write(2333,"(i3,1x)",advance='no') j
            do j = vlo(2),vhi(2)-1
               write(2333,"(1x,es10.2)",advance='no') v_g(i,j,k)
            end do
            j=vhi(2)
            write(2333,"(1x,es10.2)",advance='yes') v_g(i,j,k)
         end do

         i = 0
         write(2433,"(2/,'Pass=',i2)") 0
         do k = whi(3)-1,wlo(3)+1,-1
            write(2433,"(i3,1x)",advance='no') j
            do j = wlo(2),whi(2)-1
               write(2433,"(1x,es10.2)",advance='no') w_g(i,j,k)
            end do
            j=whi(2)
            write(2433,"(1x,es10.2)",advance='yes') w_g(i,j,k)
         end do

      endif











      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)
              p_g(i,j,k) = p_g(i,j,k) + ur_fac(1)*pp_g(i,j,k)
          enddo
        enddo
      enddo

      do k = ulo(3),uhi(3)
        do j = ulo(2),uhi(2)
          do i = ulo(1)+1,uhi(1)-1
            u_g(i,j,k) = u_g(i,j,k) - d_e(i,j,k)*(pp_g(i+1,j,k)-pp_g(i,j,k))
          enddo
        enddo
     enddo

      do k = vlo(3),  vhi(3)
        do j = vlo(2)+1, vhi(2)-1
          do i = vlo(1),  vhi(1)
            v_g(i,j,k) = v_g(i,j,k) - d_n(i,j,k)*(pp_g(i,j+1,k)-pp_g(i,j,k))
          enddo
        enddo
      enddo

      do k = wlo(3)+1,whi(3)-1
        do j =wlo(2),  whi(2)
          do i =wlo(1),  whi(1)
            w_g(i,j,k) = w_g(i,j,k) - d_t(i,j,k)*(pp_g(i,j,k+1) - pp_g(i,j,k))
          enddo
        enddo
     enddo


      if(vhi(3) == 1) then
         k = 0
         write(2333,"(2/,'Pass=',i2)") 1
         do j = vhi(2)-1,vlo(2)+1,-1
            write(2333,"(i3,1x)",advance='no') j
            do i = vlo(1),vhi(1)-1
               write(2333,"(1x,es10.2)",advance='no') v_g(i,j,k)
            end do
            i=vhi(1)
            write(2333,"(1x,es10.2)",advance='yes') v_g(i,j,k)
         end do

      else if(vhi(1) == 1) then
         i = 0
         write(2333,"(2/,'Pass=',i2)") 1
         do k = vhi(3)-1,vlo(3)+1,-1
            write(2333,"(i3,1x)",advance='no') j
            do j = vlo(2),vhi(2)-1
               write(2333,"(1x,es10.2)",advance='no') v_g(i,j,k)
            end do
            j=vhi(2)
            write(2333,"(1x,es10.2)",advance='yes') v_g(i,j,k)
         end do
      endif


      end subroutine correct_0

end module correct_0_module
