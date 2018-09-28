module correct_0_module

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: correct_p_0/correct_u_0/correct_v_0/correct_w_0         C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
!                                                                      C
!  Purpose: Correct the fluid pressure and gas velocities              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine correct_p_0(lo,hi,slo,shi,p_g,pp_g) &
           bind(C, name="correct_p_0")

      use ur_facs  , only: ur_fac

      implicit none

      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: slo(3),shi(3)

      real(rt), intent(in   ) :: pp_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int) :: i,j,k

      ! Underrelax pressure correction.  Velocity corrections should not be
      ! underrelaxed, so that the continuity eq. is satisfied.
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               p_g(i,j,k) = p_g(i,j,k) + ur_fac(1)*pp_g(i,j,k)
            enddo
         enddo
      enddo

      end subroutine correct_p_0

      subroutine correct_u_0(lo,hi,ulo,uhi,slo,shi,pp_g,u_g,d_e)&
           bind(C, name="correct_u_0")

      implicit none

      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3)

      real(rt), intent(in   ) :: pp_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(in   ) :: d_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(rt), intent(inout) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      integer(c_int) :: i,j,k

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            u_g(i,j,k) = u_g(i,j,k) - d_e(i,j,k)*(pp_g(i,j,k)-pp_g(i-1,j,k))
          enddo
        enddo
     enddo

     end subroutine correct_u_0

     subroutine correct_v_0(lo,hi,vlo,vhi,slo,shi,pp_g,v_g,d_n)&
           bind(C, name="correct_v_0")

      implicit none

      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: vlo(3),vhi(3)

      real(rt), intent(in   ) :: pp_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(in   ) :: d_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(rt), intent(inout) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      integer(c_int) :: i,j,k

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
           v_g(i,j,k) = v_g(i,j,k) - d_n(i,j,k)*(pp_g(i,j,k)-pp_g(i,j-1,k))
         enddo
       enddo
     enddo

     end subroutine correct_v_0

     subroutine correct_w_0(lo,hi,wlo,whi,slo,shi,pp_g,w_g,d_t)&
           bind(C, name="correct_w_0")

      implicit none

      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: wlo(3),whi(3)

      real(rt), intent(in   ) :: pp_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(in   ) :: d_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(rt), intent(inout) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer(c_int) :: i,j,k

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
           w_g(i,j,k) = w_g(i,j,k) - d_t(i,j,k)*(pp_g(i,j,k)-pp_g(i,j,k-1))
         enddo
       enddo
     enddo

     end subroutine correct_w_0

end module correct_0_module
