MODULE CORRECT_0_MODULE

   use amrex_fort_module, only : c_real => amrex_real
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
      SUBROUTINE CORRECT_0(slo,shi,ulo,uhi,vlo,vhi,wlo,whi,lo,hi,domlo,domhi,&
           p_g,pp_g,u_g,v_g,w_g,d_e,d_n,d_t)&
           bind(C, name="correct_0")

      use ur_facs  , only: ur_fac

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

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
      integer :: i,j,k
      integer :: llo(3), lhi(3)
      integer :: nlft, nrgt, nbot, ntop, nup, ndwn
!-----------------------------------------------

      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))

      llo(1) = slo(1) + nlft
      llo(2) = slo(2) + nbot
      llo(3) = slo(3) + ndwn

      lhi(1) = shi(1) - nrgt
      lhi(2) = shi(2) - ntop
      lhi(3) = shi(3) - nup

! Underrelax pressure correction.  Velocity corrections should not be
! underrelaxed, so that the continuity eq. is satisfied.
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               p_g(i,j,k) = p_g(i,j,k) + ur_fac(1)*pp_g(i,j,k)
            enddo
         enddo
      enddo


      do k = llo(3), lhi(3)
        do j = llo(2), lhi(2)
          do i = llo(1)+1, lhi(1)
            u_g(i,j,k) = u_g(i,j,k) - d_e(i,j,k)*(pp_g(i,j,k)-pp_g(i-1,j,k))
          enddo
        enddo
     enddo

     do k = llo(3), lhi(3)
       do j = llo(2)+1, lhi(2)
         do i = llo(1), lhi(1)
           v_g(i,j,k) = v_g(i,j,k) - d_n(i,j,k)*(pp_g(i,j,k)-pp_g(i,j-1,k))
         enddo
       enddo
     enddo

     do k = llo(3)+1, lhi(3)
       do j = llo(2), lhi(2)
         do i = llo(1), lhi(1)
           w_g(i,j,k) = w_g(i,j,k) - d_t(i,j,k)*(pp_g(i,j,k)-pp_g(i,j,k-1))
         enddo
       enddo
     enddo

     end subroutine correct_0

end module correct_0_module
