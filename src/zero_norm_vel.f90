!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: zero_norm_vel                                           C
!  Purpose: Set the velocity component normal to a wall to zero        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

   subroutine zero_norm_vel(slo,shi, &
                            ulo,uhi,vlo,vhi,wlo,whi,&
                            u_g,v_g,w_g, &
                            bc_ilo_type, bc_ihi_type, &
                            bc_jlo_type, bc_jhi_type, &
                            bc_klo_type, bc_khi_type) &
         bind(C, name="zero_norm_vel")

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use param1  , only: zero
      use geometry, only: domlo,domhi
      use ic      , only: NSW_, FSW_, PSW_

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

      real(c_real), intent(inout) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer(c_int), intent(in   ) :: bc_ilo_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jlo_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_klo_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)
      integer(c_int), intent(in   ) :: bc_khi_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)

      integer :: i,j,k

      ! Lo i side
      if (slo(1) .lt. domlo(1)) then
        do k = slo(3),shi(3)
        do j = slo(2),shi(2)
           if (bc_ilo_type(j,k,1) == NSW_ .or. &
               bc_ilo_type(j,k,1) == FSW_ .or. &
               bc_ilo_type(j,k,1) == PSW_) then
               u_g(domlo(1)-1,j,k) = zero
           end if
        end do
        end do
      end if

      ! Hi i side
      if (shi(1) .gt. domhi(1)) then
        do k = slo(3),shi(3)
        do j = slo(2),shi(2)
           if (bc_ihi_type(j,k,1) == NSW_ .or. &
               bc_ihi_type(j,k,1) == FSW_ .or. &
               bc_ihi_type(j,k,1) == PSW_) then
               u_g(domhi(1),j,k) = zero
           end if
        end do
        end do
      end if

      ! Lo j side
      if (slo(2) .lt. domlo(2)) then
        do k = slo(3),shi(3)
        do i = slo(1),shi(1)
           if (bc_jlo_type(i,k,1) == NSW_ .or. &
               bc_jlo_type(i,k,1) == FSW_ .or. &
               bc_jlo_type(i,k,1) == PSW_) then
               v_g(i,domlo(2)-1,k) = zero
           end if
        end do
        end do
      end if

      ! Hi j side
      if (shi(2) .gt. domhi(2)) then
        do k = slo(3),shi(3)
        do i = slo(1),shi(1)
           if (bc_jhi_type(i,k,1) == NSW_ .or. &
               bc_jhi_type(i,k,1) == FSW_ .or. &
               bc_jhi_type(i,k,1) == PSW_) then
               v_g(i,domhi(2),k) = zero
           end if
        end do
        end do
      end if

      ! Lo k side
      if (slo(3) .lt. domlo(3)) then
        do j = slo(2),shi(2)
        do i = slo(1),shi(1)
           if (bc_klo_type(i,j,1) == NSW_ .or. &
               bc_klo_type(i,j,1) == FSW_ .or. &
               bc_klo_type(i,j,1) == PSW_) then
               w_g(i,j,domlo(3)-1) = zero
           end if
        end do
        end do
      end if

      ! Hi k side
      if (shi(3) .gt. domhi(3)) then
        do j = slo(2),shi(2)
        do i = slo(1),shi(1)
           if (bc_khi_type(i,j,1) == NSW_ .or. &
               bc_khi_type(i,j,1) == FSW_ .or. &
               bc_khi_type(i,j,1) == PSW_) then
               w_g(i,j,domhi(3)) = zero
           end if
        end do
        end do
      end if
 
   end subroutine zero_norm_vel
