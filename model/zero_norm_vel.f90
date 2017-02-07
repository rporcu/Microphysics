module zero_norm_vel_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
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

      subroutine zero_norm_vel(slo,shi,u_g,v_g,w_g,&
                               bc_ilo_type, bc_ihi_type, &
                               bc_jlo_type, bc_jhi_type, &
                               bc_klo_type, bc_khi_type)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1  , only: zero
      use geometry, only: domlo,domhi
      use ic      , only: NSW_, FSW_, PSW_

      IMPLICIT NONE

      integer(c_int), intent(in) :: slo(3), shi(3)

      real(c_real), intent(inout) ::  u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) ::  v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) ::  w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int), intent(in   ) :: bc_ilo_type&
         (slo(2):shi(2),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
         (slo(2):shi(2),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_jlo_type&
         (slo(1):shi(1),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
         (slo(1):shi(1),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_klo_type&
         (slo(1):shi(1),slo(2):shi(2),2)
      integer(c_int), intent(in   ) :: bc_khi_type&
         (slo(1):shi(1),slo(2):shi(2),2)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      integer :: i,j,k

      ! Lo i side
      if (slo(1) .lt. domlo(1)) then
        i = domlo(1)
        do k = slo(3),shi(3)
        do j = slo(2),shi(2)
           if (bc_ilo_type(j,k,1) == NSW_ .or. &
               bc_ilo_type(j,k,1) == FSW_ .or. &
               bc_ilo_type(j,k,1) == PSW_) then
               u_g(i-1,j,k) = zero
           end if
        end do
        end do
      end if

      ! Hi i side
      if (shi(1) .gt. domhi(1)) then
        i = domhi(1)
        do k = slo(3),shi(3)
        do j = slo(2),shi(2)
           if (bc_ihi_type(j,k,1) == NSW_ .or. &
               bc_ihi_type(j,k,1) == FSW_ .or. &
               bc_ihi_type(j,k,1) == PSW_) then
               u_g(i,j,k) = zero
           end if
        end do
        end do
      end if

      ! Lo j side
      if (slo(2) .lt. domlo(2)) then
        j = domlo(2)
        do k = slo(3),shi(3)
        do i = slo(1),shi(1)
           if (bc_jlo_type(i,k,1) == NSW_ .or. &
               bc_jlo_type(i,k,1) == FSW_ .or. &
               bc_jlo_type(i,k,1) == PSW_) then
               v_g(i,j-1,k) = zero
           end if
        end do
        end do
      end if

      ! Hi j side
      if (shi(2) .gt. domhi(2)) then
        j = domhi(2)
        do k = slo(3),shi(3)
        do i = slo(1),shi(1)
           if (bc_jhi_type(i,k,1) == NSW_ .or. &
               bc_jhi_type(i,k,1) == FSW_ .or. &
               bc_jhi_type(i,k,1) == PSW_) then
               v_g(i,j,k) = zero
           end if
        end do
        end do
      end if

      ! Lo k side
      if (slo(3) .lt. domlo(3)) then
        k = domlo(3)
        do j = slo(2),shi(2)
        do i = slo(1),shi(1)
           if (bc_klo_type(i,j,1) == NSW_ .or. &
               bc_klo_type(i,j,1) == FSW_ .or. &
               bc_klo_type(i,j,1) == PSW_) then
               w_g(i,j,k-1) = zero
           end if
        end do
        end do
      end if

      ! Hi k side
      if (shi(3) .gt. domhi(3)) then
        k = domhi(3)
        do j = slo(2),shi(2)
        do i = slo(1),shi(1)
           if (bc_khi_type(i,j,1) == NSW_ .or. &
               bc_khi_type(i,j,1) == FSW_ .or. &
               bc_khi_type(i,j,1) == PSW_) then
               w_g(i,j,k) = zero
           end if
        end do
        end do
      end if
 
      end subroutine zero_norm_vel

end module zero_norm_vel_module
