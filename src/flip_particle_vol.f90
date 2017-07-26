!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: flip_particle_vol                                      C
!                                                                      C
!  Purpose: This subroutine does the initial setting of all boundary   C
!  conditions. The user specifications of the boundary conditions are  C
!  checked for veracity in various check_data/ routines:               C
!  (e.g., check_boundary_conditions).                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine flip_particle_vol(slo, shi, ep_g, &
                             bc_ilo_type, bc_ihi_type, &
                             bc_jlo_type, bc_jhi_type, &
                             bc_klo_type, bc_khi_type, domlo, domhi) &
     bind(C, name="flip_particle_vol")

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  use bc, only: nsw_, fsw_, psw_, pinf_, pout_, minf_

  implicit none

  integer(c_int), intent(in   ) :: slo(3),shi(3)
  integer(c_int), intent(in   ) :: domlo(3),domhi(3)

  real(c_real), intent(inout) :: &
       ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

  integer(c_int), intent(in   ) :: &
       bc_ilo_type(domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2), &
       bc_ihi_type(domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2), &
       bc_jlo_type(domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2), &
       bc_jhi_type(domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2), &
       bc_klo_type(domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2), &
       bc_khi_type(domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)

! Local variables
!--------------------------------------------------------------------//
  integer(c_int) :: i,j,k
  integer(c_int) :: ilo, ihi, jlo, jhi, klo, khi
  real(c_real)   :: vol_outside

!--------------------------------------------------------------------//


  if (slo(1).lt.domlo(1)) then
     ilo = domlo(1)
     do k=slo(3),shi(3)
        do j=slo(2),shi(2)
           if (bc_ilo_type(j,k,1) == NSW_ .or. &
               bc_ilo_type(j,k,1) == PSW_ .or. &
               bc_ilo_type(j,k,1) == FSW_ .or. &
               bc_ilo_type(j,k,1) == PINF_ .or. &
               bc_ilo_type(j,k,1) == MINF_ .or. &
               bc_ilo_type(j,k,1) == POUT_) then
 
               vol_outside = 1.0 - ep_g(ilo-1,j,k)
               ep_g(ilo,j,k) = ep_g(ilo,j,k) - vol_outside

           end if
        end do
     end do
  endif

  if (shi(1).gt.domhi(1)) then
     ihi = domhi(1)
     do k=slo(3),shi(3)
        do j=slo(2),shi(2)
           if (bc_ihi_type(j,k,1) == NSW_ .or. &
               bc_ihi_type(j,k,1) == PSW_ .or. &
               bc_ihi_type(j,k,1) == FSW_ .or. &
               bc_ihi_type(j,k,1) == PINF_ .or. &
               bc_ihi_type(j,k,1) == MINF_ .or. &
               bc_ihi_type(j,k,1) == POUT_) then

               vol_outside = 1.0 - ep_g(ihi+1,j,k)
               ep_g(ihi,j,k) = ep_g(ihi,j,k) - vol_outside

           end if
        end do
     end do
  endif

  if (slo(2).lt.domlo(2)) then
     jlo = domlo(2)
     do k=slo(3),shi(3)
        do i=slo(1),shi(1)
           if (bc_jlo_type(i,k,1) == NSW_ .or. &
               bc_jlo_type(i,k,1) == PSW_ .or. &
               bc_jlo_type(i,k,1) == FSW_ .or. &
               bc_jlo_type(j,k,1) == PINF_ .or. &
               bc_jlo_type(j,k,1) == MINF_ .or. &
               bc_jlo_type(j,k,1) == POUT_) then

               vol_outside = 1.0 - ep_g(i,jlo-1,k)
               ep_g(i,jlo,k) = ep_g(i,jlo,k) - vol_outside

           end if
        end do
     end do
  endif

  if (shi(2).gt.domhi(2)) then
     jhi = domhi(2)
     do k=slo(3),shi(3)
        do i=slo(1),shi(1)
           if (bc_jhi_type(i,k,1) == NSW_ .or. &
               bc_jhi_type(i,k,1) == PSW_ .or. &
               bc_jhi_type(i,k,1) == FSW_ .or. &
               bc_jhi_type(j,k,1) == PINF_ .or. &
               bc_jhi_type(j,k,1) == MINF_ .or. &
               bc_jhi_type(j,k,1) == POUT_) then

               vol_outside = 1.0 - ep_g(i,jhi+1,k)
               ep_g(i,jhi,k) = ep_g(i,jhi,k) - vol_outside

           end if
        end do
     end do
  endif

  if (slo(3).lt.domlo(3)) then
     klo = domlo(3)
     do j=slo(2),shi(2)
        do i=slo(1),shi(1)
           if (bc_klo_type(i,j,1) == NSW_ .or. &
               bc_klo_type(i,j,1) == PSW_ .or. &
               bc_klo_type(i,j,1) == FSW_) then

               vol_outside = 1.0 - ep_g(i,j,klo-1)
               ep_g(i,j,klo) = ep_g(i,j,klo) - vol_outside

           end if
        end do
     end do
  endif

  if (shi(3).gt.domhi(3)) then
     khi = domhi(3)
     do j=slo(2),shi(2)
        do i=slo(1),shi(1)
           if (bc_khi_type(i,j,1) == NSW_ .or. &
               bc_khi_type(i,j,1) == PSW_ .or. &
               bc_khi_type(i,j,1) == FSW_) then

               vol_outside = 1.0 - ep_g(i,j,khi+1)
               ep_g(i,j,khi) = ep_g(i,j,khi) - vol_outside

           end if
        end do
     end do
  endif

end subroutine flip_particle_vol
