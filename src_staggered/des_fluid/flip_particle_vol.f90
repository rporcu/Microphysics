!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: flip_particle_vol                                       C
!                                                                      C
!  Purpose: This subroutine is used to maintain total particle volume  C
!           in the domain by "flipping" any volume deposited outside   C
!           the domain back into the domain.                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine flip_particle_vol(slo, shi, vol, &
                             bc_ilo_type, bc_ihi_type, &
                             bc_jlo_type, bc_jhi_type, &
                             bc_klo_type, bc_khi_type, domlo, domhi, ng) &
     bind(C, name="flip_particle_vol")

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

  use bc, only: nsw_, fsw_, psw_, pinf_, pout_, minf_

  implicit none

  integer(c_int), intent(in   ) :: slo(3),shi(3)
  integer(c_int), intent(in   ) :: domlo(3),domhi(3),ng

  real(rt), intent(inout) :: &
       vol(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

  integer(c_int), intent(in   ) :: &
       bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
       bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
       bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
       bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
       bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
       bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

! Local variables
!--------------------------------------------------------------------//
  integer(c_int) :: i,j,k
  integer(c_int) :: ilo, ihi, jlo, jhi, klo, khi

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

               vol(ilo,j,k) = vol(ilo,j,k) + vol(ilo-1,j,k)
               vol(ilo-1,j,k) = 0.d0

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

               vol(ihi,j,k) = vol(ihi,j,k) + vol(ihi+1,j,k)
               vol(ihi+1,j,k) = 0.d0

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
               bc_jlo_type(i,k,1) == PINF_ .or. &
               bc_jlo_type(i,k,1) == MINF_ .or. &
               bc_jlo_type(i,k,1) == POUT_) then

               vol(i,jlo,k) = vol(i,jlo,k) + vol(i,jlo-1,k)
               vol(i,jlo-1,k) = 0.d0

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
               bc_jhi_type(i,k,1) == PINF_ .or. &
               bc_jhi_type(i,k,1) == MINF_ .or. &
               bc_jhi_type(i,k,1) == POUT_) then

               vol(i,jhi,k) = vol(i,jhi,k) + vol(i,jhi+1,k)
               vol(i,jhi+1,k) = 0.d0

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
               bc_klo_type(i,j,1) == FSW_ .or. &
               bc_klo_type(i,j,1) == PINF_ .or. &
               bc_klo_type(i,j,1) == MINF_ .or. &
               bc_klo_type(i,j,1) == POUT_) then

               vol(i,j,klo) = vol(i,j,klo) + vol(i,j,klo-1)
               vol(i,j,klo-1) = 0.d0

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
               bc_khi_type(i,j,1) == FSW_ .or. &
               bc_khi_type(i,j,1) == PINF_ .or. &
               bc_khi_type(i,j,1) == MINF_ .or. &
               bc_khi_type(i,j,1) == POUT_) then

               vol(i,j,khi) = vol(i,j,khi) + vol(i,j,khi+1)
               vol(i,j,khi+1) = 0.d0

           end if
        end do
     end do
  endif

end subroutine flip_particle_vol

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: cap_eps                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine cap_eps(slo, shi, ep_g)  bind(C, name="mfix_cap_eps")

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

  use amrex_paralleldescriptor_module, only : amrex_pd_ioprocessor

  implicit none

  integer(c_int), intent(in   ) :: slo(3),shi(3)

  real(rt), intent(inout) :: &
       ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))


! Local variables
!--------------------------------------------------------------------//
  integer(c_int) :: i,j,k

! real(rt), parameter :: max_pack = 0.42_rt
  real(rt), parameter :: max_pack = 0.21_rt

!--------------------------------------------------------------------//

  return

  if(amrex_pd_ioprocessor()) write(*,*) 'WARNING: Applying maximum &
       &packing volume fraction: ', max_pack

  do k=slo(3),shi(3)
     do j=slo(2),shi(2)
        do i=slo(1),shi(1)
           ep_g(i,j,k) = max(max_pack, ep_g(i,j,k))
        end do
     end do
  enddo

end subroutine cap_eps
