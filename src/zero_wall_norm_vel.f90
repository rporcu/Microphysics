!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: zero_wall_norm_vel                                      C
!                                                                      C
!  Purpose: This subroutine does the initial setting of all boundary   C
!  conditions. The user specifications of the boundary conditions are  C
!  checked for veracity in various check_data/ routines:               C
!  (e.g., check_boundary_conditions).                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine zero_wall_norm_vel(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
     u_g, v_g, w_g, bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
     bc_klo_type, bc_khi_type, domlo, domhi) &
     bind(C, name="zero_wall_norm_vel")

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  use bc, only: NSW_, FSW_, PSW_

  implicit none

  integer(c_int), intent(in   ) :: slo(3),shi(3)
  integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
  integer(c_int), intent(in   ) :: domlo(3),domhi(3)

  real(c_real), intent(inout) :: &
       u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
       v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
       w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

  integer(c_int), intent(in   ) :: &
       bc_ilo_type(domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2), &
       bc_ihi_type(domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2), &
       bc_jlo_type(domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2), &
       bc_jhi_type(domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2), &
       bc_klo_type(domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2), &
       bc_khi_type(domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)

! Local variables
!--------------------------------------------------------------------//
! local index for boundary condition
  integer :: bcv, i,j,k

  integer    nlft, nrgt, nbot, ntop, nup, ndwn
  integer    ilo, ihi, jlo, jhi, klo, khi

!--------------------------------------------------------------------//

  nlft = max(0,domlo(1)-slo(1))
  nbot = max(0,domlo(2)-slo(2))
  ndwn = max(0,domlo(3)-slo(3))

  nrgt = max(0,shi(1)-domhi(1))
  ntop = max(0,shi(2)-domhi(2))
  nup  = max(0,shi(3)-domhi(3))

  if (nlft .gt. 0) then
     ilo = domlo(1)
     do i = 0, nlft
        do k=slo(3),shi(3)
           do j=slo(2),shi(2)
              bcv = bc_ilo_type(j,k,2)
              if(bc_ilo_type(j,k,1) == NSW_ .or. &
                 bc_ilo_type(j,k,1) == PSW_ .or. &
                 bc_ilo_type(j,k,1) == FSW_) then

                 u_g(ilo-i,j,k) = 0.0d0

              end if
           end do
        end do
     end do
  endif

  if (nrgt .gt. 0) then
     ihi = domhi(1)
     do i = 0, nrgt
        do k=slo(3),shi(3)
           do j=slo(2),shi(2)
              bcv = bc_ihi_type(j,k,2)
              if(bc_ihi_type(j,k,1) == NSW_ .or. &
                 bc_ihi_type(j,k,1) == PSW_ .or. &
                 bc_ihi_type(j,k,1) == FSW_) then

                 u_g(ihi+1+i,j,k) = 0.0d0

              end if
           end do
        end do
     end do
  endif

  if (nbot .gt. 0) then
     jlo = domlo(2)
     do j = 0, nbot
        do k=slo(3),shi(3)
           do i=slo(1),shi(1)
              bcv = bc_jlo_type(i,k,2)
              if(bc_jlo_type(i,k,1) == NSW_ .or. &
                 bc_jlo_type(i,k,1) == PSW_ .or. &
                 bc_jlo_type(i,k,1) == FSW_) then

                 v_g(i,jlo-j,k) = 0.0d0

              end if
           end do
        end do
     end do
  endif

  if (ntop .gt. 0) then
     jhi = domhi(2)
     do j = 0, ntop
        do k=slo(3),shi(3)
           do i=slo(1),shi(1)
              bcv = bc_jhi_type(i,k,2)
              if(bc_jhi_type(i,k,1) == NSW_ .or. &
                 bc_jhi_type(i,k,1) == PSW_ .or. &
                 bc_jhi_type(i,k,1) == FSW_) then

                 v_g(i,jhi+1+j,k) = 0.0d0

              end if
           end do
        end do
     end do
  endif

  if (ndwn .gt. 0) then
     klo = domlo(3)
     do k = 0, ndwn
        do j=slo(2),shi(2)
           do i=slo(1),shi(1)
              bcv = bc_klo_type(i,j,2)
              if(bc_klo_type(i,j,1) == NSW_ .or. &
                 bc_klo_type(i,j,1) == PSW_ .or. &
                 bc_klo_type(i,j,1) == FSW_) then

                 w_g(i,j,klo-k) = 0.0d0

              end if
           end do
        end do
     end do
  endif

  if (nup .gt. 0) then
     khi = domhi(3)
     do k = 0, nup
        do j=slo(2),shi(2)
           do i=slo(1),shi(1)
              bcv = bc_khi_type(i,j,2)
              if(bc_khi_type(i,j,1) == NSW_ .or. &
                 bc_khi_type(i,j,1) == PSW_ .or. &
                 bc_khi_type(i,j,1) == FSW_) then

                 w_g(i,j,khi+1+k) = 0.0d0

              end if
           end do
        end do
     end do
  endif

end subroutine zero_wall_norm_vel
