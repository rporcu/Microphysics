module matvec_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   implicit none

contains

   subroutine leq_scale(rhs, rlo, rhi, A_m, alo, ahi) &
      bind(C, name = "leq_scale")
   use param1, only: small_number, one

   integer(c_int), intent(in   ) :: rlo(3),rhi(3),alo(3),ahi(3)

   real(c_real)  , intent(inout) :: rhs&
      (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
   real(c_real)  , intent(inout) :: A_m&
      (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

   integer :: i,j,k
   real(c_real) :: aijmax, oam

   do k = alo(3),ahi(3)
      do j = alo(2),ahi(2)
         do i = alo(1),ahi(1)
            aijmax = maxval(abs(a_m(i,j,k,:)) )
            if(aijmax > small_number) then
               oam = one/aijmax
               a_m(i,j,k,:) = a_m(i,j,k,:)*oam
               rhs(i,j,k)   = rhs(i,j,k)*oam
            endif
         enddo
      enddo
   enddo

  end subroutine leq_scale

  ! returns Am*Var
  subroutine leq_matvec(var, vlo, vhi, A_m, alo, ahi, res, rlo, rhi)&
    bind(C, name = "leq_matvec")

    integer(c_int), intent(in   ) :: vlo(3),vhi(3)
    integer(c_int), intent(in   ) :: alo(3),ahi(3)
    integer(c_int), intent(in   ) :: rlo(3),rhi(3)

    real(c_real)  , intent(in   ) :: var&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

    real(c_real)  , intent(in   ) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

    real(c_real)  , intent(  out) :: res&
         (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Variable
    integer :: I, J, K
!-----------------------------------------------

    do k = alo(3),ahi(3)
       do j = alo(2),ahi(2)
          do i = alo(1),ahi(1)
             res(i,j,k) = &
                          ( A_m(i,j,k,-3) * Var(i,j,k-1)    &
                          + A_m(i,j,k,-2) * Var(i,j-1,k)    &
                          + A_m(i,j,k,-1) * Var(i-1,j,k)    &
                          + A_m(i,j,k, 0) * Var(i,j,k)      &
                          + A_m(i,j,k, 1) * Var(i+1,j,k)    &
                          + A_m(i,j,k, 2) * Var(i,j+1,k)    &
                          + A_m(i,j,k, 3) * Var(i,j,k+1)  )
               enddo
            enddo
         enddo

  end subroutine leq_matvec

! returns (Rhs - Am*Var)
  subroutine leq_residual(rhs, hlo, hhi, var, vlo, vhi, A_m, alo, ahi, res, rlo, rhi) &
    bind(C, name = "leq_residual")

    integer(c_int), intent(in   ) :: hlo(3),hhi(3),vlo(3),vhi(3)
    integer(c_int), intent(in   ) :: alo(3),ahi(3),rlo(3),rhi(3)

    real(c_real)  , intent(in   ) :: rhs&
         (hlo(1):hhi(1),hlo(2):hhi(2),hlo(3):hhi(3))

    real(c_real)  , intent(in   ) :: var&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

    real(c_real)  , intent(in   ) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

    real(c_real)  , intent(  out) :: res&
         (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))

    integer :: I, J, K

    do k = alo(3),ahi(3)
       do j = alo(2),ahi(2)
          do i = alo(1),ahi(1)
             res(i,j,k) = rhs(i,j,k) -  &
                ( A_m(i,j,k,-3) * Var(i,j,k-1)    &
                + A_m(i,j,k,-2) * Var(i,j-1,k)    &
                + A_m(i,j,k,-1) * Var(i-1,j,k)    &
                + A_m(i,j,k, 0) * Var(i,j,k)      &
                + A_m(i,j,k, 1) * Var(i+1,j,k)    &
                + A_m(i,j,k, 2) * Var(i,j+1,k)    &
                + A_m(i,j,k, 3) * Var(i,j,k+1)  )
          enddo
       enddo
    enddo

  end subroutine leq_residual

end module matvec_module
