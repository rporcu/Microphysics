module matvec_module

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   implicit none

contains

   subroutine leq_scale(lo, hi, rhs, rlo, rhi, A_m, alo, ahi) &
      bind(C, name = "leq_scale")
   use param, only: small_number, one

   integer(c_int), intent(in   ) ::  lo(3), hi(3)
   integer(c_int), intent(in   ) :: rlo(3),rhi(3)
   integer(c_int), intent(in   ) :: alo(3),ahi(3)

   real(rt)  , intent(inout) :: rhs&
      (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
   real(rt)  , intent(inout) :: A_m&
      (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

   integer :: i,j,k
   real(rt) :: aijmax, oam

   do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            aijmax = maxval(abs(a_m(i,j,k,:)) )
            if (aijmax > small_number) then
               oam = one/aijmax
               a_m(i,j,k,:) = a_m(i,j,k,:)*oam
               rhs(i,j,k)   = rhs(i,j,k  )*oam
            endif
         enddo
      enddo
   enddo

  end subroutine leq_scale

  ! returns Am*Var
  subroutine leq_matvec(lo, hi, var, vlo, vhi, A_m, alo, ahi, res, rlo, rhi)&
    bind(C, name = "leq_matvec")

    integer(c_int), intent(in   ) ::  lo(3), hi(3)
    integer(c_int), intent(in   ) :: vlo(3),vhi(3)
    integer(c_int), intent(in   ) :: alo(3),ahi(3)
    integer(c_int), intent(in   ) :: rlo(3),rhi(3)

    real(rt)  , intent(in   ) :: var&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

    real(rt)  , intent(in   ) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

    real(rt)  , intent(  out) :: res&
         (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Variable
    integer :: I, J, K
!-----------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
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
  subroutine leq_residual(lo, hi, rhs, hlo, hhi, var, vlo, vhi, A_m, alo, ahi, res, rlo, rhi) &
    bind(C, name = "leq_residual")

    integer(c_int), intent(in   ) ::  lo(3), hi(3)
    integer(c_int), intent(in   ) :: hlo(3),hhi(3),vlo(3),vhi(3)
    integer(c_int), intent(in   ) :: alo(3),ahi(3),rlo(3),rhi(3)

    real(rt)  , intent(in   ) :: rhs&
         (hlo(1):hhi(1),hlo(2):hhi(2),hlo(3):hhi(3))

    real(rt)  , intent(in   ) :: var&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

    real(rt)  , intent(in   ) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

    real(rt)  , intent(  out) :: res&
         (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))

    integer :: I, J, K

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1),hi(1)
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






  subroutine out_array(array, lo, hi) &
    bind(C, name = "out_array")

    integer(c_int), intent(in   ) :: lo(3),hi(3)

    real(rt)  , intent(in   ) :: array&
         (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    integer :: i, j, k

    do j = lo(2),hi(2)
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             write(*,*) i,j,k,array(i,j,k)
          enddo
       enddo
    enddo

  end subroutine out_array


  subroutine out_matrix(matrix, lo, hi) &
       bind(C, name = "out_matrix")

    integer(c_int), intent(in   ) :: lo(3),hi(3)

    real(rt)  , intent(in   ) :: matrix&
         (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),-3:3)

    integer :: i, j, k

    do j = lo(2),hi(2)
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             write(*,*) i,j,k,matrix(i,j,k,:)
          enddo
       enddo
    enddo

  end subroutine out_matrix






end module matvec_module
