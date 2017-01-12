module matvec_module

contains



   subroutine leq_scale(rhs, A_m, slo, shi, lo, hi)&
      bind(C, name = "leq_scale")

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int
   use param1, only: small_number, one

   implicit none
   integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

   real(c_real)  , intent(inout) :: rhs&
      (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
   real(c_real)  , intent(inout) :: A_m&
      (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

   integer :: i,j,k
   real(c_real) :: aijmax, oam

   do k = lo(3)-1, hi(3)+1
      do i = lo(1)-1, hi(1)+1
         do j = lo(2)-1, hi(2)+1
            aijmax = maxval(abs(a_m(i,j,k,:)) )
            if(aijmax > small_number) then
               oam = one/aijmax
               a_m(i,j,k,:) = a_m(i,j,k,:)*oam
               rhs(i,j,k)   = rhs(i,j,k)*oam
            endif
         enddo
      enddo
   enddo

   1000 format(a1,8(2x,es12.4))
end subroutine leq_scale


  ! returns Am*Var
  subroutine leq_matvec(var, A_m, res, slo, shi, lo, hi)&
    bind(C, name = "leq_matvec")

    use bl_fort_module, only : c_real
    use iso_c_binding , only: c_int
    use functions, only: iminus,iplus,jminus,jplus,kminus,kplus

    IMPLICIT NONE

    integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

    real(c_real)  , intent(in   ) :: var&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(c_real)  , intent(in   ) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

    real(c_real)  , intent(  out) :: res&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Variable
    INTEGER :: I, J, K
!-----------------------------------------------

    do k = slo(3),shi(3)
       do j = slo(2),shi(2)
          do i = slo(1),shi(1)
             res(i,j,k) = &
                          ( A_m(i,j,k,-3) * Var(i,j,kminus(i,j,k))   &
                          + A_m(i,j,k,-2) * Var(i,jminus(i,j,k),k)   &
                          + A_m(i,j,k,-1) * Var(iminus(i,j,k),j,k)   &
                          + A_m(i,j,k, 0) * Var(i,j,k)               &
                          + A_m(i,j,k, 1) * Var(iplus(i,j,k),j,k)    &
                          + A_m(i,j,k, 2) * Var(i,jplus(i,j,k),k)    &
                          + A_m(i,j,k, 3) * Var(i,j,kplus(i,j,k))  )
               enddo
            enddo
         enddo

  end subroutine leq_matvec





! returns (Rhs - Am*Var)
  subroutine leq_residual(rhs, var, A_m, res, slo, shi, lo, hi) &
    bind(C, name = "leq_residual")

    use bl_fort_module, only : c_real
    use iso_c_binding , only: c_int
    use functions, only: iminus,iplus,jminus,jplus,kminus,kplus

    IMPLICIT NONE

    integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

    real(c_real)  , intent(in   ) :: rhs&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(c_real)  , intent(in   ) :: var&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(c_real)  , intent(in   ) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

    real(c_real)  , intent(  out) :: res&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Variable
    INTEGER :: I, J, K
!-----------------------------------------------

    do k = slo(3),shi(3)
       do i = slo(1),shi(1)
          do j = slo(2),shi(2)
             res(i,j,k) = rhs(i,j,k) -  &
                ( A_m(i,j,k,-3) * Var(i,j,kminus(i,j,k))   &
                + A_m(i,j,k,-2) * Var(i,jminus(i,j,k),k)   &
                + A_m(i,j,k,-1) * Var(iminus(i,j,k),j,k)   &
                + A_m(i,j,k, 0) * Var(i,j,k)               &
                + A_m(i,j,k, 1) * Var(iplus(i,j,k),j,k)    &
                + A_m(i,j,k, 2) * Var(i,jplus(i,j,k),k)    &
                + A_m(i,j,k, 3) * Var(i,j,kplus(i,j,k))  )

!             write(6,1000) i,j,k,res(i,j,k)
          enddo
       enddo
    enddo

  end subroutine leq_residual

end module matvec_module
