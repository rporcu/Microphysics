module matvec_module

  contains 

  SUBROUTINE LEQ_MATVEC(rhs, var, A_m, res, slo, shi, lo, hi) &
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
    INTEGER :: I, J, K, IJK
!-----------------------------------------------

    do k = lo(3),hi(3)
       do i = lo(2),hi(2)
          do j = lo(1),hi(1)
             res(i,j,k) = rhs(i,j,k) -  & 
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

  END SUBROUTINE LEQ_MATVEC

end module matvec_module
