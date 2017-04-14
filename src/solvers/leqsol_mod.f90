module leqsol

   use compar, only: mype
   use error_manager, only: ival, flush_err_msg, err_msg
   use exit_mod, only: mfix_exit
   use param, only: DIM_EQS
   use param, only: zero

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

! Maximum number of outer iterations
  integer :: MAX_NIT

! Automatic adjustment of leq parameters possible (set in iterate after
! the completion of first iteration).
  logical :: LEQ_ADJUST

! Maximum number of linear equation solver iterations
  integer :: LEQ_IT(DIM_EQS)

! Total Iterations
  integer :: ITER_TOT(DIM_EQS) = 0

! Linear equation solver sweep direction
  CHARACTER(LEN=4) :: LEQ_SWEEP(DIM_EQS)

! Linear equation solver tolerance
  real(c_real) :: LEQ_TOL(DIM_EQS)

! Preconditioner option
  CHARACTER(LEN=4) :: LEQ_PC(DIM_EQS)

! Option to minimize dot products
  logical :: MINIMIZE_DOTPRODUCTS

! Option to transpose A_m
  logical :: DO_TRANSPOSE

! Frequency of convergence check in BiCGStab
  integer :: ICHECK_BICGS

! Optimize for massively parallel machine
  logical :: OPT_PARALLEL

! Linear and non-linear solver statistics
  logical :: SOLVER_STATISTICS

  logical :: USE_DOLOOP

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MSOLVE0                                             C
!  Notes: do nothing or no preconditioning (leq_pc='none')             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE LEQ_MSOLVE0(slo, shi, B_m, Var) &
     bind(C, name = "leq_msolve0")

    IMPLICIT NONE

  integer(c_int), intent(in) :: slo(3),shi(3)

! Vector b_m
    real(c_real), INTENT(IN) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
! Variable
    real(c_real), INTENT(OUT) :: Var&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    integer :: i,j,k
!-----------------------------------------------

! do nothing or no preconditioning
    do k = slo(3),shi(3)
       do j = slo(2),shi(2)
          do i = slo(1),shi(1)
             var(i,j,k) = b_m(i,j,k)
          enddo
       enddo
    enddo

  end subroutine leq_msolve0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MSOLVE1                                             C
!  Notes: diagonal scaling (leq_pc='diag')                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  subroutine leq_msolve1(b_m, blo, bhi, A_m, alo, ahi, var, vlo, vhi) &
     bind(C, name = "leq_msolve1")
   use param, only: small_number

    IMPLICIT NONE

    integer(c_int), intent(in) :: blo(3),bhi(3),alo(3),ahi(3),vlo(3),vhi(3)

    ! Vector b_m
    real(c_real), INTENT(IN) :: b_m&
         (blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))

    ! Septadiagonal matrix A_m
    real(c_real), INTENT(IN) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

    real(c_real), INTENT(OUT) :: var&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

    integer :: i,j,k

    do k = alo(3),ahi(3)
       do j = alo(2),ahi(2)
          do i = alo(1),ahi(1)

             if(abs(A_m(i,j,k,0))>small_number) &
                var(i,j,k) = b_m(i,j,k)/A_m(i,j,k,0)

          enddo
       enddo
    enddo

  end subroutine leq_msolve1

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  real(c_real) function dot_product_par(r1,r2,slo,shi)

    use amrex_fort_module, only : c_real => amrex_real

    implicit none
    integer, intent(in) :: slo(3),shi(3)

    real(c_real), intent(in) :: r1&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(c_real), intent(in) :: r2&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(c_real) :: prod
    integer :: i, j, k

    prod = 0.0d0

    do k = slo(3),shi(3)
       do i = slo(1),shi(1)
          do j = slo(2),shi(2)
               prod = prod + r1(i,j,k)*r2(i,j,k)
          enddo
       enddo
    enddo

    dot_product_par = prod

  end function dot_product_par

end module leqsol
