module leqsol

   use error_manager, only: ival, flush_err_msg, err_msg

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

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: leq_msolve1                                             C
!  Notes: diagonal scaling (leq_pc='diag')                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  subroutine leq_msolve1(lo, hi, b_m, blo, bhi, A_m, alo, ahi, var, vlo, vhi) &
     bind(C, name = "leq_msolve1")
   use param, only: small_number

    IMPLICIT NONE

    integer(c_int), intent(in) ::  lo(3), hi(3)
    integer(c_int), intent(in) :: blo(3),bhi(3)
    integer(c_int), intent(in) :: alo(3),ahi(3)
    integer(c_int), intent(in) :: vlo(3),vhi(3)

    ! Vector b_m
    real(c_real), intent(in) :: b_m&
         (blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))

    ! Septadiagonal matrix A_m
    real(c_real), intent(in) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

    real(c_real), intent(out) :: var&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

    integer :: i,j,k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(abs(A_m(i,j,k,0))>small_number) &
                var(i,j,k) = b_m(i,j,k)/A_m(i,j,k,0)

          enddo
       enddo
    enddo

  end subroutine leq_msolve1

end module leqsol
