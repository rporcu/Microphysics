MODULE SOLVE_LIN_EQ_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_LIN_EQ                                            C
!  Purpose: Interface for linear equation solver                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOLVE_LIN_EQ(Vno, VAR, A_M, B_M)&
         bind(C, name="mfix_solve_lin_eq")

         use compar, only: istart3, iend3
         use compar, only: jstart3, jend3
         use compar, only: kstart3, kend3
         USE leqsol  , only: leq_it, leq_sweep, leq_tol, leq_pc
         use bl_fort_module, only : c_real
         use iso_c_binding , only: c_int

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable number
!     1 = pressure correction equation
!     3 = gas u-momentum
!     4 = gas v-momentum
!     5 = gas w-momentum
      integer(c_int), intent(in   ) :: vno
! variable
      real(c_real), intent(inout) :: var&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! septadiagonal matrix a_m
      real(c_real), intent(inout) :: a_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! vector b_m
      real(c_real), intent(inout) :: b_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer :: ier

      call leq_bicgs('---', vno, var, a_m, b_m, leq_sweep(vno), &
         leq_tol(vno), leq_pc(vno), leq_it(vno), ier)


      RETURN
      END SUBROUTINE SOLVE_LIN_EQ
END MODULE SOLVE_LIN_EQ_MODULE
