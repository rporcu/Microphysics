!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RESET_NEW                                              C
!  Purpose: Reset the new variables with the stored previous-time-step C
!           values of field variables.                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE RESET_NEW


! Modules
!-----------------------------------------------
      USE fldvar, only: ep_g, ep_go
      USE fldvar, only: P_g, P_go
      USE fldvar, only: RO_g, RO_go
      USE fldvar, only: ROP_g, ROP_go

      USE fldvar, only: U_g, U_go
      USE fldvar, only: V_g, V_go
      USE fldvar, only: W_g, W_go

      use functions, only: funijk

      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3

      IMPLICIT NONE

! Local Variables
!-----------------------------------------------
      integer :: i, j, k, ijk, ier

      do k=kstart3, kend3
         do j=jstart3,jend3
            do i=istart3, iend3
               ijk = funijk(i,j,k)
               Ep_g(ijk) = Ep_go(i,j,k)
               P_g(ijk) = P_go(i,j,k)
               Ro_g(ijk) = Ro_go(i,j,k)
               Rop_g(i,j,k) = Rop_go(i,j,k)
               U_g(i,j,k) = U_go(i,j,k)
               V_g(i,j,k) = V_go(i,j,k)
               W_g(i,j,k) = W_go(i,j,k)
            enddo
         enddo
      enddo

! Recalculate all coefficients
      CALL CALC_COEFF_ALL (0, IER)

      RETURN
      END SUBROUTINE RESET_NEW
